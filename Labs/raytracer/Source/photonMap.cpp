#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <limits.h>
#include "random"
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>



using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

float maxFloat = std::numeric_limits<float>::max();


struct Photon {
    vec4 hitPos;
    vec3 color;
    char phi, theta;
    short flag;
};

struct Light{
    vec4 position;
    vec3 color;
};

struct Intersection
{
  vec4 position;
  float distance;
  int triangleIndex;
};

struct Camera{
  vec4 position;
  mat4 basis;
  vec3 center;
};

const int SPECULAR = 0;
const int ABSORPTION = 1;
const int DIFFUSE = 2;
const int TRANSMISSION = 3;

#define PI 3.14159
#define SCREEN_WIDTH 1080
#define SCREEN_HEIGHT 1080
#define FULLSCREEN_MODE true


mat4 R;
bool escape = false;


float focal_length = SCREEN_HEIGHT / 2;


Light light = {
  .position = vec4(0, -0.7, -0.7, 1.0),
  .color =  50.f * vec3(1,1,1),
};

Camera camera = {
  .position = vec4(0,0,0, 1.0),
  .basis = mat4(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), vec4(0,0, 2, 1.0)),
  // .center = vec3(0.003724, 0.929729, 0.07459)
  .center = vec3(0,0,0)
};

float yaw = 2 * PI / 180;

// DEFINE DISTRIBUTIONS
std::default_random_engine generator;
std::uniform_real_distribution<float> distributionx(light.position[0]-1, light.position[0]+1);
std::uniform_real_distribution<float> distributiony(light.position[1]-1, light.position[1]+1);
std::uniform_real_distribution<float> distributionz(light.position[2]-1, light.position[2]+1);
std::uniform_real_distribution<float> distribution(0 , 1);


void firePhoton(vec4 s, vec4 d, vector<Photon>& photons, Photon& photon, const vector<Triangle>& triangles);
bool ClosestIntersection(vec4 s, vec4 d, const vector<Triangle>& triangles, Intersection& closestIntersection, int index=-1);
void startPhotons(vector<Photon>& photons, const vector<Triangle>& triangles);
vec4 calculateLightDirection();
int bounceType(const Photon &photon, const Triangle& triangle);
vec3 sampleDirectionVector(const float &r1, const float &r2);
void createCoordinateSystem(const vec3 &N, vec3 &Nt, vec3 &Nb);
vec4 toVec4(const vec3 x);
vec3 toVec3(const vec4 x);
vec4 diffuseDirection(Triangle& triangle);
vec4 specularReflection(Triangle& triangle, vec4& incoming);
void Update();
void Draw(screen* screen, const vector<Triangle>& triangles, const vector<Photon>& photons);
mat4 generateRotation(vec3 a);

int main(){
    screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

    vector<Triangle> triangles;
    LoadTestModel(triangles);
    vector<Photon> photons;
    startPhotons(photons, triangles);
    printf("STORED PHOTONS: %d\n", photons.size());

      while( !escape )
    {
      Update();
      Draw(screen, triangles, photons);
      SDL_Renderframe(screen);
    }

}

void Draw(screen* screen, const vector<Triangle>& triangles, const vector<Photon>& photons)
{
  /* Clear buffer */
    memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
    for (uint photon = 0; photon < photons.size(); photon++){
        vec4 photonPos = camera.basis * photons[photon].hitPos; 
        float u = focal_length * photonPos.x/photonPos.z + SCREEN_HEIGHT/2;
        float v = focal_length * photonPos.y/photonPos.z + SCREEN_HEIGHT/2;
        PutPixelSDL(screen, u, v, photons[photon].color);
    }

}



vec4 calculateLightDirection(){
    float x = distributionx(generator);
    float y = distributiony(generator);
    float z = distributionz(generator);

    vec4 randomPoint (x,y,z, 1);

    vec4 direction = randomPoint - light.position;

    return normalize(direction);

}

void startPhotons(vector<Photon>& photons, const vector<Triangle>& triangles) {
    int numberOfPhotons = 200000;

    for (int i = 0; i < numberOfPhotons; i++){
        Photon photon;
        photon.color = light.color;

        vec4 direction = calculateLightDirection();
        firePhoton(light.position, direction, photons, photon, triangles);

    }

}

void scalePower(const int& reflectionType, const Photon& incomingPhoton, Photon& outwardsPhoton, const Triangle& triangle){
    switch (reflectionType)
    {
        case DIFFUSE:
            /* code */
            outwardsPhoton.color = incomingPhoton.color * triangle.material.diffuse;
            break;
        case SPECULAR:
            outwardsPhoton.color = incomingPhoton.color * triangle.material.specular;
        default:
            break;
    }
}


void firePhoton(vec4 s, vec4 d, vector<Photon>& photons, Photon& photon, const vector<Triangle>& triangles){
    Intersection i;
    if (photon.color.x < 0.1 && photon.color.y < 0.1 && photon.color.z < 0.1) {
        printf("DEAD PHOTON\n");
        return;
    }
    if (ClosestIntersection(s, d, triangles, i)){
        std::cout<<glm::to_string(photon.color)<<std::endl;

        Triangle hitTriangle = triangles[i.triangleIndex];
        // HERE IS WHERE THE DIFFERENT BOUNCES WILL HAPPEN ONLY STORE WHEN PHOTON HITS A DIFFUSE
        switch (bounceType(photon, hitTriangle))
        {
            case DIFFUSE: 
            {
                photon.color *= hitTriangle.color;

                int theta = int(acos(d[2])*(256.0/PI) );
                if (theta>255)
                    photon.theta = 255;
                else
                    photon.theta = (unsigned char) theta;

                int phi = int(atan2(d[1],d[0])*(256.0/(2.0*PI)) );
                if (phi>255)
                    photon.phi = 255;
                else if (phi<0)
                    photon.phi = (unsigned char)(phi+256);
                else
                    photon.phi = (unsigned char)phi;
                printf("DIFFUSE\n");
                photon.hitPos = i.position;
                photons.push_back(photon);
                Photon newPhoton;
                scalePower(DIFFUSE, photon, newPhoton, hitTriangle);
                vec4 outDir = diffuseDirection(hitTriangle);
                firePhoton(i.position + toVec4(outDir) * 0.001f, outDir, photons, newPhoton, triangles);
                /* code */
                break;
            }
            case SPECULAR:
            {
                Photon newPhoton;
                scalePower(SPECULAR, photon, newPhoton, hitTriangle);
                vec4 outDir = specularReflection(hitTriangle, d);
                firePhoton(i.position, outDir, photons, newPhoton, triangles);
                printf("SPECULAR\n");
                break;
            }
            case ABSORPTION:
                printf("ABSORPTION\n");
                break;
            default:
                break;
        }

    } else {
        return;
        printf("NO HIT\n");
    }

}

int bounceType(const Photon &photon, const Triangle& triangle) {
    float randomNumber = distribution(generator);
    if (randomNumber < triangle.material.specular) {
        return SPECULAR;
    } else if (randomNumber < triangle.material.specular + triangle.material.refraction){
        return TRANSMISSION;
    } else if (randomNumber < triangle.material.specular + triangle.material.refraction + triangle.material.absorption){
        return ABSORPTION;
    } else if (randomNumber <= triangle.material.specular + triangle.material.refraction + triangle.material.absorption + triangle.material.diffuse){
        return DIFFUSE;
    } else {
        printf("THERE WAS AN ERROR AT BOUNCE TYPE\n");
        return -1;
    }
    return 0;
}




bool ClosestIntersection(vec4 s, vec4 d, const vector<Triangle>& triangles, Intersection& closestIntersection, int index){

  closestIntersection.distance = maxFloat;

  for(uint i = 0; i < triangles.size(); i++){
    if (index > -1) {
      if (index == i) continue;
      if (dot(normalize(triangles[index].normal), normalize(d)) < 0) continue;
      // float angle = acos(dot(normalize(triangles[i].normal),normalize(d)))/abs(dot(normalize(triangles[i].normal),normalize(d)));
      // if (angle >=  1.5708) continue;
    }
    Triangle triangle = triangles[i];
    vec4 v0 = triangle.v0;
    vec4 v1 = triangle.v1;
    vec4 v2 = triangle.v2;

    vec3 e1 = vec3(v1.x-v0.x,v1.y-v0.y,v1.z-v0.z);
    vec3 e2 = vec3(v2.x-v0.x,v2.y-v0.y,v2.z-v0.z);
    vec3 b = vec3(s.x-v0.x,s.y-v0.y,s.z-v0.z);

    vec3 direc(d[0], d[1], d[2]);

    mat3 A( -direc, e1, e2);

    float detA = glm::determinant(A);

    mat3 T(b, e1, e2);

    float detT = glm::determinant(T);

    float t = detT / detA;

    if (t < closestIntersection.distance && t > 0) {
      mat3 U(-direc, b, e2);
      mat3 V(-direc, e1, b);
      float detU = glm::determinant(U);
      float detV = glm::determinant(V);

      float u = detU / detA;
      float v = detV / detA;

      if (u >= 0 && v >= 0 && u + v <= 1) {
        vec3 m = vec3(v0.x, v0.y, v0.z) + u*e1 + v*e2;
        vec4 r = vec4(m.x, m.y, m.z, 1);

        closestIntersection.distance = t;
        closestIntersection.position = r;
        closestIntersection.triangleIndex = i;
      }
    }
  }

  if (closestIntersection.distance == maxFloat) return false;


  return true;
}

vec4 diffuseDirection(Triangle& triangle){
    vec3 hitNormal = toVec3(triangle.normal);
    vec3 Nt, Nb;
    createCoordinateSystem(hitNormal, Nt, Nb);
    float r1 = distribution(generator);
    float r2 = distribution(generator);
    
    vec3 sampledVector = sampleDirectionVector(r1, r2);
    vec3 sampleWorld(
      sampledVector.x * Nb.x + sampledVector.y * hitNormal.x + sampledVector.z * Nt.x,
      sampledVector.x * Nb.y + sampledVector.y * hitNormal.y + sampledVector.z * Nt.y,
      sampledVector.x * Nb.z + sampledVector.y * hitNormal.z + sampledVector.z * Nt.z);
    
    sampleWorld = normalize(sampleWorld);
    return toVec4(sampleWorld);

}

vec3 toVec3(const vec4 x){
  return vec3(x[0], x[1], x[2]);
}

vec4 toVec4(const vec3 x) {
  return vec4(x.x, x.y, x.z, 1);
}

void createCoordinateSystem(const vec3 &N, vec3 &Nt, vec3 &Nb){
  if (std::fabs(N[0]) > std::fabs(N[1])){
    Nt = vec3(N[2], 0, -N[0]) / sqrt(N[0] * N[0] + N[2] * N[2]);
  } else {
    Nt = vec3(0, -N[2], N[1]) / sqrt(N[1] * N[1] + N[2] * N[2]);
  }
  Nb = (glm::cross(N, Nt));
}

vec3 sampleDirectionVector(const float &r1, const float &r2){
  float sinTheta = sqrtf(1 - r1 * r1);
  float phi = 2 * PI * r2;
  float x = sinTheta * cosf(phi);
  float z = sinTheta * sinf(phi);
  return vec3(x, r1, z);
}

vec4 specularReflection(Triangle& triangle, vec4& incoming) {
    vec4 n = triangle.normal;
    vec4 reflectedDirection = (incoming - 2.f * (glm::dot(incoming, n)) * n);
    return normalize(reflectedDirection);
}

mat4 generateRotation(vec3 a){
  vec4 uno = vec4(cos(a[1])*cos(a[2]), cos(a[1])*sin(a[2]), -sin(a[1]), 0);
  vec4 dos = vec4(-cos(a[0])*sin(a[2]) + sin(a[0]*sin(a[1])*cos(a[2])), cos(a[0])*cos(a[2])+sin(a[0])*sin(a[1])*sin(a[2]), sin(a[0])*cos(a[1]), 0);
  vec4 tres = vec4(sin(a[0])*sin(a[2])+cos(a[0])*sin(a[1])*cos(a[2]), -sin(a[0])*cos(a[2])+cos(a[0])*sin(a[1])*sin(a[2]), cos(a[0])*cos(a[1]), 0);
  vec4 cuatro = vec4(0,0,0,1);
  return (mat4(uno, dos, tres, cuatro));
}



void Update()
{
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;
  /*Good idea to remove this*/
  std::cout << "Render time: " << dt << " ms." << std::endl;
  /* Update variables*/
  SDL_Event e;
  mat4 translation(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), vec4(0,0,0,1));

  while(SDL_PollEvent(&e))
  {
    if (e.type == SDL_KEYDOWN){
        switch (e.key.keysym.sym) {
            case SDLK_ESCAPE:
                escape = true;
                break;
            case SDLK_UP:
                // Move camera forward
                camera.basis[3][2] += 0.1;
                break;
            case SDLK_DOWN:
            // Move camera backward
                camera.basis[3][2] += -0.1;
                break;
            case SDLK_LEFT:
            // Move camera to the left
                camera.basis[3][0] += -0.1;
                break;
            case SDLK_RIGHT:
            // Move camera to the right
                camera.basis[3][0] += 0.1;
                break;
            case SDLK_n:
                camera.basis[3][1] += -0.1;
                break;
            case SDLK_m:
                camera.basis[3][1] += 0.1;
                break;
            case SDLK_d:
                //Rotate camera right;
                camera.basis =  generateRotation(vec3(0, yaw, 0)) * camera.basis;
                break;
            case SDLK_a:
                //Rotate camera left;
                camera.basis =  generateRotation(vec3(0, -yaw, 0)) * camera.basis;

                break;
            case SDLK_w:
                //Rotate camera top;
                camera.basis =  generateRotation(vec3(yaw, 0, 0)) * camera.basis;

                break;
            case SDLK_s:
                //Rotate camera down;
                camera.basis =  generateRotation(vec3(-yaw, 0, 0)) * camera.basis;
                break;
            case SDLK_q:
                // camera.basis =  generateRotation(vec3(0, 0, -yaw)) * camera.basis;
                break;
            case SDLK_e:
                // camera.basis =  generateRotation(vec3(0, 0, yaw)) * camera.basis;
                break;

            case SDLK_u:
                light.position += vec4(0, 0, 0.1, 0);
                break;
            case SDLK_j:
                light.position += vec4(0,0,-0.1,0);
                break;
            case SDLK_h:
                light.position += vec4(-0.1, 0, 0, 0);
                break;
            case SDLK_k:
                light.position += vec4(0.1, 0, 0, 0);
                break;
            default:
                break;
            }
   }
 }
}
