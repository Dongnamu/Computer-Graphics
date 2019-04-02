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
float minFloat = std::numeric_limits<float>::min();



struct Photon {
    vec4 hitPos;
    vec3 color;
    char phi, theta; 
    vec4 incomingDirection;
    int flag;
};

struct Light{
    vec4 position;
    vec3 color;
};

struct Intersection
{
  vec4 position;
  float distance;
  bool isTriangle;
  int circleIndex;
  vec4 circleNormal;
  int triangleIndex;
  bool isDummy;
};

struct Camera{
  vec4 position;
  mat4 basis;
  vec3 center;
};

struct Node{
    Photon photon;
    struct Node *left, *right;
};

struct WhichPhoton {
    unsigned long index;
    float distance;
};

const int SPECULAR = 0;
const int ABSORPTION = 1;
const int DIFFUSE = 2;
const int TRANSMISSION = 3;

#define PI 3.14159
#define SCREEN_WIDTH 300
#define SCREEN_HEIGHT 300
#define FULLSCREEN_MODE true
mat4 R;
bool escape = false;


float focal_length = SCREEN_HEIGHT / 2;


Light light = {
  .position = vec4(0, -0.99, -0.4, 1.0),
  .color = 100000.f * vec3(1,1,1),
};

Camera camera = {
  .position = vec4(0,0,0, 1.0),
  .basis = mat4(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), vec4(0,0, 2, 1.0)),
  // .center = vec3(0.003724, 0.929729, 0.07459)
  .center = vec3(0,0,0)
};

float yaw = 2 * PI / 180;

// DEFINE DISTRIBUTIONS
float md = 0.5;
std::default_random_engine generator;
std::uniform_real_distribution<float> distributionx(light.position[0]-md, light.position[0]+md);
std::uniform_real_distribution<float> distributiony(light.position[1], light.position[1]+md);
std::uniform_real_distribution<float> distributionz(light.position[2]-md, light.position[2]+md);
std::uniform_real_distribution<float> distribution(0 , 1);



void firePhoton(vec4 s, vec4 d, vector<Photon>& photons, Photon& photon, const vector<Triangle>& triangles, const vector<Circle>& circles, int n_fire);
bool ClosestIntersection(vec4 s, vec4 d, const vector<Triangle>& triangles, const vector<Circle>& circles, Intersection& closestIntersection, int index);
bool circleIntersection(vec4 s, vec4 d, const vector<Circle>& circles, Intersection& closestIntersection, int within_index = -1, bool isInCircle = false);
void startPhotons(vector<Photon>& photons, const vector<Triangle>& triangles, const vector<Circle>& circles, bool Caustics);
vec4 calculateLightDirection();
int bounceType(const Photon &photon, const Material& material);
vec3 sampleDirectionVector(const float &r1, const float &r2);
void createCoordinateSystem(const vec3 &N, vec3 &Nt, vec3 &Nb);
vec4 toVec4(const vec3 x);
vec3 toVec3(const vec4 x);
vec4 diffuseDirection(vec4& normal, float &random);
vec4 specularReflection(vec4& normal, vec4& incoming);
vec3 circleRefract(const vec4& d, const Intersection& incident, const float& indexOfRefraction);



void Update();
void Draw(screen* screen, const vector<Triangle>& triangles,  const vector<Photon>& kdTree, const vector<Circle>& circles, uint32_t  *  indices);
void scalePower(const int& reflectionType, const Photon& incomingPhoton, Photon& outwardsPhoton, const Material& material);


vec3 processingPart(int row, int col, const vector<Triangle>& triangles, const vector<Circle>& circles, const vector<Photon>& kdTree, uint32_t  *  indices);
vec3 castRay(vec4 s, vec4 d, const vector<Triangle>& triangles, const vector<Circle>& circles, int bounces,  int index, const vector<Photon> kdTree, uint32_t  *  indices);
vec3 DirectLight(const Intersection& i, const vector<Triangle>& triangles, const vector<Circle>& circles);


mat4 generateRotation(vec3 a);




int main(){
    screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );


    vector<Triangle> triangles;
    vector<Photon> photons;
    vector<Photon> caustics;
    vector<Photon> full;

    vector<Circle> circles;
    LoadTestModel(triangles);
    LoadCircles(circles);
    startPhotons(caustics, triangles, circles, true);
    // startPhotons(photons, triangles, circles, false);

    

    uint32_t indices[photons.size()];

    for (int i = 0; i < photons.size(); i++){
        full.push_back(photons[i]);
    }

    for (int i = 0; i < caustics.size(); i++){
        full.push_back(caustics[i]);
    }

    printf("STORED PHOTONS: %d\n", full.size());


      while( !escape )
    {
      Update();
      Draw(screen, triangles, full, circles, indices);
      SDL_Renderframe(screen);
      Update();
    //   break;
    }

    SDL_SaveImage(screen, "screenshot.png");
}

void Draw(screen* screen, const vector<Triangle>& triangles, const vector<Photon>& kdTree, const vector<Circle>& circles, uint32_t * indices)
{
  /* Clear buffer */
        memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
        for (uint photon = 0; photon < kdTree.size(); photon++){
            vec4 photonPos = camera.basis * kdTree[photon].hitPos; 
            float u = focal_length * photonPos.x/photonPos.z + SCREEN_HEIGHT/2;
            float v = focal_length * photonPos.y/photonPos.z + SCREEN_HEIGHT/2;
            PutPixelSDL(screen, u, v, kdTree[photon].color);
         }


    // memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
    
    // for (int row=0; row < SCREEN_HEIGHT; row++){
    //     for (int col = 0; col < SCREEN_WIDTH; col++){
            
    //         vec3 c = processingPart(row, col, triangles, circles, kdTree, indices);
    //         // printf("Aperture: %f\n, Focal Distance %f\n", aperture, focalDistance);
    //         PutPixelSDL(screen, row, col, c);

    //     }
    //     // SDL_Renderframe(screen);
    //     // printf("ROW: %d\n", row);
        
    // }

}

vec3 processingPart(int row, int col, const vector<Triangle>& triangles, const vector<Circle>& circles, const vector<Photon>& kdTree, uint32_t * indices){
    Intersection intersect;
    vec3 color(0.f, 0.f, 0.f);
    vec3 shadow(0.f, 0.f, 0.f);
    vec4 d = camera.basis * vec4(row - SCREEN_WIDTH/2, col - SCREEN_HEIGHT/2, focal_length, 1);
    Intersection inter;

    vec3 mainColor = castRay(camera.basis[3], d, triangles, circles, 0, -1, kdTree, indices);

    return mainColor;
}

void linearNN(const vector<Photon>& kdTree, const vec4 targetPoint, vector <int>& indexes, const float radius){
    for (uint i = 0; i < kdTree.size(); i++){
        if (glm::distance(kdTree[i].hitPos, targetPoint) < radius) indexes.push_back(i);
    }

}




vec3 castRay(vec4 s, vec4 d, const vector<Triangle>& triangles, const vector<Circle>& circles, int bounces, int index, const vector<Photon> kdTree, uint32_t * indices){
    if (bounces > 5) return vec3(0,0,0);
    Intersection i;
    vec3 hitColor(0,0,0);

    if (ClosestIntersection(s, d, triangles, circles, i, index)){

        vec4 normal;
        vec3 objectColor;
        bool isMirror;
        bool isGlass = false;
        Material material;
        if (i.isTriangle){
            objectColor = triangles[i.triangleIndex].color;  
            isMirror = triangles[i.triangleIndex].isMirror;
            normal = triangles[i.triangleIndex].normal;
            material = triangles[i.triangleIndex].material;
        }
        else {
            objectColor = circles[i.circleIndex].color;
            isMirror = circles[i.circleIndex].isMirror;
            normal = i.circleNormal;
            material = circles[i.circleIndex].material;
            isGlass = circles[i.circleIndex].isGlass;
        }

        if (isMirror){
            vec4 reflectedDirection = specularReflection(normal, d);
            hitColor += 0.8f * castRay(i.position, reflectedDirection, triangles, circles, bounces+1, -1, kdTree, indices);
        } else if (isGlass ) {
            vec4 refractedDirection = toVec4(circleRefract(d, i, 2));
            vec4 start_position = i.position - toVec4(vec3(0.01f * i.circleNormal));
            
            Intersection new_intersect;

            if (circleIntersection(start_position, refractedDirection, circles, new_intersect, i.circleIndex, true)) {
                new_intersect.circleNormal = -normalize(new_intersect.position - circles[i.circleIndex].center);
                vec4 new_startPosition = new_intersect.position + toVec4(vec3(0.01f * new_intersect.circleNormal));
                vec4 new_refract = toVec4(circleRefract(refractedDirection, new_intersect, 1.5));
                hitColor += castRay(new_startPosition, new_refract, triangles, circles, bounces, i.circleIndex, kdTree, indices);
            }
        } 
        else {

            vec3 indirect;
            vec3 directLight = DirectLight(i, triangles, circles);
            vec3 acc (0,0,0);
            const float radius = 0.1f;

            vector<Photon> copy = kdTree;
            vector<int> indexes;
            

            linearNN(kdTree, i.position, indexes, radius);
            int n = 10000;

            if (indexes.size() < n){
                n = indexes.size();
            }
            for (uint j = 0; j < indexes.size(); j++) {   
                float weight;
                if (i.isTriangle){
                    weight = (max(glm::dot(abs(kdTree[indexes[j]].incomingDirection), abs(normal)), 0.f));
                } else {
                    weight = abs(glm::dot((kdTree[indexes[j]].incomingDirection), normal));
                }
                acc += kdTree[indexes[j]].color * weight;
            }


            float area =   float(PI) * pow(radius, 2);
            acc /= area;
            
            hitColor =  acc  ;
            
        } 
        
    } else {
        return vec3(0,0,0);
    }

    return hitColor;
}

void firePhoton(vec4 s, vec4 d, vector<Photon>& photons, Photon& photon, const vector<Triangle>& triangles, const vector<Circle>& circles, int n_fire){
    Intersection i;
    if (n_fire > 3) return;
    // if (photon.color.x < 0.1 && photon.color.y < 0.1 && photon.color.z < 0.1) {
    //     // printf("DEAD PHOTON\n");
    //     return;
    // }
    if (ClosestIntersection(s, d, triangles, circles, i, -1)){
        vec4 objectNormal;
        vec3 objectColor;
        Material objectMaterial;
        bool isCircle = false;
        // if (i.isTriangle) return;
        if (i.isTriangle) {
            objectNormal = triangles[i.triangleIndex].normal;
            objectColor = triangles[i.triangleIndex].color;
            objectMaterial = triangles[i.triangleIndex].material;
        }
        else {
            isCircle = true;
            objectNormal = i.circleNormal;
            objectColor = circles[i.circleIndex].color;
            objectMaterial = circles[i.circleIndex].material;
        }

        // if (isCircle) std::cout<<glm::to_string(objectNormal)<<std::endl;

        // HERE IS WHERE THE DIFFERENT BOUNCES WILL HAPPEN ONLY STORE WHEN PHOTON HITS A DIFFUSE
        switch (bounceType(photon, objectMaterial))
        {
            case DIFFUSE: 
            {
                // photon.color *= objectColor;
                // printf("DIFFUSE\n");
                photon.incomingDirection = normalize((d));
                photon.hitPos = i.position;
                photons.push_back(photon);
                    
                Photon newPhoton;
                float random;
                float pdf = 1 / (2 * PI); 

                scalePower(DIFFUSE, photon, newPhoton, objectMaterial);
                vec4 outDir = diffuseDirection(objectNormal, random);
                newPhoton.color *=  objectColor ;
                
                firePhoton(i.position , outDir, photons, newPhoton, triangles, circles, n_fire + 1);

                /* code */
                break;
            }
            case SPECULAR:
            {
                Photon newPhoton;
                printf("SPECULAR\n");
                scalePower(SPECULAR, photon, newPhoton, objectMaterial);
                vec4 outDir = specularReflection(objectNormal, d);
                firePhoton(i.position, outDir, photons, newPhoton, triangles, circles, n_fire+1);
                break;
            }
            case ABSORPTION:
                photon.incomingDirection = normalize(d);
                photon.hitPos = i.position;
                photons.push_back(photon);
                // printf("ABSORPTION\n");
                break;
            case TRANSMISSION:
            {
                vec4 refractedDirection = toVec4( circleRefract(d, i, 1.5));
                vec4 offset = toVec4(vec3(0.001f * i.circleNormal));
                vec4 start_position = i.position - offset;

                Photon newPhoton;
                scalePower(TRANSMISSION, photon, newPhoton, objectMaterial);

                Intersection new_intersect;
                if (circleIntersection(start_position, refractedDirection, circles, new_intersect, i.circleIndex, true)) {
                    scalePower(TRANSMISSION, newPhoton, newPhoton, objectMaterial);
                    new_intersect.circleNormal = normalize(new_intersect.position - circles[i.circleIndex].center);   
                    vec4 new_startPosition = new_intersect.position + toVec4(vec3(0.001f * new_intersect.circleNormal));
                    vec4 new_refract = toVec4(circleRefract(refractedDirection, new_intersect, 1.5));
                    firePhoton(new_startPosition, new_refract, photons, newPhoton, triangles, circles, n_fire+1);
                } 

                break;

            }


            default:
                break;
        }

    } else {
        return;
        
        printf("NO HIT\n");
    }

}

vec3 DirectLight(const Intersection& i, const vector<Triangle>& triangles, const vector<Circle>& circles){
    float r = glm::distance(light.position, i.position);
    float area = 4 * PI * pow(r, 2);
    vec4 normal; 
    if (i.isTriangle) {
        normal = triangles[i.triangleIndex].normal;
    }
    else {
        normal = i.circleNormal;
    }

    vec4 direction = normalize(light.position - i.position);
    float r_n = dot(direction, normal);
    vec3 d = (light.color * max((r_n), 0.f))/area;
    
    Intersection intersect;
    if (ClosestIntersection(i.position, direction, triangles, circles, intersect, -1)){
        if (glm::distance(i.position, intersect.position) < r && glm::distance(i.position, intersect.position) >= 1e-10) {
        return vec3(0.1,0.1,0.1);
        }
    }

    return d;

    
}

vec3 circleRefract(const vec4& d, const Intersection& incident, const float& indexOfRefraction) {

  vec3 normal = incident.circleNormal;

  vec3 i = normalize(vec3(d));

  float NdotI = glm::clamp(dot(normal, i), -1.f, 1.f);

  float airMedium = 1;
  float materialMedium = indexOfRefraction;

  if (NdotI < 0) {
    NdotI = -NdotI;
  } else {
    normal = -incident.circleNormal;
    // std::swap(airMedium, materialMedium);
  }

  float mediumDiv = airMedium / materialMedium;

  float k = 1 - mediumDiv * mediumDiv * (1 - NdotI * NdotI);

  if (k < 0) return vec3(0,0,0);
  else return mediumDiv * i + (mediumDiv * NdotI - sqrt(k)) * normal;
  
}


vec4 calculateLightDirection(){
    float x = distributionx(generator);
    float y = distributiony(generator);
    float z = distributionz(generator);

    vec4 randomPoint (x,y,z, 1);

    vec4 direction =   randomPoint - light.position ;

    return normalize(direction);

}

vec4 calculateCausticsDirection(const vector<Circle>& circles) {
    vec4 center;
    for (uint i = 0; i < circles.size(); i++ ){
        if (circles[i].isGlass) {
            center = circles[i].center;
        }
    }
    // center = circles[0].center;
    float noise = 0.001;
    std::uniform_real_distribution<float> distributionCx(center[0]-noise, center[0]+noise);
    std::uniform_real_distribution<float> distributionCy(center[1]-noise, center[1]+noise);
    std::uniform_real_distribution<float> distributionCz(center[2]-noise, center[2]+noise);

    vec4 randomNoise (distributionCx(generator), distributionCy(generator), distributionCz(generator), 1);

    vec4 direction = randomNoise - light.position;

    return normalize(direction);
}

void startPhotons(vector<Photon>& photons, const vector<Triangle>& triangles, const vector<Circle>& circles, bool Caustics) {
    float numberOfPhotons = 10000;

    for (int i = 0; i < numberOfPhotons; i++){
        Photon photon;
        photon.color = light.color/numberOfPhotons;
        vec4 direction;

        if (Caustics) {
            direction = calculateCausticsDirection(circles);
            photon.color /= 2;

        } else {

            direction = calculateLightDirection();
        }


        firePhoton(light.position, direction, photons, photon, triangles, circles, 1);
    }

}

void scalePower(const int& reflectionType, const Photon& incomingPhoton, Photon& outwardsPhoton, const Material& material){
    switch (reflectionType)
    {
        case DIFFUSE:
            /* code */
            outwardsPhoton.color = incomingPhoton.color * material.diffuse;
            break;
        case SPECULAR:
            outwardsPhoton.color = incomingPhoton.color * material.specular;
            break;
        case TRANSMISSION:
            outwardsPhoton.color = incomingPhoton.color * material.refraction;
            break;
        default:
            break;
        
    }
}




int bounceType(const Photon &photon, const Material& material) {
    float randomNumber = distribution(generator);
    if (randomNumber < material.specular) {
        return SPECULAR;
    } else if (randomNumber < material.specular + material.refraction){
        return TRANSMISSION;
    } else if (randomNumber < material.specular + material.refraction + material.absorption){
        return ABSORPTION;
    } else if (randomNumber <= material.specular + material.refraction + material.absorption + material.diffuse){
        return DIFFUSE;
    } else {
        printf("THERE WAS AN ERROR AT BOUNCE TYPE\n");
        return -1;
    }
    return 0;
}



bool ClosestIntersection(vec4 s, vec4 d, const vector<Triangle>& triangles, const vector<Circle>& circles, Intersection& closestIntersection, int index){

  closestIntersection.distance = maxFloat;
  // float current = closestIntersection.distance;
  for(uint i = 0; i < triangles.size(); i++){
    // if (index > -1) {
    //   if (index == i) continue;
    //   if (dot(normalize(triangles[index].normal), normalize(d)) < 0) continue;
    //   // float angle = acos(dot(normalize(triangles[i].normal),normalize(d)))/abs(dot(normalize(triangles[i].normal),normalize(d)));
    //   // if (angle >=  1.5708) continue;
    // }
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
        closestIntersection.isTriangle = true;
      }
    }
  }

//   vec4 offset = vec4(0.001f, 0.001f, 0.001f, 1.f) * closestIntersection.circleNormal;
//   vec4 start_position = closestIntersection.position - offset;

//   if (!circleIntersection(s, d, circles, closestIntersection)) return false;
  
    circleIntersection(s, d, circles, closestIntersection, index, false);

  if (closestIntersection.distance == maxFloat) return false;
  return true;
}

bool circleIntersection(vec4 s, vec4 d, const vector<Circle>& circles, Intersection& closestIntersection, int within_index, bool isInCircle) {

  if (isInCircle) closestIntersection.distance = maxFloat;

  float current = closestIntersection.distance;
  for (uint i = 0; i < circles.size(); i++) {

    if (isInCircle && (i != within_index)) continue;
    if (!isInCircle && within_index != -1 && (i == within_index)) continue;
   
    float t0 = maxFloat;
    float t1 = maxFloat;

    vec4 L4 = s - circles[i].center;

    vec3 L = vec3(L4);
    vec3 d3 = vec3(d);

    float a = dot(d3 ,d3);
    float b = 2 * dot(d3, L);
    float c = dot(L, L) - pow(circles[i].radius, 2);

    if ((pow(b,2) - 4 * a * c) < 0) continue;

    // if (within_index == 2 && (pow(b,2) - 4 * a * c) > 0) printf("I'm here");
    
    

    t0 = (-b + sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);
    t1 = (-b - sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);

    if (t0 > t1) std::swap(t0, t1);

    if (t0 < 0) {
      t0 = t1;

      if (isInCircle && t0 < 0) printf("Can't find t value\n");
      if (t0 < 0) continue;
    }

    if (t0 < closestIntersection.distance) {
      // printf("I'm here\n");
      closestIntersection.distance = t0;
      closestIntersection.circleNormal = normalize((s + t0 * d) - circles[i].center);
      closestIntersection.position = (s + t0 * d) + closestIntersection.circleNormal * 0.00001f;
      closestIntersection.isTriangle = false;
      closestIntersection.circleIndex = i;
    }
  }

  if (closestIntersection.distance == current) return false;
  return true;
}


vec4 diffuseDirection(vec4& normal, float &random){
    vec3 hitNormal = toVec3(normal);
    vec3 Nt, Nb;
    createCoordinateSystem(hitNormal, Nt, Nb);
    float r1 = distribution(generator);
    float r2 = distribution(generator);
    
    vec3 sampledVector = sampleDirectionVector(r1, r2);
    vec3 sampleWorld(
      sampledVector.x * Nb.x + sampledVector.y * hitNormal.x + sampledVector.z * Nt.x,
      sampledVector.x * Nb.y + sampledVector.y * hitNormal.y + sampledVector.z * Nt.y,
      sampledVector.x * Nb.z + sampledVector.y * hitNormal.z + sampledVector.z * Nt.z);
    random = r1;
    sampleWorld = (sampleWorld);
    return normalize(toVec4(sampleWorld));

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

vec4 specularReflection(vec4& normal, vec4& incoming) {
    vec4 n = normal;
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
//   std::cout << "Render time: " << dt << " ms." << std::endl;

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