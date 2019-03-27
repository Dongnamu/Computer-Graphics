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





#define SCREEN_WIDTH 256
#define SCREEN_HEIGHT 256

#define FULLSCREEN_MODE true
#define PI 3.14159

float maxFloat = std::numeric_limits<float>::max();

float yaw = 2 * PI / 180;

float focal_length = SCREEN_HEIGHT/2;

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

struct Light{
    vec4 position;
    vec3 color;
};
struct Options{
  float bias;
  vec3 indirectLight;
};
Camera camera = {
  .position = vec4(0,0,0, 1.0),
  .basis = mat4(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), vec4(0,0, -2, 1.0)),
  // .center = vec3(0.003724, 0.929729, 0.07459)
  .center = vec3(0,0,0)
};

Light light = {
  .position = vec4(0, -0.7, -0.7, 1.0),
  .color = 50.f * vec3(1,1,1),
};

Options options = {
  .bias = 0.1f,
  .indirectLight = 0.01f*vec3(0,0,0)
};

// vec4 cameraPos(0, 0, -2, 1.0);
mat4 R;
bool escape = false;
bool is_lookAt = false;

float focalDistance = 0.058f;


float aperture = 0.0003f;
std::default_random_engine generator;
std::uniform_real_distribution<float> distributionx(-0.005, 0.005);
std::uniform_real_distribution<float> distributiony(-aperture, aperture);
std::uniform_real_distribution<float> distribution(0, 1); 



/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw(screen* screen, const vector<Triangle>& triangles);
bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closestIntersection, int index=-1);
float calA(float radius);
vec3 calB(vec3 power, float radius);
vec3 calD(vec3 r, vec3 n, vec3 power, float radius);
vec3 DirectLight(const Intersection& i, const vector<Triangle>& triangles);
vec3 fadedShadows(const Intersection& i, const vector<Triangle>& triangles);
vec3 focusGaussian(const vector<Triangle>& triangles, int row, int col, vec4 principalDirection);
vec3 processingPart(int row, int col, const vector<Triangle>& triangles);
vec3 IndirectLight(const Intersection& i, const vector<Triangle>& triangles, int depth);
vec3 sampleDirectionVector(const float &r1, const float &r2);
void createCoordinateSystem(const vec3 &N, vec3 &Nt, vec3 &Nb);
vec4 toVec4(const vec3 x);
vec3 toVec3(const vec4 x);
vec3 castRay(vec4 s, vec4 d, const vector<Triangle>& triangles, int bounces);



int main( int argc, char* argv[] )
{

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  vector<Triangle> triangles;
  LoadTestModel(triangles);

  while( !escape )
    {
      Update();
      Draw(screen, triangles);
      SDL_Renderframe(screen);
      break;
    }

  SDL_SaveImage( screen, "screenshot.bmp" );
  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen* screen, const vector<Triangle>& triangles)
{
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  vec3 emergencyColor(0.f,0.f,0.f);
  for (int row=SCREEN_HEIGHT-1; row > SCREEN_HEIGHT/2; row--){
    for (int col = SCREEN_HEIGHT-1; col > SCREEN_HEIGHT/2; col--){
      vec3 c = processingPart(row, col, triangles);
      // printf("Aperture: %f\n, Focal Distance %f\n", aperture, focalDistance);
      PutPixelSDL(screen, row, col, c);
      SDL_Renderframe(screen);

    }
    // printf("ROWS: %d\n", row);
  }
}


vec3 processingPart(int row, int col, const vector<Triangle>& triangles) {
  Intersection intersect;
  vec3 color(0.f, 0.f, 0.f);
  vec3 shadow(0.f, 0.f, 0.f);
  vec4 d = camera.basis * vec4(row - SCREEN_WIDTH/2, col - SCREEN_HEIGHT/2, focal_length, 1);

  // color = castRay(camera.basis[3],d, triangles, 0);

  vec3 blurrColor = focusGaussian(triangles, row, col, normalize(d));
  color += blurrColor;
  return color/2.f;
}

vec3 focusGaussian(const vector<Triangle>& triangles, int row, int col, vec4 principalDirection) {
  vec3 color(0.f,0.f,0.f);
  float hitNumber = 0.f;
  vec3 shadow(0.f, 0.f, 0.f);
  vec4 target = camera.basis[3] + focalDistance * principalDirection;
  for (float x = -aperture; x <= aperture; x+= aperture){
    for (float y = -aperture; y <= aperture; y+= aperture){
      // if (x == 0 && y == 0) printf("centered\n");
      vec4 randomPoint = vec4(camera.basis[3][0] + x, camera.basis[3][1] + y, camera.basis[3][2], camera.basis[3][3]);
      vec4 direction = target - randomPoint;
      color += castRay(randomPoint, direction, triangles, 0);
      Intersection inter;
    }
  }
  // printf("%f\n", hitNumber);
  if (hitNumber == 0.f) {
    return color;
  }
  else return color/hitNumber;
}


vec3 fadedShadows(const Intersection& i, const vector<Triangle>& triangles){
  // vec3 light_power = IndirectLight(i, triangles, 0);
  // vec3 light_power = triangles[i.triangleIndex].color;
  // vec4 lightDirection;
  // vec3 light_power = DirectLight(i, triangles, lightDirection, 1);
  vec3 light_power (0,0,0);
  // vec3 hitColor = vec3(0.18f, 0.18f, 0.18f)/float(PI) * light_power * triangles[i.triangleIndex].color * std::max(0.f, dot(triangles[i.triangleIndex].normal, lightDirection));
  

  float rayCounter = 0.0f;
  for (float j = 0.0f; j <= 0.1; j+=0.01) {
    Intersection inter = i;
    inter.position = inter.position + j*normalize(triangles[i.triangleIndex].normal);
    light_power += DirectLight(inter, triangles);
    rayCounter++;
  }
  return light_power/rayCounter;
}

mat4 lookAt(vec3 from, vec3 to) {
  vec3 forward = normalize(from - to);
  vec3 right = normalize(cross(vec3(0,1,0), forward));
  vec3 up = cross(forward, right);

  vec4 forward4(forward.x, forward.y, forward.z, 0);
  vec4 right4(right.x, right.y, right.z, 0);
  vec4 up4(up.x, up.y, up.z, 0);
  vec4 id(0, 0, 0, 1);
  mat4 camToWorld(right4, up4, forward4, id);

  // vec4 from4(-from.x, -from.y, -from.z, 1);

  // mat4 position(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), from4);

  return camToWorld;
}

vec3 castRay(vec4 s, vec4 d, const vector<Triangle>& triangles, int bounces){
  if (bounces > 2) return vec3(0,0,0);
  Intersection i;
  vec3 hitColor(0,0,0);
  if (ClosestIntersection(s, d, triangles, i)){
    if (triangles[i.triangleIndex].isMirror){
      vec4 n = triangles[i.triangleIndex].normal;
      vec4 reflectedDirection = (d - 2.f * (glm::dot(d, n)) * n);
      hitColor += 0.8f * castRay(i.position, reflectedDirection, triangles, bounces+1);
    } else {
      vec3 directLight = DirectLight(i, triangles);
      vec3 faded = fadedShadows(i, triangles);
      directLight = faded;
      vec3 indirectLight = IndirectLight(i, triangles, bounces);
      hitColor += (directLight/float(PI) + 2.f * indirectLight) * triangles[i.triangleIndex].color/float(PI);
    }
  } else {
    return vec3(0,0,0);
  }
  return hitColor;
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


/*Place updates of parameters here*/
mat4 generateRotation(vec3 a){
  vec4 uno = vec4(cos(a[1])*cos(a[2]), cos(a[1])*sin(a[2]), -sin(a[1]), 0);
  vec4 dos = vec4(-cos(a[0])*sin(a[2]) + sin(a[0]*sin(a[1])*cos(a[2])), cos(a[0])*cos(a[2])+sin(a[0])*sin(a[1])*sin(a[2]), sin(a[0])*cos(a[1]), 0);
  vec4 tres = vec4(sin(a[0])*sin(a[2])+cos(a[0])*sin(a[1])*cos(a[2]), -sin(a[0])*cos(a[2])+cos(a[0])*sin(a[1])*sin(a[2]), cos(a[0])*cos(a[1]), 0);
  vec4 cuatro = vec4(0,0,0,1);
  return (mat4(uno, dos, tres, cuatro));
}

vec3 DirectLight(const Intersection& i, const vector<Triangle>& triangles){
  float r = glm::distance(light.position, i.position);

  

  float area = 4 * PI * pow(r, 2);

  vec4 normal = normalize(triangles[i.triangleIndex].normal);
  vec4 direction = normalize(light.position - i.position) ;
  float r_n = dot(direction, normal);

  vec3 d = (light.color * max((r_n), 0.f))/area;


  Intersection intersect;


  if (ClosestIntersection(i.position, direction, triangles, intersect, i.triangleIndex)){
    if (glm::distance(i.position, intersect.position) < r && glm::distance(i.position, intersect.position) >= 1e-10) {
      return vec3(0,0,0);
    }
  }


  return d;
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

// PATH TRACER
vec3 IndirectLight(const Intersection& i, const vector<Triangle>& triangles, int depth) {
  vec3 indirectLight (0,0,0);
  vec3 hitColor = triangles[i.triangleIndex].color;

  uint32_t N = 256;
  vec3 hitNormal = toVec3(triangles[i.triangleIndex].normal);
  vec3 Nt, Nb;
  createCoordinateSystem(hitNormal, Nt, Nb);
  float pdf = 1 / (2 * PI); 
  float sample = distribution(generator);
  // if (sample < pdf){
  //   directLight = vec3(0,0,0);
  // }
  for(uint32_t n = 0; n < N; n++){
    float r1 = distribution(generator);
    float r2 = distribution(generator);
    
    vec3 sampledVector = sampleDirectionVector(r1, r2);
    vec3 sampleWorld(
      sampledVector.x * Nb.x + sampledVector.y * hitNormal.x + sampledVector.z * Nt.x,
      sampledVector.x * Nb.y + sampledVector.y * hitNormal.y + sampledVector.z * Nt.y,
      sampledVector.x * Nb.z + sampledVector.y * hitNormal.z + sampledVector.z * Nt.z);
    sampleWorld = normalize(sampleWorld);
    Intersection intersect;

    indirectLight += r1 * castRay(i.position + toVec4(sampleWorld) * options.bias, toVec4(sampleWorld), triangles, depth + 1)/pdf;
    // if (ClosestIntersection(i.position + toVec4(sampleWorld)*options.bias, toVec4(sampleWorld), triangles, intersect, i.triangleIndex)){
    //   castRay(i.position + toVec4(sampleWorld) * options.bias, toVec4(sampleWorld), triangles, depth + 1);
    //   indirectLight += r1 * (IndirectLight(intersect, triangles, depth+1)/pdf) ;
    //   // printf("%f\n", r2);
    // }  
  }
  indirectLight /= float(N);

  return indirectLight; 

 

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
        // Rotate camera right;
        camera.basis =  generateRotation(vec3(0, yaw, 0)) * camera.basis;
        if (is_lookAt) camera.position =  generateRotation(vec3(0, yaw, 0)) * camera.position;
        break;
      case SDLK_a:
        // Rotate camera left;
        camera.basis =  generateRotation(vec3(0, -yaw, 0)) * camera.basis;
        if (is_lookAt) camera.position =  generateRotation(vec3(0, -yaw, 0)) * camera.position;

        break;
      case SDLK_w:
        // Rotate camera top;
        camera.basis =  generateRotation(vec3(yaw, 0, 0)) * camera.basis;
        if (is_lookAt) camera.position =  generateRotation(vec3(yaw, 0, 0)) * camera.position;

        break;
      case SDLK_s:
        // Rotate camera down;
        camera.basis =  generateRotation(vec3(-yaw, 0, 0)) * camera.basis;
        if (is_lookAt) camera.position =  generateRotation(vec3(-yaw, 0, 0)) * camera.position;
        break;
      case SDLK_q:
        camera.basis =  generateRotation(vec3(0, 0, -yaw)) * camera.basis;
        if (is_lookAt) camera.position =  generateRotation(vec3(0, 0, -yaw)) * camera.position;
        break;
      case SDLK_e:
        camera.basis =  generateRotation(vec3(0, 0, yaw)) * camera.basis;
        if (is_lookAt) camera.position =  generateRotation(vec3(0, 0, yaw)) * camera.position;
        break;
      case SDLK_l:
        if (is_lookAt) is_lookAt = false;
        else is_lookAt = true;
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
      case SDLK_t:
        focalDistance += 0.001;
        break;
      case SDLK_g:
        focalDistance -= 0.001;
        break;
      case SDLK_c:
        aperture -= 0.0001;
        break;
      case SDLK_v:
        aperture += 0.0001;
        break;
      default:
        break;
    }

    if (is_lookAt) {
      vec3 position = vec3(camera.position[0], camera.position[1], camera.position[2]);
      camera.basis = lookAt(camera.center, position);
    }
   }
 }
}