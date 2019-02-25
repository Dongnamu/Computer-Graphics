#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <limits.h>
#include "random"

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
  .position = vec4(0,0,-2, 1.0),
  .basis = mat4(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), vec4(0,0,0,1)),
  // .center = vec3(0.003724, 0.929729, 0.07459)
  .center = vec3(0,0,0)
};

Light light = {
  .position = vec4(0, -0.5, -0.7, 1.0),
  .color = 20.f * vec3(1,1,1),
};

Options options = {
  .bias = 1e-3,
  .indirectLight = 0.01f*vec3(1,1,1)
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
  for (int row=0; row< SCREEN_HEIGHT; row++){
    for (int col = 0; col<SCREEN_WIDTH; col++ ){
      vec3 c = processingPart(row, col, triangles);
      printf("Aperture: %f\n, Focal Distance %f\n", aperture, focalDistance);
      PutPixelSDL(screen, row, col, c);
    }
  }
}


vec3 processingPart(int row, int col, const vector<Triangle>& triangles) {
  Intersection intersect;
  vec3 color(0.f, 0.f, 0.f);
  vec3 shadow(0.f, 0.f, 0.f);
  vec4 d = camera.basis * vec4(row - SCREEN_WIDTH/2, col - SCREEN_HEIGHT/2, focal_length, 1);
  vec3 blurrColor = focusGaussian(triangles, row, col, normalize(d));
  color += blurrColor;
  return color + options.indirectLight;
}

vec3 focusGaussian(const vector<Triangle>& triangles, int row, int col, vec4 principalDirection) {
  vec3 color(0.f,0.f,0.f);
  float hitNumber = 0.f;
  vec3 shadow(0.f, 0.f, 0.f);
  vec4 target = camera.position + focalDistance * principalDirection;
  for (float x = -aperture; x <= aperture; x+= aperture){
    for (float y = -aperture; y <= aperture; y+= aperture){
      if (x == 0 && y == 0) printf("centered\n");
      vec4 randomPoint = vec4(camera.position[0] + x, camera.position[1] + y, camera.position[2], camera.position[3]);
      vec4 direction = target - randomPoint;
      Intersection inter;
      if (ClosestIntersection(randomPoint, direction, triangles, inter)){
        hitNumber += 1.0f;
        color += triangles[inter.triangleIndex].color;
        shadow+= fadedShadows(inter, triangles);
      }
    }
  }
  // printf("%f\n", hitNumber);
  if (hitNumber == 0.f) return color;
  else return color/hitNumber * shadow/hitNumber;
}


vec3 fadedShadows(const Intersection& i, const vector<Triangle>& triangles){
  vec3 light_power = DirectLight(i, triangles);
  float rayCounter = 0.0f;
  for (float j = 0.01f; j <= 0.02; j+=0.01) {
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


bool ClosestIntersection(vec4 s, vec4 d, const vector<Triangle>& triangles, Intersection& closestIntersection, int index){
  closestIntersection.distance = maxFloat;

  for(uint i = 0; i < triangles.size(); i++){
    if (index > -1) {
      if (index == i) continue;
      if (dot(normalize(triangles[index].normal), normalize(d)) < 0) continue;
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
    vec3 x = glm::inverse( A ) * b;

    vec3 m = vec3(v0.x, v0.y, v0.z) + x[1]*e1 + x[2]*e2;
    vec4 r = vec4(m.x, m.y, m.z, 1);

    if (x[0] < closestIntersection.distance && x[0]>0 && x[1] >= 0 && x[2] >= 0 && x[1]+x[2] <= 1){
      closestIntersection.distance = x[0];
      closestIntersection.position = r;
      closestIntersection.triangleIndex = i;
    }
  }
  if (closestIntersection.distance == maxFloat) {
    return false;
  }

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
  vec4 direction = normalize(light.position + vec4(distributionx(generator), distributionx(generator), distributionx(generator), 0) - i.position) ;

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
        translation[3][2] = 0.1;
        camera.position = translation*camera.position;
        break;
      case SDLK_DOWN:
      // Move camera backward
        translation[3][2] = -0.1;
        camera.position = translation*camera.position;
        break;
      case SDLK_LEFT:
      // Move camera to the left
        translation[3][0] = -0.1;
        camera.position = translation*camera.position;
        break;
      case SDLK_RIGHT:
      // Move camera to the right
        translation[3][0] = 0.1;
        camera.position = translation*camera.position;
        break;
      case SDLK_n:
        translation[3][1] = -0.1;
        camera.position = translation*camera.position;
        break;
      case SDLK_m:
        translation[3][1] = 0.1;
        camera.position = translation * camera.position;
        break;
      case SDLK_d:
        // Rotate camera right;
        camera.basis = translation * generateRotation(vec3(0, yaw, 0)) * camera.basis;
        if (is_lookAt) camera.position = translation * generateRotation(vec3(0, yaw, 0)) * camera.position;
        break;
      case SDLK_a:
        // Rotate camera left;
        camera.basis = translation * generateRotation(vec3(0, -yaw, 0)) * camera.basis;
        if (is_lookAt) camera.position = translation * generateRotation(vec3(0, -yaw, 0)) * camera.position;

        break;
      case SDLK_w:
        // Rotate camera top;
        camera.basis = translation * generateRotation(vec3(yaw, 0, 0)) * camera.basis;
        if (is_lookAt) camera.position = translation * generateRotation(vec3(yaw, 0, 0)) * camera.position;

        break;
      case SDLK_s:
        // Rotate camera down;
        camera.basis = translation * generateRotation(vec3(-yaw, 0, 0)) * camera.basis;
        if (is_lookAt) camera.position = translation * generateRotation(vec3(-yaw, 0, 0)) * camera.position;
        break;
      case SDLK_q:
        camera.basis = translation * generateRotation(vec3(0, 0, -yaw)) * camera.basis;
        if (is_lookAt) camera.position = translation * generateRotation(vec3(0, 0, -yaw)) * camera.position;
        break;
      case SDLK_e:
        camera.basis = translation * generateRotation(vec3(0, 0, yaw)) * camera.basis;
        if (is_lookAt) camera.position = translation * generateRotation(vec3(0, 0, yaw)) * camera.position;
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
