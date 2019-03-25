#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <limits.h>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;


#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 800
#define FULLSCREEN_MODE true
#define PI 3.14159

float maxFloat = std::numeric_limits<float>::max();

float yaw = 2 * PI / 180;

float focal_length = SCREEN_HEIGHT / 2;

struct Intersection
{
  vec4 position;
  float distance;
  bool isTriangle;
  int circleIndex;
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

Camera camera = {
  .position = vec4(0,0,-2, 1.0),
  .basis = mat4(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), vec4(0,0,0,1)),
  // .center = vec3(0.003724, 0.929729, 0.07459)
  .center = vec3(0,0,0)
};

Light light = {
  .position = vec4(0, -0.5, -0.7, 1.0),
  .color = 14.f * vec3(1,1,1),
};

// vec4 cameraPos(0, 0, -2, 1.0);
mat4 R;
bool escape = false;
bool is_lookAt = false;




/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw(screen* screen, const vector<Triangle>& triangles, const vector<Circle>& circles);
bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closestIntersection);
bool circleIntersection(vec4 start, vec4 dir, const vector<Circle>& circles, Intersection& closestIntersection);
float calA(float radius);
vec3 calB(vec3 power, float radius);
vec3 calD(vec3 r, vec3 n, vec3 power, float radius);
vec3 DirectLight(const Intersection& i, const vector<Triangle>& triangles);
vec3 circleDirectLight(const Intersection& i, const vector<Circle>& circle);




int main( int argc, char* argv[] )
{

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  vector<Triangle> triangles;
  vector<Circle> circles;

  LoadTestModel(triangles);
  LoadCircles(circles);

  while( !escape )
    {
      Update();
      Draw(screen, triangles, circles);
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );
  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen* screen, const vector<Triangle>& triangles, const vector<Circle>& circles)
{
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  vec3 colour(1.0,0.0,0.0);

  for (int row=0; row< SCREEN_HEIGHT; row++){
    for (int col = 0; col<SCREEN_WIDTH; col++ ){
      vec4 d = camera.basis * vec4(row - SCREEN_WIDTH/2, col - SCREEN_HEIGHT/2, focal_length, 1);
      // vec4 d = vec4(row - SCREEN_WIDTH/2, col - SCREEN_HEIGHT/2, focal_length, 1);
      Intersection intersect;
      Intersection lightIntersect;
      if (ClosestIntersection(camera.position, d, triangles, intersect)){
      // if(ClosestIntersection(light.position, d, triangles, lightIntersect)) {
        // vec3 light_power = DirectLight(lightIntersect, triangles);
        vec3 light_power = DirectLight(intersect, triangles);
        // PutPixelSDL(screen, row, col, light_power);
        PutPixelSDL(screen, row, col, triangles[intersect.triangleIndex].color * light_power);
      }
      if (circleIntersection(camera.position, d, circles, intersect)) {
        vec3 light_power = circleDirectLight(intersect, circles);
        PutPixelSDL(screen, row, col, circles[intersect.circleIndex].color * light_power);
      }

    }
  }

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

bool ClosestIntersection(vec4 s, vec4 d, const vector<Triangle>& triangles, Intersection& closestIntersection){
  closestIntersection.distance = maxFloat;
  // float current = closestIntersection.distance;
  for(uint i = 0; i < triangles.size(); i++){

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

    //
    // vec3 m = vec3(v0.x, v0.y, v0.z) + x[1]*e1 + x[2]*e2;
    // vec4 r = vec4(m.x, m.y, m.z, 1);

    if (x[0] < closestIntersection.distance && x[0]>0 && x[1] >= 0 && x[2] >= 0 && x[1]+x[2] <= 1){
      closestIntersection.distance = x[0];
      closestIntersection.position = s + x[0] * d;
      closestIntersection.triangleIndex = i;
    }
  }
  if (closestIntersection.distance == maxFloat) return false;

  return true;
}

bool circleIntersection(vec4 s, vec4 d, const vector<Circle>& circles, Intersection& closestIntersection) {

  // closestIntersection.distance = maxFloat;

  float current = closestIntersection.distance;

  for (uint i = 0; i < circles.size(); i++) {

    float t0, t1;

    vec4 L4 = circles[i].center - s;

    vec3 L = vec3(L4.x, L4.y, L4.z);
    vec3 d3 = vec3(d.x, d.y, d.z);

    float a = dot(d3, d3);
    float b = 2 * dot(d3, L);
    float c = dot(L, L) - pow(circles[i].radius, 2);

    if ((pow(b,2) - 4 * a * c) < 0) continue;

    t0 = abs((-b + sqrt(pow(b, 2) - 4 * a * c)) / (2 * a));
    t1 = abs((-b - sqrt(pow(b, 2) - 4 * a * c)) / (2 * a));


    if (t0 > t1) std::swap(t0, t1);

    if (t0 < 0) {
      t0 = t1;

      if (t0 < 0) continue;
    }

    if (t0 < closestIntersection.distance) {
      closestIntersection.distance = t0;
      closestIntersection.position = s + t0 * d;
      closestIntersection.isTriangle = false;
      closestIntersection.circleIndex = i;
    }
  }

  if (closestIntersection.distance == current) return false;
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

float calA(float radius) {
  return 4 * PI * (radius * radius);
}
vec3 calB(vec3 power, float radius) {
  float a = calA(radius);

  vec3 b = power / a;
  return b;
}
vec3 calD(vec3 r, vec3 n, vec3 power, float radius) {

  float r_n = dot(r, n);

  vec3 D = calB(power, radius) * glm::max(r_n, 0.f);

  return D;

}
vec3 DirectLight(const Intersection& i, const vector<Triangle>& triangles) {

  vec3 n(triangles[i.triangleIndex].normal.x, triangles[i.triangleIndex].normal.y, triangles[i.triangleIndex].normal.z);

  vec3 x = normalize(n);

  vec3 r = normalize(light.position - i.position);

  float radius = glm::distance(light.position, i.position);

  return calD(r, x, light.color, radius);
}

vec3 circleDirectLight(const Intersection& i, const vector<Circle>& circles) {

  vec4 dir = i.position - camera.position;
  vec4 normal4 = (i.position - circles[i.circleIndex].center);
  vec3 normal = normalize(vec3(normal4));
  // circles[i.circleIndex].normal = normal;

  vec3 r = normalize(light.position - i.position);

  float radius = glm::distance(light.position, i.position);

  return calD(r, normal, light.color, radius);
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