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


#define SCREEN_WIDTH 50
#define SCREEN_HEIGHT 50
#define FOCAL_LENGTH 25
#define FULLSCREEN_MODE true
#define PI 3.141592

float maxFloat = std::numeric_limits<float>::max();

struct Intersection
{
  vec4 position;
  float distance;
  int triangleIndex;
};

struct Camera{
  vec4 position;
  mat4 basis;
};

struct Light{
  vec4 lightPos;
  vec3 lightColor;
};

Camera camera = {
  .position = vec4(0,0,-2, 1.0),
  .basis = mat4(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), vec4(0,0,0,1))
};

Light light = {
  .lightPos = vec4(0, -0.5, -0.7, 1.0),
  .lightColor = 14.f * vec3(1,1,1)

};
// vec4 cameraPos(0, 0, -2, 1.0);
float rototo = 0.25;
bool escape = false;



/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw(screen* screen, const vector<Triangle>& triangles);
bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closestIntersection);
void Inter(vec4 s, const Triangle triangle, Intersection& intersect);
vec3 DirectLight(const Intersection& i, const vector<Triangle>& triangles);



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

  vec3 colour(1.0,0.0,0.0);

  for (int row=0; row< SCREEN_HEIGHT; row++){
    for (int col = 0; col<SCREEN_WIDTH; col++ ){
      vec4 d = camera.basis * vec4(row - SCREEN_WIDTH/2, col - SCREEN_HEIGHT/2, FOCAL_LENGTH, 1);
      Intersection intersect;
      if (ClosestIntersection(camera.position, d, triangles, intersect)){
        vec3 color = DirectLight(intersect, triangles);
        PutPixelSDL(screen, row, col, color*triangles[intersect.triangleIndex].color);
      }
    }
  }

}

vec3 DirectLight(const Intersection& i, const vector<Triangle>& triangles){
  float r = glm::distance(light.lightPos, i.position);
  float area = 4 * PI * pow(r, 2);

  vec4 normal = normalize(triangles[i.triangleIndex].normal);
  vec4 direction = normalize(light.lightPos - i.position );

  float r_n = dot(direction, normal);

  vec3 d = (light.lightColor * max((r_n), 0.f))/area;

  return d;
}

bool ClosestIntersection(vec4 s, vec4 d, const vector<Triangle>& triangles, Intersection& closestIntersection){
  closestIntersection.distance = maxFloat;
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


    vec3 m = vec3(v0.x, v0.y, v0.z) + x[1]*e1 + x[2]*e2;
    vec4 r = vec4(m.x, m.y, m.z, 1);

    if (x[1] >= 0 && x[2] >= 0 && x[1]+x[2] <= 1){
      if (x[0] < closestIntersection.distance && x[0]>0){
        closestIntersection.distance = x[0];
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

mat4 lookAt(vec3 from, vec3 to){
  vec3 forward = normalize(from - to);
  vec3 right   = cross(normalize(vec3(0,1,0)), forward);
  vec3 up      = cross(forward, right);

  vec4 f = vec4(forward[0], forward[1], forward[2], 0);
  vec4 r = vec4(right[0], right[1], right[2], 0);
  vec4 u  = vec4(up[0], up[1], up[2], 0);

  vec4 fr = vec4(from[0], from[1], from[2], 1);


  mat4  camToWorld = mat4(r, u, f, fr);
  return camToWorld;

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
      case SDLK_UP:
        // Move camera forward
        translation[3][2] = 0.5;
        camera.position = translation*camera.position;
        break;
      case SDLK_DOWN:
      // Move camera backward
        translation[3][2] = -0.5;
        camera.position = translation*camera.position;
        break;
      case SDLK_LEFT:
      // Move camera to the left
        translation[3][0] = -0.5;
        camera.position = translation*camera.position;
        break;
      case SDLK_RIGHT:
      // Move camera to the right
        translation[3][0] = 0.5;
        camera.position = translation*camera.position;
        break;
      case SDLK_p:
        translation[3][1] = 0.5;
        camera.position = translation*camera.position;
        break;
      case SDLK_l:
        translation[3][1] = -0.5;
        camera.position = translation*camera.position;
        break;
      case SDLK_d:
        // Rotate camera up;
        camera.basis = translation * generateRotation(vec3(0, rototo, 0)) * camera.basis;
        break;
      case SDLK_a:
        // Rotate camera up;
        camera.basis = translation * generateRotation(vec3(0, -rototo, 0)) * camera.basis;
        break;
      case SDLK_w:
        // Rotate camera up;
        camera.basis = translation * generateRotation(vec3(rototo, 0, 0)) * camera.basis;
        break;
      case SDLK_s:
        // Rotate camera up;
        camera.basis = translation * generateRotation(vec3(-rototo, 0, 0)) * camera.basis;
        break;
      default:
        break;
    }
   }
 }
 camera.basis = lookAt(vec3(0, 0, -1), vec3(camera.position[0], camera.position[1], camera.position[2]));

}
