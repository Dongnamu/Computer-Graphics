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
using glm::vec2;
using glm::mat4;
using glm::ivec2;


#define SCREEN_WIDTH 500
#define SCREEN_HEIGHT 500
#define FULLSCREEN_MODE false
#define PI 3.14159

float maxFloat = std::numeric_limits<float>::max();

float yaw = 2 * PI / 180;

float focal_length = SCREEN_HEIGHT / 2;



struct Camera{
  vec4 position;
  mat4 basis;
  vec3 center;
};

struct Pixel {
  ivec2 position;
  // float zinv;
};

Camera camera = {
  .position = vec4(0,0,-3.001, 1.0),
  .basis = mat4(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), vec4(0,0,0,1)),
  // .center = vec3(0.003724, 0.929729, 0.07459)
  .center = vec3(0,0,0)
};




// vec4 cameraPos(0, 0, -2, 1.0);
mat4 R;
bool escape = false;
bool is_lookAt = false;


/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw(screen* screen, const vector<Triangle>& triangles);
void VertexShader(const vec4& v, Pixel& p);
// void DrawLineSDL(screen* surface, ivec2 a, ivec2 b, vec3 color);
// void DrawPolygonEdges(screen* screen, const vector<vec4>& vertices);
void ComputePolygonRows(const vector<Pixel>& vertextPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawRows(screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels, vec3 color);
void DrawPolygon( screen* screen, const vector<vec4>& vertices, vec3 color);




int main( int argc, char* argv[] )
{
  vector<Pixel> vertexPixels(3);
  Pixel v1 = {.position = ivec2(10, 5)};
  Pixel v2 = {.position = ivec2( 5,10)};
  Pixel v3 = {.position = ivec2( 15,15)};

  vertexPixels[0] = v1;
  vertexPixels[1] = v2;
  vertexPixels[2] = v3;
  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;
  ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
  for( int row=0; row<leftPixels.size(); ++row )
  {
    cout << "Start: ("
    << leftPixels[row].position.x << ","
    << leftPixels[row].position.y << "). "
    << "End: ("
    << rightPixels[row].position.x << ","
    << rightPixels[row].position.y << "). " << endl;
  }
  //
  // screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  // vector<Triangle> triangles;
  // LoadTestModel(triangles);
  //
  // while( !escape )
  //   {
  //     Update();
  //     Draw(screen, triangles);
  //     SDL_Renderframe(screen);
  //   }
  //
  // SDL_SaveImage( screen, "screenshot.bmp" );
  // KillSDL(screen);
  // return 0;
}

/*Place your drawing here*/
void Draw(screen* screen, const vector<Triangle>& triangles)
{
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  for (uint32_t i=0; i<triangles.size(); ++i){
    vector<vec4> vertices(3);
    vertices[0] = triangles[i].v0;
    vertices[1] = triangles[i].v1;
    vertices[2] = triangles[i].v2;

    DrawPolygon(screen, vertices, triangles[i].color);
    // DrawPolygonEdges(screen, vertices);

  }
}

void VertexShader(const vec4& v, Pixel& p){
  vec4 n = camera.basis *(v - camera.position);
  // p.zinv = 1/n[2];
  p.position.x = focal_length*(n[0]/n[2]) + SCREEN_WIDTH/2;
  p.position.y = focal_length*(n[1]/n[2]) + SCREEN_HEIGHT/2;
}

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result){
  int N = result.size();
  Pixel step;
  // step.position = vec2(b.position - a.position)/float(max(N-1, 1));
  step.position.x = (b.position.x - a.position.x)/float(max(N-1, 1));
  // step.zinv = (b.zinv - a.zinv)/float(max(N-1, 1));
  printf("STEP X: %d\n", step.position.x);

  Pixel current = a;
  for (int i=0; i<N; i++){
    result[i] = current;
    current.position += step.position;
    printf("%d\n", current.position.x);
    // current.zinv += step.zinv;
    }
}

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels) {
  int largestVal = -1000000;
  int smallestVal = 1000000;

  int largeIndex;
  int smallIndex;
  int otherIndex;

  for (int v = 0; v < 3; ++v){
    if (vertexPixels[v].position.y > largestVal) {
      largestVal = vertexPixels[v].position.y;
      largeIndex = v;
    }
    if (vertexPixels[v].position.y < smallestVal) {
      smallestVal = vertexPixels[v].position.y;
      smallIndex = v;
    }
  }

  for (int i = 0; i < 3; i++) {
    if (i != largeIndex && i != smallIndex){
      otherIndex = i;
      break;
    }
  }
  printf("%d\n", largeIndex);
  printf("%d\n", smallIndex);
  printf("%d\n", otherIndex);

  int rows = largestVal - smallestVal + 1;
  leftPixels.resize(rows);
  rightPixels.resize(rows);

  for (int i = 0; i < rows; i++){
    leftPixels[i].position = (ivec2(+numeric_limits<int>::max(),0));
    rightPixels[i].position = (ivec2(-numeric_limits<int>::max(),0));
  }

  int toprows = vertexPixels[largeIndex].position.y  - vertexPixels[otherIndex].position.y + 1;
  int botrows = vertexPixels[otherIndex].position.y  - vertexPixels[smallIndex].position.y + 1;

  vector<Pixel> edge1(toprows);
  vector<Pixel> edge2(botrows);
  vector<Pixel> bigEdge(rows);

  Interpolate(vertexPixels[otherIndex],vertexPixels[largeIndex], edge1);
  Interpolate(vertexPixels[smallIndex], vertexPixels[otherIndex], edge2);
  Interpolate(vertexPixels[smallIndex],vertexPixels[largeIndex], bigEdge);

  if (edge1[0].position.x > bigEdge[botrows-1].position.x){
    for (int i = 0; i < rows; i++){
      if (i < botrows) rightPixels[i] = edge2[i];
      if (i >= botrows) rightPixels[i] = edge1[i-botrows+1];
      if (leftPixels[i].position.x >= edge2[i].position.x ) leftPixels[i] = bigEdge[i];
    }
  } else {
    printf("HERE");
    for (int i = 0; i < rows; i++){
      if (i < botrows) leftPixels[i] = edge2[i];
      if (i >= botrows) leftPixels[i] = edge1[i-botrows+1];
      if (rightPixels[i].position.x <= edge2[i].position.x ) rightPixels[i] = bigEdge[i];
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


/*Place updates of parameters here*/
mat4 generateRotation(vec3 a){
  vec4 uno = vec4(cos(a[1])*cos(a[2]), cos(a[1])*sin(a[2]), -sin(a[1]), 0);
  vec4 dos = vec4(-cos(a[0])*sin(a[2]) + sin(a[0]*sin(a[1])*cos(a[2])), cos(a[0])*cos(a[2])+sin(a[0])*sin(a[1])*sin(a[2]), sin(a[0])*cos(a[1]), 0);
  vec4 tres = vec4(sin(a[0])*sin(a[2])+cos(a[0])*sin(a[1])*cos(a[2]), -sin(a[0])*cos(a[2])+cos(a[0])*sin(a[1])*sin(a[2]), cos(a[0])*cos(a[1]), 0);
  vec4 cuatro = vec4(0,0,0,1);
  return (mat4(uno, dos, tres, cuatro));
}

void DrawRows(screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels, vec3 color) {

  for (uint i = 0; i < leftPixels.size(); i++) {
      int left = leftPixels[i].position.x;
      int right = rightPixels[i].position.x;

      for (int j = left; j < right; j++) {
        PutPixelSDL(screen, j, leftPixels[i].position.y, color);
      }
  }
}

void DrawPolygon(screen* screen, const vector<vec4>& vertices, vec3 color) {
  int V = vertices.size();

  vector<Pixel> vertexPixels(V);

  for (int i = 0; i < V; i++) {
    VertexShader(vertices[i], vertexPixels[i]);
  }

  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;
  ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
  DrawRows(screen, leftPixels, rightPixels, color);
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
      // case SDLK_u:
      //   light.position += vec4(0, 0, 0.1, 0);
      //   break;
      // case SDLK_j:
      //   light.position += vec4(0,0,-0.1,0);
      //   break;
      // case SDLK_h:
      //   light.position += vec4(-0.1, 0, 0, 0);
      //   break;
      // case SDLK_k:
      //   light.position += vec4(0.1, 0, 0, 0);
      //   break;
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
