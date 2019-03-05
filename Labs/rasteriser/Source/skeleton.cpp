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


#define SCREEN_WIDTH 1080
#define SCREEN_HEIGHT 1080
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
  int x;
  int y;
  float zinv;
  vec4 pos3d;
};

struct Vertex {
  vec4 position;
};

struct Light {
  vec4 position;
  vec3 color;
  vec3 indirectLight;
};

struct Current {
  vec4 currentNormal;
  vec3 currentReflectance;
};

Current current = {
  .currentNormal = vec4(0,0,0,0),
  .currentReflectance = vec3(0,0,0)
};

Camera camera = {
  .position = vec4(0,0,-3.001, 1.0),
  .basis = mat4(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), vec4(0,0,0,1)),
  // .center = vec3(0.003724, 0.929729, 0.07459)
  .center = vec3(0,0,0)
};

Light light = {
  .position = vec4(0, -0.5, -0.7, 1.0),
  .color = 14.0f * vec3(1, 1, 1),
  .indirectLight = 0.5f * vec3(1, 1, 1)
};

float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];



// vec4 cameraPos(0, 0, -2, 1.0);
bool escape = false;
bool is_lookAt = false;


/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw(screen* screen, const vector<Triangle>& triangles);
void VertexShader(const vec4& v, Pixel& p, Vertex& vertex);
// void DrawLineSDL(screen* surface, ivec2 a, ivec2 b, vec3 color);
// void DrawPolygonEdges(screen* screen, const vector<vec4>& vertices);
void ComputePolygonRows(const vector<Pixel>& vertextPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawRows(screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
void DrawPolygon( screen* screen, const vector<vec4>& vertices, vector<Vertex>& vertex);
void PixelShader( screen* screen, const Pixel& p);


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

  for (int y = 0; y < SCREEN_HEIGHT; y++) {
    for (int x = 0; x < SCREEN_WIDTH; x++) {
      depthBuffer[y][x] = 0.0f;
    }
  }

  for (uint32_t i=0; i<triangles.size(); ++i){
    vector<vec4> vertices(3);
    vector<Vertex> vertex(3);

    vertices[0] = triangles[i].v0;
    vertices[1] = triangles[i].v1;
    vertices[2] = triangles[i].v2;

    current.currentNormal = triangles[i].normal;
    current.currentReflectance = triangles[i].color;

    DrawPolygon(screen, vertices, vertex);
    // DrawPolygonEdges(screen, vertices);

  }
}

void VertexShader(const vec4& v, Pixel& p, Vertex& vertex){
  vertex.position = camera.basis *(v - camera.position);
  p.zinv = 1/vertex.position[2];
  p.x = focal_length*(vertex.position[0]/vertex.position[2]) + SCREEN_WIDTH/2;
  p.y = focal_length*(vertex.position[1]/vertex.position[2]) + SCREEN_HEIGHT/2;
  p.pos3d = v;
}

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result){
  int N = result.size();

  float stepX = (b.x - a.x)/float(max(N-1, 1));
  float stepY = (b.y - a.y)/float(max(N-1, 1));
  float stepZ = (b.zinv - a.zinv)/float(max(N-1, 1));
  vec4 stepPos = ((b.pos3d * b.zinv) - (a.pos3d * a.zinv))/float(max(N-1, 1));

  float currentX = a.x;
  float currentY = a.y;
  float currentZ = a.zinv;
  vec4 currentStepPos = a.pos3d * a.zinv;

  for (int i=0; i<N; i++){
    result[i].x = currentX;
    result[i].y = currentY;
    result[i].zinv = currentZ;
    result[i].pos3d = (currentStepPos / result[i].zinv);
    currentX += stepX;
    currentY += stepY;
    currentZ += stepZ;
    currentStepPos += stepPos;
  }
}

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels) {
  int maxY = -numeric_limits<int>::max();
  int minY = numeric_limits<int>::max();
  int midY;

  Pixel highest;
  Pixel lowest;
  Pixel midPoint;
  int h_i;
  int l_i;
  bool is_left = false;

  for (int i = 0; i < 3; i++) {
    if (vertexPixels[i].y > maxY) {
      maxY = vertexPixels[i].y;
      highest = vertexPixels[i];
      h_i = i;
    }
    if (vertexPixels[i].y < minY) {
      minY = vertexPixels[i].y;
      lowest = vertexPixels[i];
      l_i = i;
    }
  }

  for (int i = 0; i < 3; i++) {
    if (i != h_i && i != l_i) {
      midPoint = vertexPixels[i];
      midY = vertexPixels[i].y;
    }
  }

  int numberOfRows = maxY - minY + 1;
  int numberOfhigh2mid = maxY - midY + 1;
  int numberOfmid2low = midY - minY + 1;

  leftPixels.resize(numberOfRows);
  rightPixels.resize(numberOfRows);

  for (int i = 0; i < numberOfRows; ++i) {
    leftPixels[i].x = numeric_limits<int>::max();
    rightPixels[i].x = -numeric_limits<int>::max();
  }

  vector<Pixel> longest;
  vector<Pixel> mid2high;
  vector<Pixel> low2mid;

  longest.resize(numberOfRows);
  mid2high.resize(numberOfhigh2mid);
  low2mid.resize(numberOfmid2low);

  Interpolate(lowest, highest, longest);
  Interpolate(midPoint, highest, mid2high);
  Interpolate(lowest, midPoint, low2mid);


  if (longest[numberOfmid2low - 1].x < mid2high[0].x) {
    is_left = true;
  }

  if (is_left) {
    leftPixels = longest;

    if (midY == maxY) {
      rightPixels = low2mid;
    } else {

      for (int i = 0; i < numberOfmid2low; i++) {
        rightPixels[i] = low2mid[i];
      }

      for (int i = numberOfmid2low + 1; i < numberOfRows + 1; i++) {
        rightPixels[i - 1] = mid2high[i - numberOfmid2low];
      }
    }

  } else {
    rightPixels = longest;

    if (midY == maxY) {
      leftPixels = low2mid;
    } else {
      for (int i = 0; i < numberOfmid2low; i++) {
        leftPixels[i] = low2mid[i];
      }

      for (int i = numberOfmid2low + 1; i < numberOfRows + 1; i++) {
        leftPixels[i - 1] = mid2high[i - numberOfmid2low];
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


/*Place updates of parameters here*/
mat4 generateRotation(vec3 a){
  vec4 uno = vec4(cos(a[1])*cos(a[2]), cos(a[1])*sin(a[2]), -sin(a[1]), 0);
  vec4 dos = vec4(-cos(a[0])*sin(a[2]) + sin(a[0]*sin(a[1])*cos(a[2])), cos(a[0])*cos(a[2])+sin(a[0])*sin(a[1])*sin(a[2]), sin(a[0])*cos(a[1]), 0);
  vec4 tres = vec4(sin(a[0])*sin(a[2])+cos(a[0])*sin(a[1])*cos(a[2]), -sin(a[0])*cos(a[2])+cos(a[0])*sin(a[1])*sin(a[2]), cos(a[0])*cos(a[1]), 0);
  vec4 cuatro = vec4(0,0,0,1);
  return (mat4(uno, dos, tres, cuatro));
}

void DrawRows(screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels) {

  for (uint i = 0; i < leftPixels.size(); i++) {
    int left = leftPixels[i].x;
    int right = rightPixels[i].x;

    Pixel depth_left = leftPixels[i];
    Pixel depth_right = rightPixels[i];

    int distance = abs(right - left) + 1;
    vector<Pixel> drawRow(distance);
    // drawRow.resize(distance);

    Interpolate(depth_left, depth_right, drawRow);
    int k = 0;
    for (int j = left; j < right; j++) {
      if (j >= 0 && j < SCREEN_WIDTH) {
        if (leftPixels[i].y >= 0 && leftPixels[i].y < SCREEN_HEIGHT) {
          PixelShader(screen, drawRow[k]);
        }
      }
      k++;
    }
  }
}

void PixelShader(screen* screen, const Pixel& p) {

  int x = p.x;
  int y = p.y;

  if (p.zinv >= 0) {
    if (depthBuffer[y][x] < p.zinv) {
      depthBuffer[y][x] = p.zinv;
      float r = glm::distance(light.position, p.pos3d);
      float area = 4 * PI * pow(r, 2);

      vec4 normal = current.currentNormal;
      vec4 direction = normalize(light.position - p.pos3d);

      float r_n  = dot(direction, normal);

      vec3 illumination = ((light.color * max((r_n), 0.f))/area) + light.indirectLight;
      PutPixelSDL(screen, x, y, current.currentReflectance * illumination);
    }
  }
}

void DrawPolygon(screen* screen, const vector<vec4>& vertices, vector<Vertex>& vertex) {
  int V = vertices.size();

  vector<Pixel> vertexPixels(V);

  for (int i = 0; i < V; i++) {
    VertexShader(vertices[i], vertexPixels[i], vertex[i]);
  }

  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;
  ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
  DrawRows(screen, leftPixels, rightPixels);
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
