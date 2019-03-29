#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <limits.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

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


float depth = 0.55;

const float epsilon = 0.1;


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

struct Clipper {
  vec4 leftNormal;
  vec4 rightNormal;
  vec4 topNormal;
  vec4 botNormal;
  vec4 nearNormal;
  vec4 nearPoint;
  vec4 leftPoint;
  vec4 rightPoint;
  vec4 topPoint;
  vec4 botPoint;
};

Current current = {
  .currentNormal = vec4(0,0,0,0),
  .currentReflectance = vec3(0,0,0)
};

Camera camera = {
  .position = vec4(0,0,-2,1.0),
  .basis = mat4(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), vec4(0,0,0,1)),  // .center = vec3(0.003724, 0.929729, 0.07459)
  .center = vec3(0,0,0)
};

Light light = {
  .position = vec4(0, -0.5, -0.7, 1.0),
  .color = 14.0f * vec3(1, 1, 1),
  .indirectLight = 0.5f * vec3(1, 1, 1)
};

Clipper clipper = {
  .leftNormal = vec4(0,0,0,0),
  .rightNormal = vec4(0,0,0,0),
  .topNormal = vec4(0,0,0,0),
  .botNormal = vec4(0,0,0,0),
  .nearNormal = vec4(0,0,0,0),
  .nearPoint = vec4(0,0,0,0),
  .leftPoint = vec4(0,0,0,0),
  .rightPoint = vec4(0,0,0,0),
  .topPoint = vec4(0,0,0,0),
  .botPoint = vec4(0,0,0,0)
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
mat4 generateRotation(vec3 a);
void ClipTriangles(const vector<Triangle>& triangles, vector<Triangle>& clippedTriangles, vec4 normal, vec4 point);
void calSign(vec4 trianglePoint, float& value, vec4 normal, vec4 point);
void countSign(float v0, float v1, float v2, vector<int>& numSigns);
bool isBoundary(vector<int> numSigns);
bool isNegative(vector<int> numSigns);
void organiseData(vec4 point, float v, vector<vec4>& in, vector<vec4>& boundary, vector<vec4>& out);
void updateClippers();
mat4 getRotation();

int main( int argc, char* argv[] )
{

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  vector<Triangle> triangles;
  LoadTestModel(triangles);

  int size = triangles.size();

  // camera.position = mat4(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), light.position) * camera.position;
  // camera.basis =  generateRotation(vec3((90 * PI / 180), 0, 0)) * camera.basis;
  while( !escape )
    {
      vector<Triangle> nearTriangles;
      vector<Triangle> leftTriangles;
      vector<Triangle> rightTriangles;
      vector<Triangle> topTriangles;
      vector<Triangle> bottomTriangles;
      Update();
      updateClippers();
      // ClipTriangles(triangles, nearTriangles, clipper.nearNormal, clipper.nearPoint);
      // ClipTriangles(triangles, leftTriangles, clipper.leftNormal, clipper.leftPoint);
      // ClipTriangles(triangles, rightTriangles, clipper.rightNormal, clipper.rightPoint);
      // ClipTriangles(triangles, topTriangles, clipper.topNormal, clipper.topPoint);
      ClipTriangles(triangles, bottomTriangles, clipper.botNormal, clipper.botPoint);
      Draw(screen, bottomTriangles);
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

void updateClippers() {

    // These are the directions towards the four corners of the img plane
<<<<<<< HEAD
    vec4 leftUpCorner = normalize(camera.basis * vec4(SCREEN_WIDTH/2, -SCREEN_HEIGHT/2, focal_length / 4, 1));
    vec3 leftUp = vec3(leftUpCorner[0],  leftUpCorner[1], leftUpCorner[2]);
   
    vec4 leftBotCorner = normalize(camera.basis * vec4(SCREEN_WIDTH/2, SCREEN_HEIGHT/2, focal_length / 4, 1));
    vec3 leftBot = vec3(leftBotCorner[0],  leftBotCorner[1], leftBotCorner[2]);

    vec4 rightTopCorner = normalize(camera.basis * vec4(-SCREEN_WIDTH/2, -SCREEN_HEIGHT/2, focal_length / 4, 1));
    vec3 rightUp = vec3(rightTopCorner[0],  rightTopCorner[1], rightTopCorner[2]);

    vec4 rightBotCorner = normalize(camera.basis * vec4(-SCREEN_WIDTH/2, SCREEN_HEIGHT/2, focal_length / 4, 1));
=======
    vec4 leftUpCorner = normalize(getRotation() * vec4(-SCREEN_WIDTH/2, -SCREEN_HEIGHT/2, focal_length / 4, 1));
    vec3 leftUp = vec3(leftUpCorner[0],  leftUpCorner[1], leftUpCorner[2]);
   
    vec4 leftBotCorner = normalize(getRotation() * vec4(-SCREEN_WIDTH/2, SCREEN_HEIGHT/2, focal_length / 4, 1));
    vec3 leftBot = vec3(leftBotCorner[0],  leftBotCorner[1], leftBotCorner[2]);

    vec4 rightTopCorner = normalize(getRotation() * vec4(SCREEN_WIDTH/2, -SCREEN_HEIGHT/2, focal_length / 4, 1));
    vec3 rightUp = vec3(rightTopCorner[0],  rightTopCorner[1], rightTopCorner[2]);

    vec4 rightBotCorner = normalize(getRotation() * vec4(SCREEN_WIDTH/2, SCREEN_HEIGHT/2, focal_length / 4, 1));
>>>>>>> 1fd5e7724470c739227011bcd01c66d7a0e4c171
    vec3 rightBot = vec3(rightBotCorner[0],  rightBotCorner[1], rightBotCorner[2]);

    // We use those directions to get the normal of the plane
    vec3 leftNormal = (glm::cross(leftUp, leftBot));
    vec3 rightNormal =(glm::cross(rightUp, rightBot));
    vec3 topNormal = (glm::cross(leftUp, rightUp));
    vec3 botNormal = (glm::cross(leftBot, rightBot));

    
    
    vec4 leftN = vec4(leftNormal.x, leftNormal.y, leftNormal.z, 1);
    vec4 rightN = vec4(rightNormal.x, rightNormal.y, rightNormal.z, 1);
    vec4 topN = vec4(topNormal.x, topNormal.y, topNormal.z, 1);
    vec4 botN = vec4(botNormal.x, botNormal.y, botNormal.z, 1);

    // We use the directions from the top to find points in the near plane and then find the normal of this plane
    vec4 leftTopPoint = camera.position + (leftUpCorner * depth);
    vec4 leftBottomPoint = camera.position + (leftBotCorner * depth);
    vec4 rightBottomPoint = camera.position + (rightBotCorner * depth);
    vec4 rightTopPoint = camera.position + (rightTopCorner * depth);
    vec4 leftToRight = rightBottomPoint - leftBottomPoint;
    vec4 leftToTop = leftTopPoint - leftBottomPoint;
    vec4 topTotop = rightTopPoint - leftTopPoint;

    vec3 leftToRight3 = vec3(leftToRight.x, leftToRight.y, leftToRight.z);
    vec3 leftToTop3 = vec3(leftToTop.x, leftToTop.y, leftToTop.z);

    vec3 nearNormal = (glm::cross(leftToRight3, leftToTop3));
    vec4 nearN = normalize(vec4(nearNormal.x, nearNormal.y, nearNormal.z, 1));

    clipper.nearNormal = nearN;
    clipper.nearPoint = leftBottomPoint;
    clipper.leftNormal = leftN;
    clipper.rightNormal = rightN;
    clipper.topNormal = topN;
    clipper.botNormal = botN;
    clipper.leftPoint = leftTopPoint;
    clipper.rightPoint = rightTopPoint;
    clipper.botPoint = leftBottomPoint + depth * leftToRight;
    clipper.topPoint = leftTopPoint + depth * topTotop;

}

mat4 getRotation() {
  return mat4(camera.basis[0], camera.basis[1], camera.basis[2], vec4(0,0,0,1));
}


void ClipTriangles(const vector<Triangle>& triangles, vector<Triangle>& clippedTriangles, vec4 normal, vec4 point) {
  for (int i = 0; i < triangles.size(); i++){
    vec4 v0 = triangles[i].v0;
    vec4 v1 = triangles[i].v1;
    vec4 v2 = triangles[i].v2;
    float v0S;
    float v1S;
    float v2S;
    calSign(v0, v0S, normal, point);
    calSign(v1, v1S, normal, point);
    calSign(v2, v2S, normal, point);
    
    // Count how many vertices are positive, zero and negative
    vector<int> numSigns(3);
    countSign(v0S, v1S, v2S, numSigns);

    
    //If all of them are positive, add to the clipped triange
    if (numSigns[0] == 3) {
      // printf("Positive\n");
      clippedTriangles.push_back(triangles[i]);
      continue;
    }

    // If all of them are negative, skip to next triangle
    if (isNegative(numSigns)) {
      // printf("Negative\n");
      continue;
    }


    //If two points are at the boundary and one point is positive,
    // or one point is at the boundary and two points are positive,
    //add to the clipped triange
    if (isBoundary(numSigns)) {
      clippedTriangles.push_back(triangles[i]);
      continue;
    }

    vector<vec4> in;
    vector<vec4> boundary;
    vector<vec4> out;

    organiseData(v0, v0S, in, boundary, out);
    organiseData(v1, v1S, in, boundary, out);
    organiseData(v2, v2S, in, boundary, out);
    
    // printf("%d %d %d\n", numSigns[0], numSigns[1], numSigns[2]);
    std::cout<<glm::to_string(in[0])<<std::endl;
    
    if (numSigns[0] == 1) {
      vec4 inside = in[0];

      if (numSigns[1] == 1) {
        vec4 outside = out[0];

        float d1;
        float d2;
        calSign(inside, d1, normal, point);
        calSign(outside, d2, normal, point);
        
        float t = d1 / (d1 - d2);

        vec4 intersect = inside + t * (outside - inside);
        boundary.push_back(intersect);
      } else {
        vec4 outside0 = out[0];
        vec4 outside1 = out[1];

        float d1;
        float d2; 
        float d3;
        calSign(inside, d1, normal, point);
        calSign(outside0, d2, normal, point);
        calSign(outside1, d3, normal, point);

        float t0 = d1 / (d1 - d2);
        float t1 = d1 / (d1 - d3);

        vec4 intersect0 = inside + t0 * (outside0 - inside);
        vec4 intersect1 = inside + t1 * (outside1 - inside);

        boundary.push_back(intersect0);
        boundary.push_back(intersect1);
      }

      Triangle new_triangle(inside, boundary[0], boundary[1], triangles[i].color);   
      new_triangle.normal = triangles[i].normal;   
      clippedTriangles.push_back(new_triangle);



    } else {
      vec4 inside0 = in[0];
      vec4 inside1 = in[1];
      vec4 outside = out[0];

      float d1;
      float d2;
      float d3;
      calSign(inside0, d1, normal, point);
      calSign(inside1, d2, normal, point);
      calSign(outside, d3, normal, point);
      
      float t0 = d1 / (d1 - d3);
      float t1 = d2 / (d2 - d3);

      vec4 intersect0 = inside0 + t0 * (outside - inside0);
      vec4 intersect1 = inside1 + t1 * (outside - inside1);

      Triangle new_triangle0(inside0, intersect1, intersect0, triangles[i].color);
      Triangle new_triangle1(inside0, inside1, intersect1, triangles[i].color);

      new_triangle0.normal = triangles[i].normal;
      new_triangle1.normal = triangles[i].normal;

      clippedTriangles.push_back(new_triangle0);
      clippedTriangles.push_back(new_triangle1);
    }
  }
}

void calSign(vec4 trianglePoint, float& value, vec4 normal, vec4 point) {
  float calculation = (dot(normal, (trianglePoint - point)));

  if ((calculation < epsilon) && (calculation > 0.0)) {
    value = 0;
  } else {
    value = calculation;
  }
}

void countSign(float v0, float v1, float v2, vector<int>& numSigns) {
  int negative = 0;
  int positive = 0;
  int zero = 0;
  
  if (v0 > 0) positive++;
  if (v0 == 0) zero++;
  if (v0 < 0) negative++;
  if (v1 > 0) positive++;
  if (v1 == 0) zero++;
  if (v1 < 0) negative++;
  if (v2 > 0) positive++;
  if (v2 == 0) zero++;
  if (v2 < 0) negative++;

  numSigns[0] = positive;
  numSigns[1] = zero;
  numSigns[2] = negative;
}

bool isNegative(vector<int> numSigns) {
  return (numSigns[2] == 3 || (numSigns[2] == 1 && numSigns[1] == 2) || (numSigns[2] == 2 && numSigns[1] == 1));
}

bool isBoundary(vector<int> numSigns) {
  int positive = numSigns[0];
  int zero = numSigns[1];

  if ((positive == 1 && zero == 2) || (positive == 2 && zero == 1) || (zero == 3)) return true;
  else return false;
}

void organiseData(vec4 point, float v, vector<vec4>& in, vector<vec4>& boundary, vector<vec4>& out) {
  if (v > 0) in.push_back(point); 
  if (v == 0) boundary.push_back(point); 
  if (v < 0) out.push_back(point); 
}

void VertexShader(const vec4& v, Pixel& p, Vertex& vertex){
  vertex.position = getRotation() * (v - (camera.position + camera.basis[3]));
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
        camera.basis[3][2] -= 0.1;
        break;
      case SDLK_LEFT:
      // Move camera to the left
        camera.basis[3][0] -= 0.1;
        break;
      case SDLK_RIGHT:
      // Move camera to the right
        camera.basis[3][0] += 0.1;
        break;
      case SDLK_n:
        camera.basis[3][1] -= 0.1;
        break;
      case SDLK_m:
        camera.basis[3][1] += 0.1;
        break;
      case SDLK_d:
        // Rotate camera right;
        camera.basis =  generateRotation(vec3(0, -yaw, 0)) * camera.basis;
        if (is_lookAt) camera.position =  generateRotation(vec3(0, -yaw, 0)) * camera.position;
        break;
      case SDLK_a:
        // Rotate camera left;
        camera.basis =  generateRotation(vec3(0, yaw, 0)) * camera.basis;
        if (is_lookAt) camera.position =  generateRotation(vec3(0, yaw, 0)) * camera.position;

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