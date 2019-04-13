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


#define SCREEN_WIDTH 1000
#define SCREEN_HEIGHT 1000
#define FULLSCREEN_MODE false
#define PI 3.14159

float maxFloat = std::numeric_limits<float>::max();

float yaw = 2 * PI / 180;

const float bias = 0.0085f;

const int antiFactor = 3;

const int antiWidth = SCREEN_WIDTH * antiFactor;
const int antiHeight= SCREEN_HEIGHT * antiFactor;

float focal_length = antiWidth / 2;

float depth = 0.55;

const float epsilon = 0.1;


struct Camera{
  vec4 position;
  mat4 basis;
  mat4 planeBasis;
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
  vec3 direction;
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
  .position = vec4(0,0,-2, 1.0),
  .basis = mat4(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), vec4(0,0,0,1)),
  .planeBasis = mat4(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), vec4(0,0,0,1)),
  // .center = vec3(0.003724, 0.929729, 0.07459)
  .center = vec3(0,0,0)
};

Light light = {
  .position = vec4(0, -0.9, -0.9, 1.0),
  .color = 14.0f * vec3(1, 1, 1),
  .indirectLight = 0.7f * vec3(1, 1, 1),
  .direction = vec3((90 * PI) / 180, 0, 0)
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


const mat4 identity = mat4(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), vec4(0,0,0,1));

float depthBuffer[antiHeight][antiWidth];
float lightDepthBuffer[antiHeight][antiWidth];
float lightDepthUpBuffer[antiHeight][antiWidth];
vec3 pixelValue[antiHeight][antiWidth];


// vec4 cameraPos(0, 0, -2, 1.0);
bool escape = false;
bool is_lookAt = false;

vec4 preLightPos = light.position;


/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Render(const vector<Triangle>& triangles);
void Draw(screen* screen);
void VertexShader(const vec4& v, Pixel& p, Vertex& vertex, const mat4 basis, const vec4 position, const float focal_length);
// void DrawLineSDL(screen* surface, ivec2 a, ivec2 b, vec3 color);
// void DrawPolygonEdges(screen* screen, const vector<vec4>& vertices);
void ComputePolygonRows(const vector<Pixel>& vertextPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
void DrawPolygon( const vector<vec4>& vertices, vector<Vertex>& vertex);
void PixelShader( const Pixel& p);
mat4 generateRotation(vec3 a);
void ClipTriangles(const vector<Triangle>& triangles, vector<Triangle>& clippedTriangles, vec4 normal, vec4 point);
void calSign(vec4 trianglePoint, float& value, vec4 normal, vec4 point);
void countSign(float v0, float v1, float v2, vector<int>& numSigns);
bool isBoundary(vector<int> numSigns);
bool isNegative(vector<int> numSigns);
void organiseData(vec4 point, float v, vector<vec4>& in, vector<vec4>& boundary, vector<vec4>& out);
void updateClippers(const mat4 basis, const vec4 position);
void fillLightBuffer(const vector<Triangle> triangles, const mat4 lightBasis, const mat4 clipBasis, const vec4 position, bool isUp = false);
void lightClip(const vector<Triangle> triangles, vector<Triangle>& finalTriangles);
void renderFromLight(const vector<Triangle>& triangles, const mat4 lightBasis, const vec4 position, bool isUp = false);
void updateBuffer(const vector<Pixel> leftPixels, const vector<Pixel> rightPixels, bool isUp = false);
void lightShader(const Pixel &p, bool isUp = false);
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result);

int main( int argc, char* argv[] )
{

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  vector<Triangle> triangles;
  LoadTestModel(triangles);

  int size = triangles.size();

  mat4 lightBasis = generateRotation(light.direction) * identity;
  mat4 clipBasis = generateRotation(-light.direction) * identity;

  fillLightBuffer(triangles, lightBasis, clipBasis, light.position);
  fillLightBuffer(triangles, clipBasis, lightBasis, light.position, true);

  while( !escape )
    {
      vector<Triangle> nearTriangles;
      vector<Triangle> leftTriangles;
      vector<Triangle> rightTriangles;
      vector<Triangle> topTriangles;
      vector<Triangle> bottomTriangles;
      Update();
      if (preLightPos != light.position) {
        preLightPos = light.position;
        fillLightBuffer(triangles, lightBasis, clipBasis, light.position);
        fillLightBuffer(triangles, clipBasis, lightBasis, light.position, true);
      }

      updateClippers(camera.planeBasis, camera.position);
      ClipTriangles(triangles, nearTriangles, clipper.nearNormal, clipper.nearPoint);
      ClipTriangles(nearTriangles, leftTriangles, clipper.leftNormal, clipper.leftPoint);
      ClipTriangles(leftTriangles, rightTriangles, clipper.rightNormal, clipper.rightPoint);
      ClipTriangles(rightTriangles, topTriangles, clipper.topNormal, clipper.topPoint);
      ClipTriangles(topTriangles, bottomTriangles, clipper.botNormal, clipper.botPoint);
      Render(bottomTriangles);
      Draw(screen);
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );
  KillSDL(screen);
  return 0;
}

void Draw(screen* screen) {

  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  
  for (int x = 0; x < SCREEN_WIDTH; x++) {
    for (int y = 0; y < SCREEN_HEIGHT; y++) {

        vec3 value = vec3(0,0,0);
        float addNumber = 0.f;

        int currentX = x * antiFactor;
        int currentY = y * antiFactor;
        int nextX = (x + 1) * antiFactor;
        int nextY = (y + 1) * antiFactor;

        for (int w = currentX; w < nextX; w++) {
          for (int v = currentY; v < nextY; v++) {
            value += pixelValue[v][w];
            addNumber += 1;
          }
        }

        value = value / addNumber;

        PutPixelSDL(screen, x, y, value);
    }
  }
}

void fillLightBuffer(const vector<Triangle> triangles, const mat4 lightBasis, const mat4 clipBasis, const vec4 position, bool isUp) {
  vector<Triangle> finalTriangles;
  updateClippers(clipBasis, position);
  lightClip(triangles, finalTriangles);
  renderFromLight(finalTriangles, lightBasis, position, isUp);
}

void lightClip(const vector<Triangle> triangles, vector<Triangle>& finalTriangles) {
  vector<Triangle> nearTriangles;
  vector<Triangle> leftTriangles;
  vector<Triangle> rightTriangles;
  vector<Triangle> topTriangles;

  ClipTriangles(triangles, nearTriangles, clipper.nearNormal, clipper.nearPoint);
  ClipTriangles(nearTriangles, leftTriangles, clipper.leftNormal, clipper.leftPoint);
  ClipTriangles(leftTriangles, rightTriangles, clipper.rightNormal, clipper.rightPoint);
  ClipTriangles(rightTriangles, topTriangles, clipper.topNormal, clipper.topPoint);
  ClipTriangles(topTriangles, finalTriangles, clipper.botNormal, clipper.botPoint);

}


void renderFromLight(const vector<Triangle>& triangles, const mat4 lightBasis, const vec4 position, bool isUp) {
  for (int y = 0; y < antiHeight; y++) {
    for (int x = 0; x < antiWidth; x++) {
      if (!isUp) lightDepthBuffer[y][x] = 0.0f;
      else lightDepthUpBuffer[y][x] = 0.0f;
    }
  }

  for (uint32_t i=0; i<triangles.size(); ++i){
    vector<vec4> vertices(3);
    vector<Vertex> vertex(3);

    vertices[0] = triangles[i].v0;
    vertices[1] = triangles[i].v1;
    vertices[2] = triangles[i].v2;

    vector<Pixel> vertexPixels(3);

    for (int i = 0; i < 3; i++) {
      VertexShader(vertices[i], vertexPixels[i], vertex[i], lightBasis, position, focal_length / 2);
    }

    vector<Pixel> leftPixels;
    vector<Pixel> rightPixels;
    ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
    updateBuffer(leftPixels, rightPixels, isUp);
  }
}

void updateBuffer(const vector<Pixel> leftPixels, const vector<Pixel> rightPixels, bool isUp) {
  
  for (uint i = 0; i < leftPixels.size(); i++) {
    int left = leftPixels[i].x;
    int right = rightPixels[i].x;

    Pixel depth_left = leftPixels[i];
    Pixel depth_right = rightPixels[i];

    int distance = abs(right - left) + 1;
    vector<Pixel> drawRow(distance);

    Interpolate(depth_left, depth_right, drawRow);
    int k = 0;
    for (int j = left; j < right; j++) {
      if (j >= 0 && j < antiWidth) {
        if (leftPixels[i].y >= 0 && leftPixels[i].y < antiHeight) {
          lightShader(drawRow[k], isUp);
        }
      }
      k++;
    }
  }
}

void lightShader(const Pixel &p, bool isUp) {
  int x = p.x;
  int y = p.y;

  if (p.zinv >= 0) {
    if (!isUp) {
      if (lightDepthBuffer[y][x] < p.zinv) {
        lightDepthBuffer[y][x] = p.zinv;
      }
    } else {
      if (lightDepthUpBuffer[y][x] < p.zinv) {
        lightDepthUpBuffer[y][x] = p.zinv;
      }
    }
  }
}

/*Place your drawing here*/
void Render(const vector<Triangle>& triangles)
{
  
  for (int y = 0; y < antiHeight; y++) {
    for (int x = 0; x < antiWidth; x++) {
      depthBuffer[y][x] = 0.0f;
      pixelValue[y][x] = vec3(0,0,0);
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
    // current.currentReflectance = vec3(1,1,1);

    DrawPolygon(vertices, vertex);
    // DrawPolygonEdges(screen, vertices);

  }
}

void updateClippers(const mat4 basis, const vec4 position) {

    // These are the directions towards the four corners of the img plane
    
    vec4 leftUpCorner = normalize(basis * vec4(-antiWidth/2, -antiHeight/2, focal_length / 2, 1));
    vec3 leftUp = vec3(leftUpCorner);
   
    vec4 leftBotCorner = normalize(basis * vec4(-antiWidth/2, antiHeight/2, focal_length / 2, 1));
    vec3 leftBot = vec3(leftBotCorner);

    vec4 rightTopCorner = normalize(basis * vec4(antiWidth/2, -antiHeight/2, focal_length / 2, 1));
    vec3 rightUp = vec3(rightTopCorner);

    vec4 rightBotCorner = normalize(basis * vec4(antiWidth/2, antiHeight/2, focal_length / 2, 1));
    vec3 rightBot = vec3(rightBotCorner);

    // We use those directions to get the normal of the plane
    vec3 leftNormal = (glm::cross(leftUp, -leftBot));
    vec3 rightNormal =(glm::cross(rightUp, rightBot));
    vec3 topNormal = (glm::cross(leftUp, rightUp));
    vec3 botNormal = (glm::cross(leftBot, -rightBot));

    
    
    vec4 leftN = vec4(leftNormal.x, leftNormal.y, leftNormal.z, 1);
    vec4 rightN = vec4(rightNormal.x, rightNormal.y, rightNormal.z, 1);
    vec4 topN = vec4(topNormal.x, topNormal.y, topNormal.z, 1);
    vec4 botN = vec4(botNormal.x, botNormal.y, botNormal.z, 1);
    // vec4 temp = leftN;
    // leftN = vec4(-rightN.x, rightN.y, rightN.z, 1);
    // rightN = vec4(-temp.x, temp.y, temp.z, 1);

    // We use the directions from the top to find points in the near plane and then find the normal of this plane
    vec4 leftTopPoint = position + (leftUpCorner * depth);
    vec4 leftBottomPoint = position + (leftBotCorner * depth);
    vec4 rightBottomPoint = position + (rightBotCorner * depth);
    vec4 rightTopPoint = position + (rightTopCorner * depth);
    vec4 leftToRight = rightBottomPoint - leftBottomPoint;
    vec4 leftToTop = leftTopPoint - leftBottomPoint;
    vec4 topTotop = rightTopPoint - leftTopPoint;

    vec3 leftToRight3 = vec3(leftToRight);
    vec3 leftToTop3 = vec3(leftToTop);

    vec3 nearNormal = (glm::cross(-leftToRight3, leftToTop3));
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
    // std::cout<<glm::to_string(in[0])<<std::endl;
    
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

void VertexShader(const vec4& v, Pixel& p, Vertex& vertex, const mat4 basis, const vec4 position, const float focal_length){
  vertex.position = basis * (v - position);
  p.zinv = 1/vertex.position[2];
  p.x = focal_length*(vertex.position[0]/vertex.position[2]) + antiWidth/2;
  p.y = focal_length*(vertex.position[1]/vertex.position[2]) + antiHeight/2;
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

void DrawRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels) {

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
      if (j >= 0 && j < antiWidth) {
        if (leftPixels[i].y >= 0 && leftPixels[i].y < antiHeight) {
          PixelShader(drawRow[k]);
        }
      }
      k++;
    }
  }
}

void PixelShader(const Pixel& p) {

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

      vec3 illumination;
      
      vec4 fromLight = (generateRotation(light.direction) * identity) * (p.pos3d - light.position);

      float zinv = 1/fromLight[2] + bias;
      int v = focal_length / 2 * (fromLight[0] / fromLight[2]) + antiWidth/2;
      int w = focal_length / 2 * (fromLight[1] / fromLight[2]) + antiHeight/2;
      
      if (v < 0 || v > antiWidth) illumination = ((light.color * max((r_n), 0.f))/area) + light.indirectLight;
      else if (w < 0 || w > antiHeight) illumination = ((light.color * max((r_n), 0.f))/area) + light.indirectLight;
      else if (lightDepthBuffer[w][v] - bias < zinv || lightDepthUpBuffer[w][v] - bias > zinv) illumination = ((light.color * max((r_n), 0.f))/area) + light.indirectLight;
      else if (lightDepthBuffer[w][v] - bias > zinv || lightDepthUpBuffer[w][v] - bias < zinv) illumination = light.indirectLight;
      // else if () illumination = light.indirectLight;
      else illumination = light.indirectLight;

      pixelValue[y][x] = current.currentReflectance * illumination;

    }
  }
}

void DrawPolygon(const vector<vec4>& vertices, vector<Vertex>& vertex) {
  int V = vertices.size();

  vector<Pixel> vertexPixels(V);

  for (int i = 0; i < V; i++) {
    VertexShader(vertices[i], vertexPixels[i], vertex[i], camera.basis, camera.position, focal_length);
  }

  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;
  ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
  DrawRows(leftPixels, rightPixels);
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
        camera.planeBasis = translation * generateRotation(vec3(0, -yaw, 0)) * camera.planeBasis;
        if (is_lookAt) camera.position = translation * generateRotation(vec3(0, yaw, 0)) * camera.position;
        break;
      case SDLK_a:
        // Rotate camera left;
        camera.basis = translation * generateRotation(vec3(0, -yaw, 0)) * camera.basis;
        camera.planeBasis = translation * generateRotation(vec3(0, yaw, 0)) * camera.planeBasis;
        if (is_lookAt) camera.position = translation * generateRotation(vec3(0, -yaw, 0)) * camera.position;

        break;
      case SDLK_w:
        // Rotate camera top;
        camera.basis = translation * generateRotation(vec3(yaw, 0, 0)) * camera.basis;
        camera.planeBasis = translation * generateRotation(vec3(-yaw, 0, 0)) * camera.planeBasis;
        if (is_lookAt) camera.position = translation * generateRotation(vec3(yaw, 0, 0)) * camera.position;

        break;
      case SDLK_s:
        // Rotate camera down;
        camera.basis = translation * generateRotation(vec3(-yaw, 0, 0)) * camera.basis;
        camera.planeBasis = translation * generateRotation(vec3(yaw, 0, 0)) * camera.planeBasis;
        if (is_lookAt) camera.position = translation * generateRotation(vec3(-yaw, 0, 0)) * camera.position;
        break;
      case SDLK_q:
        camera.basis = translation * generateRotation(vec3(0, 0, -yaw)) * camera.basis;
        camera.planeBasis = translation * generateRotation(vec3(0, 0, yaw)) * camera.planeBasis;
        if (is_lookAt) camera.position = translation * generateRotation(vec3(0, 0, -yaw)) * camera.position;
        break;
      case SDLK_e:
        camera.basis = translation * generateRotation(vec3(0, 0, yaw)) * camera.basis;
        camera.planeBasis = translation * generateRotation(vec3(0, 0, -yaw)) * camera.planeBasis;
        if (is_lookAt) camera.position = translation * generateRotation(vec3(0, 0, yaw)) * camera.position;
        break;
      case SDLK_l:
        if (is_lookAt) is_lookAt = false;
        else is_lookAt = true;
        break;
      case SDLK_u:
        light.position += vec4(0, 0, 0.1, 0);
        // depth += 0.1;
        break;
      case SDLK_j:
        light.position += vec4(0,0,-0.1,0);
        // depth -= 0.1;
        break;
      case SDLK_h:
        light.position += vec4(-0.1, 0, 0, 0);
        break;
      case SDLK_k:
        light.position += vec4(0.1, 0, 0, 0);
        break;
      case SDLK_y:
        light.position += vec4(0,0.1,0,0);
        break;
      case SDLK_i:
        light.position += vec4(0,-0.1,0,0);
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