#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <limits.h>
#include "mpi.h"

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

#define N_DIMENSION 4
#define MASTER 0
#define SCREEN_WIDTH 500
#define SCREEN_HEIGHT 500
#define FULLSCREEN_MODE false
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
  vec3 circleNormal;
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

// vec4 cameraPos(0, 0, -2, 1.0);
mat4 R;




/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update(Camera &camera, Light &light, bool &escape, bool &is_lookAt);
void Draw(screen* screen, const vec3 *pixel_light_value, int local_ncols, int local_nrows, const vec3 *pixel_color_value, int local_cols, int local_rows, const int& row_start, const int& row_end, const int& col_start, const int& col_end);
bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closestIntersection);
bool circleIntersection(vec4 start, vec4 dir, const vector<Circle>& circles, Intersection& closestIntersection, int within_index = -1, bool isInCircle = false);
float calA(float radius);
vec3 calB(vec3 power, float radius);
vec3 calD(vec3 r, vec3 n, vec3 power, float radius);
vec3 calNormal(const Intersection& i, const vector<Circle>& circles);
vec3 DirectLight(Light &light, const Intersection& i, const vector<Triangle>& triangles);
vec3 circleDirectLight(Light &light, const Intersection& i, const vector<Circle>& circle);
vec3 circleRefract(const vec4& d, const Intersection& incident, const float& indexOfRefraction);
void processPart(Camera &camera, Light &light, vec3 *pixel_light_value, int local_ncols, int local_nrows, vec3 *pixel_color_value, int local_cols, int local_rows, const vector<Triangle> & triangles, vector<Circle>& circles, const int& row_start, const int& row_end, const int& col_start, const int& col_end);
void updateValues(int col, int row, vec3* pixel_light_value, int local_ncols, vec3 *pixel_color_value, int local_cols, int local_rows, const int& row_start, const int& col_start, const vec3& color, const vec3 light_power);
void circleRecursion(Light& light, const vector<Triangle>& triangles, const vector<Circle>& circles, const vec4& start_position, const vec4& new_direction, vec3 &light_power, vec3& color, const int& within_index, const float& ior, bool isInCircle);



int main( int argc, char* argv[] )
{

  int rank;
  int size;

  int local_nrows;
  int local_ncols;
  int loop_row_start_point;
  int loop_col_start_point;
  int loop_row_end_point;
  int loop_col_end_point;
  int tag = 0;
  MPI_Status status;

  bool escape = false;
  bool is_lookAt = false;

  float *sendbuf;
  float *recvbuf;

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

  vector<Triangle> triangles;
  vector<Circle> circles;
  LoadTestModel(triangles);
  LoadCircles(circles);

  MPI_Init(&argc, &argv);

  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if (size != N_DIMENSION) {
    fprintf(stderr, "Error: number of dimension is assumed to be 4\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  local_nrows = SCREEN_HEIGHT / 2;
  local_ncols = SCREEN_WIDTH / 2;

  vec3 color(0.0, 0.0, 0.0);

  vec3 pixel_light_value[local_nrows][local_ncols];
  vec3 pixel_color_value[local_nrows][local_ncols];

  switch (rank) {
    case 0:
      loop_row_start_point = 0;
      loop_col_start_point = 0;
      loop_row_end_point = local_nrows;
      loop_col_end_point = local_ncols;
      break;
    case 1:
      loop_row_start_point = 0;
      loop_col_start_point = local_ncols;
      loop_row_end_point = local_nrows;
      loop_col_end_point = SCREEN_WIDTH;
      break;
    case 2:
      loop_row_start_point = local_nrows;
      loop_col_start_point = 0;
      loop_row_end_point = SCREEN_HEIGHT;
      loop_col_end_point = local_ncols;
      break;
    case 3:
      loop_row_start_point = local_nrows;
      loop_col_start_point = local_ncols;
      loop_row_end_point = SCREEN_HEIGHT;
      loop_col_end_point = SCREEN_WIDTH;
      break;
  }

  sendbuf = (float*) malloc(sizeof(float) * (SCREEN_WIDTH * 3));
  recvbuf = (float*) malloc(sizeof(float) * (SCREEN_HEIGHT * 3));
  MPI_Buffer_attach(sendbuf, sizeof(float) * (SCREEN_WIDTH * 3));
  MPI_Buffer_attach(recvbuf, sizeof(float) * (SCREEN_HEIGHT * 3));

  bool is_screen = false;
  screen* screen;

  while( !escape ) {

    for (int i = 0; i < local_ncols; i++) {
      for (int j = 0; j < local_nrows; j++) {
        pixel_light_value[i][j] = color;
        pixel_color_value[i][j] = color;
      }
    }

    if (rank == MASTER) {
      if (!is_screen) {
        screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
        memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
        is_screen = true;
      }

      Update(camera, light, escape, is_lookAt);

      for (int i = 0; i < 4; i++) {
        sendbuf[i] = camera.position[i];
      }

      vec4 c0 = camera.basis[0];
      vec4 c1 = camera.basis[1];
      vec4 c2 = camera.basis[2];
      vec4 c3 = camera.basis[3];

      for (int i = 4; i < 8; i++) {
        sendbuf[i] = c0[i - 4];
      }

      for (int i = 8; i < 12; i++) {
        sendbuf[i] = c1[i - 8];
      }

      for (int i = 12; i < 16; i++) {
        sendbuf[i] = c2[i - 12];
      }

      for (int i = 16; i < 20; i++) {
        sendbuf[i] = c3[i-16];
      }

      for (int i = 20; i < 24; i++) {
        sendbuf[i] = light.position[i - 20];
      }

      if (escape) sendbuf[24] = 1;
      else sendbuf[24] = 0;

      if (escape) {
        SDL_SaveImage( screen, "screenshot.bmp" );
        KillSDL(screen);
        MPI_Abort(MPI_COMM_WORLD, 0);
        exit(0);
      }

      for (int k = 1; k < size; k++) {
        MPI_Send(sendbuf, 24, MPI_FLOAT, k, tag, MPI_COMM_WORLD);
      }

      processPart(camera, light, &pixel_light_value[0][0], local_ncols, local_nrows, &pixel_color_value[0][0], local_ncols, local_nrows, triangles, circles, loop_row_start_point, loop_row_end_point, loop_col_start_point, loop_col_end_point);

      Draw(screen, &pixel_light_value[0][0], local_ncols, local_nrows, &pixel_color_value[0][0], local_ncols, local_nrows, loop_row_start_point, loop_row_end_point, loop_col_start_point, loop_col_end_point);

      for (int i = 0; i < local_ncols; i++) {
        for (int j = 0; j < local_nrows; j++) {
          pixel_light_value[i][j] = color;
          pixel_color_value[i][j] = color;
        }
      }


      int rank_row_start;
      int rank_row_end;
      int rank_col_start;
      int rank_col_end;

      for (int k = 1; k < size; k++) {
        MPI_Recv(recvbuf, 4, MPI_FLOAT, k, tag, MPI_COMM_WORLD, &status);
        rank_row_start = recvbuf[0];
        rank_row_end = recvbuf[1];
        rank_col_start = recvbuf[2];
        rank_col_end = recvbuf[3];

        for (int i = 0; i < local_nrows; i++) {
          MPI_Recv(recvbuf, (local_nrows * 3), MPI_FLOAT, k, tag, MPI_COMM_WORLD, &status);
          for (int j = 0; j < local_ncols; j++) {
            pixel_light_value[i][j][0] = recvbuf[0 + 3 * j];
            pixel_light_value[i][j][1] = recvbuf[1 + 3 * j];
            pixel_light_value[i][j][2] = recvbuf[2 + 3 * j];
          }
        }
        for (int i = 0; i < local_nrows; i++) {
          MPI_Recv(recvbuf, (local_nrows * 3), MPI_FLOAT, k, tag, MPI_COMM_WORLD, &status);
          for (int j = 0; j < local_ncols; j++) {
            pixel_color_value[i][j][0] = recvbuf[0 + 3 * j];
            pixel_color_value[i][j][1] = recvbuf[1 + 3 * j];
            pixel_color_value[i][j][2] = recvbuf[2 + 3 * j];
          }
        }

        Draw(screen, &pixel_light_value[0][0], local_ncols, local_nrows, &pixel_color_value[0][0], local_ncols, local_nrows, rank_row_start, rank_row_end, rank_col_start, rank_col_end);

      }
      SDL_Renderframe(screen);
    } else {

      MPI_Recv(recvbuf, 24, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD, &status);

      camera.position = vec4(recvbuf[0], recvbuf[1], recvbuf[2], recvbuf[3]);

      vec4 c0(recvbuf[4], recvbuf[5], recvbuf[6], recvbuf[7]);
      vec4 c1(recvbuf[8], recvbuf[9], recvbuf[10], recvbuf[11]);
      vec4 c2(recvbuf[12], recvbuf[13], recvbuf[14], recvbuf[15]);
      vec4 c3(recvbuf[16], recvbuf[17], recvbuf[18], recvbuf[19]);

      camera.basis = mat4(c0, c1, c2, c3);

      light.position = vec4(recvbuf[20], recvbuf[21], recvbuf[22], recvbuf[23]);

      if (recvbuf[24] == 1) escape = true;
      else escape = false;

      processPart(camera, light, &pixel_light_value[0][0], local_ncols, local_nrows, &pixel_color_value[0][0], local_ncols, local_nrows, triangles, circles, loop_row_start_point, loop_row_end_point, loop_col_start_point, loop_col_end_point);

      sendbuf[0] = loop_row_start_point;
      sendbuf[1] = loop_row_end_point;
      sendbuf[2] = loop_col_start_point;
      sendbuf[3] = loop_col_end_point;

      MPI_Send(sendbuf, 4, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD);

      for (int i = 0; i < local_nrows; i++) {
        for (int j = 0; j < local_ncols; j++) {
            sendbuf[0 + 3 * j] = pixel_light_value[i][j][0];
            sendbuf[1 + 3 * j] = pixel_light_value[i][j][1];
            sendbuf[2 + 3 * j] = pixel_light_value[i][j][2];
        }
        MPI_Send(sendbuf, (local_ncols * 3), MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD);
      }
      for (int i = 0; i < local_nrows; i++) {
        for (int j = 0; j < local_ncols; j++) {
            sendbuf[0 + 3 * j] = pixel_color_value[i][j][0];
            sendbuf[1 + 3 * j] = pixel_color_value[i][j][1];
            sendbuf[2 + 3 * j] = pixel_color_value[i][j][2];
        }
        MPI_Send(sendbuf, (local_ncols * 3), MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD);
      }
    }
  }
  //
  // while( !escape )
  //   {
  //     Update();
  //     Draw(screen, triangles, loop_row_start_point, loop_row_end_point, loop_col_start_point, loop_col_end_point, rank);
  //     SDL_Renderframe(screen);
  //   }
  //
  SDL_SaveImage( screen, "screenshot.bmp" );
  KillSDL(screen);

  MPI_Finalize();
  exit(0);
  return 0;
}

/*Place your drawing here*/
void Draw(screen* screen, const vec3 *pixel_light_value, int local_ncols, int local_nrows, const vec3 *pixel_color_value, int local_cols, int local_rows, const int& row_start, const int& row_end, const int& col_start, const int& col_end)
{
  /* Clear buffer */

  vec3 colour(1.0,0.0,0.0);

  for (int row=row_start; row < row_end; row++){
    for (int col = col_start; col<col_end; col++ ){
      PutPixelSDL(screen, row, col, pixel_color_value[(col - col_start) + (row - row_start) * local_ncols] * pixel_light_value[(col - col_start) + (row - row_start) * local_ncols]);

      // vec4 d = camera.basis * vec4(row - SCREEN_WIDTH/2, col - SCREEN_HEIGHT/2, focal_length, 1);
      // Intersection intersect;
      // if (ClosestIntersection(camera.position, d, triangles, intersect)){
      //   vec3 light_power = DirectLight(intersect, triangles);
      //   PutPixelSDL(screen, row, col, triangles[intersect.triangleIndex].color * light_power);
      // }
    }
  }
}

void processPart(Camera &camera, Light &light, vec3 *pixel_light_value, int local_ncols, int local_nrows, vec3 *pixel_color_value, int local_cols, int local_rows, const vector<Triangle> & triangles, vector<Circle>& circles, const int& row_start, const int& row_end, const int& col_start, const int& col_end) {
  for (int row = row_start; row < row_end; row++) {
    for (int col = col_start; col < col_end; col++) {
      vec4 d = camera.basis * vec4(row - SCREEN_WIDTH/2, col - SCREEN_HEIGHT/2, focal_length, 1);
      Intersection intersect;
      if (ClosestIntersection(camera.position, d, triangles, intersect)){
        vec3 light_power = DirectLight(light, intersect, triangles);
        updateValues(col, row, pixel_light_value, local_ncols, pixel_color_value, local_cols, local_rows, row_start, col_start, triangles[intersect.triangleIndex].color, light_power);
      }
     
      if (circleIntersection(camera.position, d, circles, intersect)) {
        intersect.circleNormal = calNormal(intersect, circles);
        
        if (!circles[intersect.circleIndex].isGlass) {
          vec3 light_power = circleDirectLight(light, intersect, circles);
          updateValues(col, row, pixel_light_value, local_ncols, pixel_color_value, local_cols, local_rows, row_start, col_start, circles[intersect.circleIndex].color, light_power);
        
        } else {

          vec3 light_power = vec3(0,0,0);
          vec3 color = vec3(0,0,0);

          vec3 offset = vec3(0.001, 0.001, 0.001) * intersect.circleNormal;
          vec4 start_position = intersect.position - vec4(offset.x, offset.y, offset.z, 1);

          vec3 refract = circleRefract(d, intersect, 1.5);
          vec4 new_direction = vec4(refract.x, refract.y, refract.z, 1);     

          circleRecursion(light, triangles, circles, start_position, new_direction, light_power, color, intersect.circleIndex, 1.5, true);
          updateValues(col, row, pixel_light_value, local_ncols, pixel_color_value, local_cols, local_rows, row_start, col_start, color, light_power);

        }
      }

    }
  }
}

void circleRecursion(Light& light, const vector<Triangle>& triangles, const vector<Circle>& circles, const vec4& start_position, const vec4& new_direction, vec3 &light_power, vec3& color, const int& within_index, const float& ior, bool isInCircle) {
  
  Intersection intersection;


  if (isInCircle) {


    if (circleIntersection(start_position, new_direction, circles, intersection, within_index, true)) {
      intersection.circleNormal = calNormal(intersection, circles);

      vec3 offset = vec3(0.001, 0.001, 0.001) * intersection.circleNormal;

      vec4 next_start_position = intersection.position + vec4(offset.x, offset.y, offset.z, 1);

      vec3 refract = circleRefract(new_direction, intersection, 1.5);
      vec4 next_new_direction = vec4(refract.x, refract.y, refract.z, 1);
      circleRecursion(light, triangles, circles, next_start_position, next_new_direction, light_power, color, intersection.circleIndex, 1.5, false);
      
    } 
      
  } else {

    if (ClosestIntersection(start_position, new_direction, triangles, intersection)) {
        light_power = DirectLight(light, intersection, triangles);
        color = triangles[intersection.triangleIndex].color;
    } 
    if (circleIntersection(start_position, new_direction, circles, intersection, within_index)) {
      intersection.circleNormal = calNormal(intersection, circles);

      if (!circles[intersection.circleIndex].isGlass) {
        light_power = circleDirectLight(light, intersection, circles);
        color = circles[intersection.circleIndex].color;
      } else {
        vec3 offset = vec3(0.001, 0.001, 0.001) * intersection.circleNormal;
        vec4 next_start_position = intersection.position - vec4(offset.x, offset.y, offset.z, 1);

        vec3 refract = circleRefract(new_direction, intersection, 1.5);
        vec4 next_new_direction = vec4(refract.x, refract.y, refract.z, 1);

        circleRecursion(light, triangles, circles, next_start_position, next_new_direction, light_power, color, intersection.circleIndex, 1.5, true);
        
      }
    }
  }
}

void updateValues(int col, int row, vec3* pixel_light_value, int local_ncols, vec3* pixel_color_value, int local_cols, int local_rows, const int& row_start, const int& col_start, const vec3& color, const vec3 light_power) {
  pixel_light_value[(col - col_start) + (row - row_start) * local_ncols] = light_power;
  pixel_color_value[(col - col_start) + (row - row_start) * local_ncols] = color;
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

    vec3 L = vec3(L4.x, L4.y, L4.z);
    vec3 d3 = vec3(d.x, d.y, d.z);

    float a = dot(d3, d3);
    float b = 2 * dot(d3, L);
    float c = dot(L, L) - pow(circles[i].radius, 2);

    if ((pow(b,2) - 4 * a * c) < 0) continue;
    
    

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

vec3 calNormal(const Intersection& i, const vector<Circle>& circles) {
  return normalize(vec3(i.position - circles[i.circleIndex].center));
}

vec3 DirectLight(Light &light, const Intersection& i, const vector<Triangle>& triangles) {

  vec3 n(triangles[i.triangleIndex].normal.x, triangles[i.triangleIndex].normal.y, triangles[i.triangleIndex].normal.z);

  vec3 x = normalize(n);

  vec3 r = normalize(light.position - i.position);

  float radius = glm::distance(light.position, i.position);

  return calD(r, x, light.color, radius);
}

vec3 circleDirectLight(Light& light, const Intersection& i, const vector<Circle>& circles) {

  vec3 normal = i.circleNormal;
  // circles[i.circleIndex].normal = normal;

  vec3 r = normalize(light.position - i.position);

  float radius = glm::distance(light.position, i.position);

  return calD(r, normal, light.color, radius);
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
    std::swap(airMedium, materialMedium);
  }

  float mediumDiv = airMedium / materialMedium;

  float k = 1 - mediumDiv * mediumDiv * (1 - NdotI * NdotI);

  if (k < 0) return vec3(0,0,0);
  else return mediumDiv * i + (mediumDiv * NdotI - sqrt(k)) * normal;
  
}


void Update(Camera &camera, Light &light, bool &escape, bool &is_lookAt)
{
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;
  /*Good idea to remove this*/
  // std::cout << "Render time: " << dt << " ms." << std::endl;
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