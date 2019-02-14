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
#define SCREEN_WIDTH 50
#define SCREEN_HEIGHT 50
#define FULLSCREEN_MODE false
#define PI 3.14159

float maxFloat = std::numeric_limits<float>::max();

float yaw = 2 * PI / 180;

float focal_length = SCREEN_HEIGHT / 2;

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

// Camera camera = {
//   .position = vec4(0,0,-2, 1.0),
//   .basis = mat4(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), vec4(0,0,0,1)),
//   // .center = vec3(0.003724, 0.929729, 0.07459)
//   .center = vec3(0,0,0)
// };
//
// Light light = {
//   .position = vec4(0, -0.5, -0.7, 1.0),
//   .color = 14.f * vec3(1,1,1),
// };

// // vec4 cameraPos(0, 0, -2, 1.0);
// bool escape = false;
// bool is_lookAt = false;




/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update(Camera& camera, Light& light);
void Draw(screen* screen, const vector<vec3>& pixel_light_value, const vector<vec3>& pixel_color_value, const int& row_start, const int& row_end, const int& col_start, const int& col_end);
bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closestIntersection);
float calA(float radius);
vec3 calB(vec3 power, float radius);
vec3 calD(vec3 r, vec3 n, vec3 power, float radius);
vec3 DirectLight(vec4 lightPos, vec3 lightCol, const Intersection& i, const vector<Triangle>& triangles);
void processPart(vec4 cameraPos, vec4 lightPos, vec3 lightCol, vector<vec3>& pixel_light_value, vector<vec3>& pixel_color_value, const vector<Triangle>& triangles, const int& row_start, const int& row_end, const int& col_start, const int& col_end);




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

  float *sendbuf;
  float *recvbuf;

  vector<Triangle> triangles;
  LoadTestModel(triangles);

  MPI_Init(&argc, &argv);

  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  vector<vec3> pixel_light_value;
  vector<vec3> pixel_color_value;

  if (size != N_DIMENSION) {
    fprintf(stderr, "ErrorL number of dimension is assumed to be 4\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  local_nrows = SCREEN_HEIGHT / 2;
  local_ncols = SCREEN_WIDTH / 2;

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

  sendbuf = (float*) malloc(sizeof(float) * (local_ncols * 3));
  recvbuf = (float*) malloc(sizeof(float) * (local_ncols * 3));

  while( true ) {

    if (rank == MASTER) {
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

      Update(camera, light);

      sendbuf[0] = camera.position[0];
      sendbuf[1] = camera.position[1];
      sendbuf[2] = camera.position[2];
      sendbuf[3] = camera.position[3];
      sendbuf[4] = light.position[0];
      sendbuf[5] = light.position[1];
      sendbuf[6] = light.position[2];
      sendbuf[7] = light.position[3];
      sendbuf[8] = light.color[0];
      sendbuf[9] = light.color[1];
      sendbuf[10] = light.color[2];

      for (int k = 1; k < size; k++) {
        MPI_Send(sendbuf, 11, MPI_FLOAT, k, tag, MPI_COMM_WORLD);
      }

      screen* screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
      memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

      processPart(camera.position, light.position, light.color, pixel_light_value, pixel_color_value, triangles, loop_row_start_point, loop_row_end_point, loop_col_start_point, loop_col_end_point);
      Draw(screen, pixel_light_value, pixel_color_value, loop_row_start_point, loop_row_end_point, loop_col_start_point, loop_col_end_point);

      vector<vec3> rank_light;
      vector<vec3> rank_color;

      int rank_row_start;
      int rank_row_end;
      int rank_col_start;
      int rank_col_end;
      float x;
      float y;
      float z;
      for (int k = 1; k < size; k++) {
        MPI_Recv(recvbuf, 4, MPI_FLOAT, k, tag, MPI_COMM_WORLD, &status);

        rank_row_start = recvbuf[0];
        rank_row_end = recvbuf[1];
        rank_col_start = recvbuf[2];
        rank_col_end = recvbuf[3];

        for (int i = 0; i < local_nrows; i++) {
          MPI_Recv(recvbuf, (local_nrows * 3), MPI_FLOAT, k, tag, MPI_COMM_WORLD, &status);
          for (int j = 0; j < local_ncols; j++) {
            x = recvbuf[0 + 3 * j];
            y = recvbuf[1 + 3 * j];
            z = recvbuf[2 + 3 * j];

            rank_light[j + i * local_ncols] = vec3(x, y, z);
          }
        }
        for (int i = 0; i < local_nrows; i++) {
          MPI_Recv(recvbuf, (local_nrows * 3), MPI_FLOAT, k, tag, MPI_COMM_WORLD, &status);
          for (int j = 0; j < local_ncols; j++) {
            x = recvbuf[0 + 3 * j];
            y = recvbuf[1 + 3 * j];
            z = recvbuf[2 + 3 * j];

            rank_color[j + i * local_ncols] = vec3(x, y, z);
          }
        }

        Draw(screen, rank_light, rank_color, rank_row_start, rank_row_end, rank_col_start, rank_col_end);
      }
      SDL_Renderframe(screen);
    } else {

      MPI_Recv(recvbuf, 11, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD, &status);

      vec4 cameraPos(sendbuf[0], sendbuf[1], sendbuf[2], sendbuf[3]);
      vec4 lightPos(sendbuf[4], sendbuf[5], sendbuf[6], sendbuf[7]);
      vec3 lightCol(sendbuf[8], sendbuf[9], sendbuf[10]);

      processPart(cameraPos, lightPos, lightCol, pixel_light_value, pixel_color_value, triangles, loop_row_start_point, loop_row_end_point, loop_col_start_point, loop_col_end_point);

      sendbuf[0] = loop_row_start_point;
      sendbuf[1] = loop_row_end_point;
      sendbuf[2] = loop_col_start_point;
      sendbuf[3] = loop_col_end_point;

      MPI_Send(sendbuf, 4, MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD);

      for (int i = 0; i < local_nrows; i++) {
        for (int j = 0; j < local_ncols; j++) {
            sendbuf[0 + 3 * j] = pixel_light_value[j + i * local_ncols].x;
            sendbuf[1 + 3 * j] = pixel_light_value[j + i * local_ncols].y;
            sendbuf[2 + 3 * j] = pixel_light_value[j + i * local_ncols].z;
        }
        MPI_Send(sendbuf, (local_ncols * 3), MPI_FLOAT, MASTER, tag, MPI_COMM_WORLD);
      }
      for (int i = 0; i < local_nrows; i++) {
        for (int j = 0; j < local_ncols; j++) {
            sendbuf[0 + 3 * j] = pixel_color_value[j + i * local_ncols].x;
            sendbuf[1 + 3 * j] = pixel_color_value[j + i * local_ncols].y;
            sendbuf[2 + 3 * j] = pixel_color_value[j + i * local_ncols].z;
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
  // SDL_SaveImage( screen, "screenshot.bmp" );
  // KillSDL(screen);

  MPI_Finalize();
  return 0;
}

/*Place your drawing here*/
void Draw(screen* screen, const vector<vec3>& pixel_light_value, const vector<vec3>& pixel_color_value, const int& row_start, const int& row_end, const int& col_start, const int& col_end)
{
  /* Clear buffer */

  vec3 colour(1.0,0.0,0.0);

  for (int row=row_start; row < row_end; row++){
    for (int col = col_start; col<col_end; col++ ){
      PutPixelSDL(screen, row, col, pixel_color_value[(col - col_start) + (row - row_start) * (col_end - col_start)] * pixel_light_value[(col - col_start) + (row - row_start) * (col_end - col_start)]);
      // vec4 d = camera.basis * vec4(row - SCREEN_WIDTH/2, col - SCREEN_HEIGHT/2, focal_length, 1);
      // Intersection intersect;
      // if (ClosestIntersection(camera.position, d, triangles, intersect)){
      //   vec3 light_power = DirectLight(intersect, triangles);
      //   PutPixelSDL(screen, row, col, triangles[intersect.triangleIndex].color * light_power);
      // }
    }
  }

}

void processPart(vec4 cameraPos, vec4 lightPos, vec3 lightCol, vector<vec3>& pixel_light_value, vector<vec3>& pixel_color_value, const vector<Triangle> & triangles, const int& row_start, const int& row_end, const int& col_start, const int& col_end) {
  for (int row = row_start; row < row_end; row++) {
    for (int col = col_start; col < col_end; col++) {
      // vec4 d = camera.basis * vec4(row - SCREEN_WIDTH/2, col - SCREEN_HEIGHT/2, focal_length, 1);
      vec4 d = vec4(row - SCREEN_WIDTH/2, col - SCREEN_HEIGHT/2, focal_length, 1);
      Intersection intersect;
      if (ClosestIntersection(cameraPos, d, triangles, intersect)){
        vec3 light_power = DirectLight(lightPos, lightCol, intersect, triangles);
        // pixel_light_value[(col - col_start) + (row - row_start) * (col_end - col_start)] = light_power;
        printf("\nI'm here\n");
        // pixel_color_value[(col - col_start) + (row - row_start) * (col_end - col_start)] = triangles[intersect.triangleIndex].color;
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
vec3 DirectLight(vec4 lightPos, vec3 lightCol, const Intersection& i, const vector<Triangle>& triangles) {

  vec3 n(triangles[i.triangleIndex].normal.x, triangles[i.triangleIndex].normal.y, triangles[i.triangleIndex].normal.z);

  vec3 x = normalize(n);

  vec3 r = normalize(lightPos - i.position);

  float radius = glm::distance(lightPos, i.position);

  return calD(r, x, lightCol, radius);
}

void Update(Camera& camera, Light& light)
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
        // escape = true;
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
        // if (is_lookAt) camera.position = translation * generateRotation(vec3(0, yaw, 0)) * camera.position;
        break;
      case SDLK_a:
        // Rotate camera left;
        camera.basis = translation * generateRotation(vec3(0, -yaw, 0)) * camera.basis;
        // if (is_lookAt) camera.position = translation * generateRotation(vec3(0, -yaw, 0)) * camera.position;

        break;
      case SDLK_w:
        // Rotate camera top;
        camera.basis = translation * generateRotation(vec3(yaw, 0, 0)) * camera.basis;
        // if (is_lookAt) camera.position = translation * generateRotation(vec3(yaw, 0, 0)) * camera.position;

        break;
      case SDLK_s:
        // Rotate camera down;
        camera.basis = translation * generateRotation(vec3(-yaw, 0, 0)) * camera.basis;
        // if (is_lookAt) camera.position = translation * generateRotation(vec3(-yaw, 0, 0)) * camera.position;
        break;
      case SDLK_q:
        camera.basis = translation * generateRotation(vec3(0, 0, -yaw)) * camera.basis;
        // if (is_lookAt) camera.position = translation * generateRotation(vec3(0, 0, -yaw)) * camera.position;
        break;
      case SDLK_e:
        camera.basis = translation * generateRotation(vec3(0, 0, yaw)) * camera.basis;
        // if (is_lookAt) camera.position = translation * generateRotation(vec3(0, 0, yaw)) * camera.position;
        break;
      case SDLK_l:
        // if (is_lookAt) is_lookAt = false;
        // else is_lookAt = true;
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

    // if (is_lookAt) {
      // vec3 position = vec3(camera.position[0], camera.position[1], camera.position[2]);
      // camera.basis = lookAt(camera.center, position);
    // }
   }
 }
}
