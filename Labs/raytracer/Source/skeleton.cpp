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
#define SCREEN_WIDTH 200
#define SCREEN_HEIGHT 200
#define FULLSCREEN_MODE true
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

struct Options{
  float bias;
  vec3 indirectLight;
};

// vec4 cameraPos(0, 0, -2, 1.0);
mat4 R;




/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update(Camera &camera, Light &light, bool &escape, bool &is_lookAt);
void Draw(screen* screen, const vec3 *pixel_light_value, int local_ncols, int local_nrows, const vec3 *pixel_color_value, int local_cols, int local_rows, const int& row_start, const int& row_end, const int& col_start, const int& col_end);
bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closestIntersection);
float calA(float radius);
vec3 calB(vec3 power, float radius);
vec3 calD(vec3 r, vec3 n, vec3 power, float radius);
vec3 DirectLight(Light &light, Options &options, const Intersection& i, const vector<Triangle>& triangles);
void processPart(Camera &camera, Light &light, Options &options, vec3 *pixel_light_value, int local_ncols, int local_nrows, vec3 *pixel_color_value, int local_cols, int local_rows, const vector<Triangle> & triangles, const int& row_start, const int& row_end, const int& col_start, const int& col_end);




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

  Options options = {
    .bias = 1e-2,
    .indirectLight = 0.5f * vec3(1,1,1)
  };

  vector<Triangle> triangles;
  LoadTestModel(triangles);

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

      processPart(camera, light, options, &pixel_light_value[0][0], local_ncols, local_nrows, &pixel_color_value[0][0], local_ncols, local_nrows, triangles, loop_row_start_point, loop_row_end_point, loop_col_start_point, loop_col_end_point);

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

      processPart(camera, light, options, &pixel_light_value[0][0], local_ncols, local_nrows, &pixel_color_value[0][0], local_ncols, local_nrows, triangles, loop_row_start_point, loop_row_end_point, loop_col_start_point, loop_col_end_point);

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

void processPart(Camera &camera, Light &light, Options &options, vec3 *pixel_light_value, int local_ncols, int local_nrows, vec3 *pixel_color_value, int local_cols, int local_rows, const vector<Triangle> & triangles, const int& row_start, const int& row_end, const int& col_start, const int& col_end) {
  for (int row = row_start; row < row_end; row++) {
    for (int col = col_start; col < col_end; col++) {
      vec4 d = camera.basis * vec4(row - SCREEN_WIDTH/2, col - SCREEN_HEIGHT/2, focal_length, 1);
      Intersection intersect;
      if (ClosestIntersection(camera.position, d, triangles, intersect)){
        vec3 light_power = DirectLight(light, options, intersect, triangles);
        pixel_light_value[(col - col_start) + (row - row_start) * local_ncols] = light_power;
        pixel_color_value[(col - col_start) + (row - row_start) * local_ncols] = triangles[intersect.triangleIndex].color;
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

vec3 DirectLight(Light &light, Options &options, const Intersection& i, const vector<Triangle>& triangles){
  float r = glm::distance(light.position, i.position);
  float area = 4 * PI * pow(r, 2);

  vec4 normal = normalize(triangles[i.triangleIndex].normal);
  vec4 direction = normalize(light.position - i.position );

  float r_n = dot(direction, normal);

  vec3 d = (light.color * max((r_n), 0.f))/area;
  Intersection intersect;

  ClosestIntersection(i.position+triangles[i.triangleIndex].normal*options.bias, direction, triangles, intersect);

  if (glm::distance(i.position, intersect.position) < r) return vec3(0,0,0);

  return d;
}

void Update(Camera &camera, Light &light, bool &escape, bool &is_lookAt)
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
