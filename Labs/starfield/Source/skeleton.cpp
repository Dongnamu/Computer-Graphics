#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <stdint.h>

using namespace std;
using glm::vec3;
using glm::mat3;

#define SCREEN_WIDTH 380
#define SCREEN_HEIGHT 240
#define FULLSCREEN_MODE false

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */
int t;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update(vector<vec3> & stars);
void Draw(screen* screen, vector<vec3> & stars);
void Interpolate( vec3 a, vec3 b, vector<vec3> & result );


int main( int argc, char* argv[] )
{

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  t = SDL_GetTicks();	/*Set start value for timer.*/

  vector<vec3> stars(1000);

  for (size_t i = 0; i < stars.size(); i++) {
    stars[i].x = (float(rand()) - float(RAND_MAX / 2))/float(RAND_MAX / 2);
    stars[i].y = (float(rand()) - float(RAND_MAX / 2))/float(RAND_MAX / 2);
    stars[i].z = float(rand())/float(RAND_MAX);
  }

  while( NoQuitMessageSDL() )
    {
      Update(stars);
      Draw(screen, stars);
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  // vector<vec3> result(4);
  // vec3 a(1, 4, 9.2);
  // vec3 b(4, 1, 9.8);
  // Interpolate( a, b, result );
  // for (uint i = 0; i < result.size(); ++i) {
  //   cout << "("
  //        << result[i].x << ","
  //        << result[i].y << ","
  //        << result[i].z << ")";
  // }

  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen* screen, vector<vec3> & stars)
{
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  float f = SCREEN_HEIGHT / 2;


  for (size_t s = 0; s < stars.size(); s++) {
    float u = f * stars[s].x / stars[s].z + SCREEN_WIDTH / 2;
    float v = f * stars[s].y / stars[s].z + f;
    vec3 colour = 0.2f * vec3(1,1,1) / (stars[s].z*stars[s].z);

    PutPixelSDL(screen, u, v, colour);
  }


  // Drawing for interpolate colors
  //
  // vec3 topLeft(1,0,0);
  // vec3 topRight(0,0,1);
  // vec3 bottomRight(0,1,0);
  // vec3 bottomLeft(1,1,0);
  //
  // vector<vec3> leftSide (SCREEN_HEIGHT);
  // vector<vec3> rightSide (SCREEN_HEIGHT);
  // Interpolate( topLeft, bottomLeft, leftSide);
  // Interpolate( topRight, bottomRight, rightSide);

  // for (int i = 0; i < SCREEN_HEIGHT; i++) {
  //   vec3 left = leftSide[i];
  //   vec3 right = rightSide[i];
  //
  //   vector<vec3> inter (SCREEN_WIDTH);
  //
  //   Interpolate(left, right, inter);
  //
  //   for (int j = 0; j < SCREEN_WIDTH; j++) {
  //     PutPixelSDL(screen, j, i, inter[j]);
  //   }
  // }

  // Drawing for plain colors
  //
  // vec3 colour(255.0,255.0,255.0);
  // for (int i = 0; i < SCREEN_WIDTH; i++) {
  //   for (int j = 0; j < SCREEN_HEIGHT; j++) {
  //     PutPixelSDL(screen, i, j, colour);
  //   }

  // for(int i=0; i<1000; i++)
  //   {
  //     uint32_t x = rand() % screen->width;
  //     uint32_t y = rand() % screen->height;
  //     PutPixelSDL(screen, x, y, colour);
  //   }
  //
}

void Interpolate( vec3 a, vec3 b, vector<vec3> & result ) {

  vec3 steps;

  if (result.size() != 0){
    for (int i = 0; i < 3; i++) {
      steps[i] = (b[i]-a[i]) / (result.size() - 1);
    }

    vec3 vectors = a;

    for (uint j = 0; j < result.size();j++) {
      result[j] = vectors;
      for (int i = 0; i < 3; i++) {
        vectors[i] += steps[i];
      }
    }
  } else {
    result[0] = a;
  }


}


/*Place updates of parameters here*/
void Update(vector<vec3> &stars)
{
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;

  for (size_t i = 0; i < stars.size(); i++) {
    stars[i].z -= 0.0002 * dt;
    if (stars[i].z <= 0) {
      stars[i].z += 1;
    }
    // if (stars[i].z > 1) {
      // stars[i].z -= 1;
    // }
  }
  /*Good idea to remove this*/
  std::cout << "Render time: " << dt << " ms." << std::endl;
  /* Update variables*/
}
