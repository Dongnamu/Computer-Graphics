#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <stdint.h>

using namespace std;
using glm::vec3;
using glm::mat3;

#define SCREEN_WIDTH 640
#define SCREEN_HEIGHT 518
#define FULLSCREEN_MODE false
#define FOCAL_LENGTH 259



/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */
//int t;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update(vector<vector<vec3>>& trails);
void Draw(screen* screen, float* offset, vector<vector<vec3>>& trails);
void Interpolate( float a, float b, vector<float>& result );
void InterpolateVector( vec3 a, vec3 b, vector<vec3>& result );

int main( int argc, char* argv[] )
{
  // vector<vec3> result( 4 );
  // vec3 a(1,4,9.2);
  // vec3 b(4,1,9.8);
  // InterpolateVector( a, b, result );
  // for( uint i=0; i<result.size(); ++i )
  // {
  //   cout << "( "
  //   << result[i].x << ", "
  //   << result[i].y << ", "
  //   << result[i].z << " ) ";
  // }

  // Starfield
  vector<vec3> stars(1000);
  vector<vector<vec3>> trails (1000);
  for (uint i = 0; i < trails.size(); i++){
    trails[i] = vector<vec3>(10);
  }
  for (uint trail = 0; trail < trails.size(); trail++){
    float x = 1.0 - (float(rand())/float(RAND_MAX) * 2.0);
    float y = 1.0 - (float(rand())/float(RAND_MAX) * 2.0);
    while (true){
      x = 1.0 - (float(rand())/float(RAND_MAX) * 2.0);
      y = 1.0 - (float(rand())/float(RAND_MAX) * 2.0);
      if (x != 0.0 && y != 0.0) break;
    }
    float z = float(rand())/float(RAND_MAX);
    for (uint star = 0; star < trails[trail].size(); star++){
      trails[trail][star][0] = x;
      trails[trail][star][1] = y;
      if (star == 0){
        trails[trail][star][2] = z;
      } else {
        trails[trail][star][2] = trails[trail][star-1][2] + 0.0025;
      }

    }
  }
  // for (uint i = 0; i < stars.size(); i++){
  //
  //   stars[i][0] = x;
  //   stars[i][1] = y;
  //   stars[i][2] = float(rand())/float(RAND_MAX);
  // }


  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  //t = SDL_GetTicks();	/*Set start value for timer.*/
  float i = 0;

  while( NoQuitMessageSDL() )
    {
      Draw(screen, &i, trails);
      Update(trails);
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}
void Interpolate(float a, float b, vector<float>& result){
  if (result.size() == 1){
    result[0] = a;
  }
  float distance = abs((a-b)/(result.size()-1));
  float current = a;
  for (uint i =0; i < result.size(); i++){
    result[i]= current;
    current += distance;
  }
}

void InterpolateVector(vec3 a, vec3 b, vector<vec3>& result){
  float stepX = ((b[0] - a[0])/(result.size()-1));
  float stepY = ((b[1] - a[1])/(result.size()-1));
  float stepZ = ((b[2] - a[2])/(result.size()-1));

  vec3 current = a;
  for (uint i =0; i < result.size(); i++){
    result[i]= current;
    current[0] += stepX;
    current[1] += stepY;
    current[2] += stepZ;
  }

}

/*Place your drawing here*/
void Draw(screen* screen, float* offset, vector<vector<vec3>>& trails)
{
  // vec3 topLeft(1,0,0);
  // vec3 topRight(0,0,1+*offset);
  // vec3 bottomRight(0,1+*offset,0);
  // vec3 bottomLeft(1,1-*offset,0);
  //
  // vector<vec3> leftSide(SCREEN_HEIGHT);
  // vector<vec3> rightSide(SCREEN_HEIGHT);
  // InterpolateVector(topLeft, bottomLeft, leftSide);
  // InterpolateVector(topRight, bottomRight, rightSide);
  //


  vec3 colour(255,255,255);


  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
  for (uint trail=0; trail< trails.size(); trail++){
    for (uint i=0; i < trails[trail].size(); i++){
      vec3 color = 0.2f * vec3(1,1,1) / (trails[trail][i].z*trails[trail][i].z);
      float u = (FOCAL_LENGTH * (trails[trail][i][0]/trails[trail][i][2])) + (screen-> width)/2;
      float v = (FOCAL_LENGTH * (trails[trail][i][1]/trails[trail][i][2])) + (screen-> height)/2;
      PutPixelSDL(screen, u, v, color);
    }
  }


  // for(int i=0; i<1000; i++)
  //   {
  //     uint32_t x = rand() % screen->width;
  //     uint32_t y = rand() % screen->height;
  //     PutPixelSDL(screen, x, y, colour);
  //   }
  // for (int i = 0; i < screen->height; i++){
  //   vector<vec3> row (SCREEN_WIDTH);
  //   InterpolateVector(leftSide[i], rightSide[i], row);
  //   for (int j = 0; j < screen -> width; j++){
  //     PutPixelSDL(screen, j, i, row[j]);
  //   }
  // }
  //*offset = (*offset - 0.001);

    // for (int i = 0; i< screen->width; i++){
    //   for(int j = 0; j < screen->height; j++){
    //
    //     PutPixelSDL(screen, i, j, colour);
    //   }
    // }
}


/*Place updates of parameters here*/
void Update(vector<vector<vec3>>& trails)
{
  /* Compute frame time */
  static int t = SDL_GetTicks();
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;
  for (uint trail = 0; trail < trails.size(); trail++){
    for (uint i = 0; i < trails[trail].size(); i++){
      trails[trail][i].z = (trails[trail][i].z - (0.0005*dt));
      if(trails[trail][i].z <= 0){
        // trails[trail][i].z += 1;
      }
    }
  }

  /*Good idea to remove this*/
  std::cout << "Render time: " << dt << " ms." << std::endl;
  /* Update variables*/
}
