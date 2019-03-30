#ifndef STRUCTURE_H
#define STRUCTURE_H

#include "functions.h"


struct Photon {
    vec4 hitPos;
    vec3 color;
    char phi, theta; 
    vec4 incomingDirection;
    int flag;
};

struct Light{
    vec4 position;
    vec3 color;
};

struct Intersection
{
  vec4 position;
  float distance;
  bool isTriangle;
  int circleIndex;
  vec4 circleNormal;
  int triangleIndex;
  bool isDummy;
};

struct Camera{
  vec4 position;
  mat4 basis;
  vec3 center;
};

struct Node{
    Photon photon;
    struct Node *left, *right;
};
