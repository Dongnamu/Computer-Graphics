#ifndef FUNCTION_H
#define FUNCTION_H

#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include <glm/gtx/string_cast.hpp>
#include <iostream>


void firePhoton(vec4 s, vec4 d, vector<Photon>& photons, Photon& photon, const vector<Triangle>& triangles, const vector<Circle>& circles, int n_fire);
bool ClosestIntersection(vec4 s, vec4 d, const vector<Triangle>& triangles, const vector<Circle>& circles, Intersection& closestIntersection, int index=-1);
bool circleIntersection(vec4 start, vec4 dir, const vector<Circle>& circles, Intersection& closestIntersection, int within_index = -1, bool isInCircle = false);
void startPhotons(vector<Photon>& photons, const vector<Triangle>& triangles, const vector<Circle>& circles);
vec4 calculateLightDirection();
int bounceType(const Photon &photon, const Material& material);
vec3 sampleDirectionVector(const float &r1, const float &r2);
void createCoordinateSystem(const vec3 &N, vec3 &Nt, vec3 &Nb);
vec4 toVec4(const vec3 x);
vec3 toVec3(const vec4 x);
vec4 diffuseDirection(vec4& normal);
vec4 specularReflection(vec4& normal, vec4& incoming);
int getMaximumDimension( vector<Photon>& photons);
void Update();
void Draw(screen* screen, const vector<Triangle>& triangles, vector<Photon>& kdTree, const vector<Circle>& circles);
void scalePower(const int& reflectionType, const Photon& incomingPhoton, Photon& outwardsPhoton, const Material& material);
vec3 processingPart(int row, int col, const vector<Triangle>& triangles, const vector<Circle>& circles, const vector<Photon>& kdTree);
vec3 castRay(vec4 s, vec4 d, const vector<Triangle>& triangles, const vector<Circle>& circles, int bounces, const vector<Photon> kdTree);
vec3 DirectLight(const Intersection& i, const vector<Triangle>& triangles, const vector<Circle>& circles);
void linearNN(vector<Photon>& kdTree, vec3& col, const vec4 targetPoint, int &number);