#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <limits.h>
#include "random"
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>
#include <math.h>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

float maxFloat = std::numeric_limits<float>::max();
float minFloat = std::numeric_limits<float>::min();

unsigned long maxLong = std::numeric_limits<unsigned long>::max();



struct Photon {
    vec4 hitPos;
    vec3 color;
    char phi, theta; 
    vec4 incomingDirection;
    int flag;
    unsigned long originalIndex;
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

struct WhichPhoton {
    unsigned long index;
    float distance;
};

const int SPECULAR = 0;
const int ABSORPTION = 1;
const int DIFFUSE = 2;
const int TRANSMISSION = 3;

#define PI 3.14159
#define SCREEN_WIDTH 500
#define SCREEN_HEIGHT 500
#define FULLSCREEN_MODE false


mat4 R;
bool escape = false;


float focal_length = SCREEN_HEIGHT / 2;


Light light = {
  .position = vec4(0, -0.7, -0.7, 1.0),
  .color =   14.f * vec3(1,1,1),
};

Camera camera = {
  .position = vec4(0,0,0, 1.0),
  .basis = mat4(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), vec4(0,0, -2, 1.0)),
  // .center = vec3(0.003724, 0.929729, 0.07459)
  .center = vec3(0,0,0)
};

float yaw = 2 * PI / 180;

// DEFINE DISTRIBUTIONS
std::default_random_engine generator;
std::uniform_real_distribution<float> distributionx(light.position[0]-1, light.position[0]+1);
std::uniform_real_distribution<float> distributiony(light.position[1]-1, light.position[1]+1);
std::uniform_real_distribution<float> distributionz(light.position[2]-1, light.position[2]+1);
std::uniform_real_distribution<float> distribution(0 , 1);

void firePhoton(vec4 s, vec4 d, vector<Photon>& photons, Photon& photon, const vector<Triangle>& triangles, const vector<Circle>& circles, int n_fire);
bool ClosestIntersection(vec4 s, vec4 d, const vector<Triangle>& triangles, const vector<Circle>& circles, Intersection& closestIntersection, int index=-1);
bool circleIntersection(vec4 start, vec4 dir, const vector<Circle>& circles, Intersection& closestIntersection, int within_index = -1, bool isInCircle = false);
void startPhotons(vector<Photon>& photons, const vector<Triangle>& triangles, const vector<Circle>& circles);
vec4 calculateLightDirection();
int bounceType(const Photon &photon, const Material& material);
vec3 sampleDirectionVector(const float &r1, const float &r2);
void createCoordinateSystem(const vec3 &N, vec3 &Nt, vec3 &Nb);
int getMaximumDimension( vector<Photon>& photons);
void swap(Photon *a, Photon *b);
int partition(vector <Photon> &photons, int dim, int low, int high);
void quickSort(vector<Photon> &photons, int dim, int low, int high);
vec4 toVec4(const vec3 x);
vec3 toVec3(const vec4 x);
vec4 diffuseDirection(vec4& normal);
vec4 specularReflection(vec4& normal, vec4& incoming);
int getMaximumDimension( vector<Photon>& photons);
void Update();
void Draw(screen* screen, const vector<Triangle>& triangles, vector<Photon>& kdTree, const vector<unsigned long> &photonIndex, const vector<int> &dimensions, const int photonSize, const vector<Circle>& circles);
void scalePower(const int& reflectionType, const Photon& incomingPhoton, Photon& outwardsPhoton, const Material& material);
vec3 processingPart(int row, int col, const vector<Triangle>& triangles, const vector<Circle>& circles, const vector<Photon>& kdTree, const vector<unsigned long> &photonIndex, const int photonSize, const vector<int> &dimensions);
vec3 castRay(vec4 s, vec4 d, const vector<Triangle>& triangles, const vector<Circle>& circles, int bounces, const vector<Photon> &kdTree, const vector<unsigned long> &photonIndex, const vector<int> &dimensions, const int photonSize);
vec3 DirectLight(const Intersection& i, const vector<Triangle>& triangles, const vector<Circle>& circles);
mat4 generateRotation(vec3 a);
void createkdTree(const unsigned long index, vector<unsigned long>& photonIndex, vector<int>& dimensions, vector<Photon>& photons, const int photonSize, bool &isSuccessful);
void searchkdTree(int index, vector<WhichPhoton> &closest, const int maxSize, int &occupied, const vector<unsigned long> &photonIndex, const vector<int> &dimensions, const vector<Photon> &photons, const int photonSize, const vec4 position, float maxDistance);
void max_heapify(vector<WhichPhoton> &closest, int i, int n);
void build_maxheap(vector<WhichPhoton> &closest, int n);

int main(){
    screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );


    vector<Triangle> triangles;
    vector<Photon> photons;
    vector<Circle> circles;
    LoadTestModel(triangles);
    LoadCircles(circles);
    startPhotons(photons, triangles, circles);
    unsigned long index = 0;
    vector<unsigned long> photonIndex;
    vector<int> dimensions;

    int size = 1 + pow(2, int(ceil(log2(photons.size()))));

    printf("Size is %d\n", size);

    for (int i = 0; i < size; i++) {
        photonIndex.push_back(maxLong);
        dimensions.push_back(-2);
    }
    bool isSuccessful = true;

    printf("STORED PHOTONS: %d\n", photons.size());

    createkdTree(index, photonIndex, dimensions, photons, size, isSuccessful);
    
    printf("KDTREE: %d\n", photonIndex.size());

    // printf("Value of photonIndex: %d, %d, %d, %d, %d\n", photonIndex[0], photonIndex[1], photonIndex[2], photonIndex[3], photonIndex[4]);
    // printf("Dimension: %d, %d, %d, %d, %d\n", dimensions[0],dimensions[1], dimensions[2], dimensions[3], dimensions[4]);

      while( !escape )
    {
      Draw(screen, triangles, photons, photonIndex, dimensions, photons.size(), circles);
      SDL_Renderframe(screen);
      Update();

    }
      SDL_SaveImage( screen, "screenshot.bmp" );
}

void Draw(screen* screen, const vector<Triangle>& triangles, vector<Photon>& kdTree, const vector<unsigned long> &photonIndex, const vector<int> &dimensions, const int photonSize, const vector<Circle>& circles)
{

    memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
    for (int row= 0; row < SCREEN_HEIGHT; row++){
        for (int col = 0; col < SCREEN_WIDTH; col++){
            
            vec3 c = processingPart(row, col, triangles, circles, kdTree, photonIndex, photonSize, dimensions);
            PutPixelSDL(screen, row, col, c);

        }

        printf("ROW: %d\n", row);
    }

}

vec3 processingPart(int row, int col, const vector<Triangle>& triangles, const vector<Circle>& circles, const vector<Photon>& kdTree, const vector<unsigned long> &photonIndex, const int photonSize, const vector<int> &dimensions){
    Intersection intersect;
    vec3 color(0.f, 0.f, 0.f);
    vec3 shadow(0.f, 0.f, 0.f);
    vec4 d = camera.basis * vec4(row - SCREEN_WIDTH/2, col - SCREEN_HEIGHT/2, focal_length, 1);
    Intersection inter;

    vec3 mainColor = castRay(camera.basis[3], d, triangles, circles, 0, kdTree, photonIndex, dimensions, photonSize);

    return mainColor;
}

void max_heapify(vector<WhichPhoton> &closest, int i, int n) {
    int j;
    WhichPhoton temp = closest[i];
    j = 2 * i;
    while (j <= n)
    {
        if (j < n && closest[j+1].distance > closest[j].distance)
            j = j + 1;
        if (temp.distance > closest[j].distance)
            break;
        else if (temp.distance <= closest[j].distance)
        {
            closest[j / 2] = closest[j];
            j = 2 * j;
        }
    }
    closest[j/2] = temp;
    return;
}

void build_maxheap(vector<WhichPhoton> &closest, int n) {
    int i;
    for(i = n/2; i >= 1; i--)
    {
        max_heapify(closest,i,n);
    }
}

void searchkdTree(int index, vector<WhichPhoton> &closest, const int maxSize, int &occupied, const vector<unsigned long> &photonIndex, const vector<int> &dimensions, const vector<Photon> &photons, const int photonSize, const vec4 position, float maxDistance) {
    int right = index * 2 + 2;
    int left = index * 2 + 1;

    if (right < photonSize) {
        int dimension = dimensions[index];
        if (dimension > -1) {
            
            float distance = photons[photonIndex[index]].hitPos[dimension] - position[dimension];

            if (distance < 0) {
                searchkdTree(left, closest, maxSize, occupied, photonIndex, dimensions, photons, photonSize, position, maxDistance);
                if (pow(distance,2) < maxDistance) {
                    searchkdTree(right, closest, maxSize, occupied, photonIndex, dimensions, photons, photonSize, position, maxDistance);
                }
            } else {
                searchkdTree(right, closest, maxSize, occupied, photonIndex, dimensions, photons, photonSize, position, maxDistance);

                if (pow(distance, 2) < maxDistance) {
                    searchkdTree(left, closest, maxSize, occupied, photonIndex, dimensions, photons, photonSize, position, maxDistance);
                }
            }
        }
    }

    if (photonIndex[index] != maxLong) {
        float distance = pow(glm::distance(photons[photonIndex[index]].hitPos, position), 2);
        if (distance < maxDistance) {
            if (occupied == maxSize) {
                if (distance < closest[0].distance) {

                    WhichPhoton newPhoton = {
                        .index = photonIndex[index],
                        .distance = distance
                    };

                    closest[0] = newPhoton;
                    
                    build_maxheap(closest, occupied);

                    maxDistance = closest[0].distance;   
                }
            } else {
                occupied++;
                WhichPhoton newPhoton = {
                    .index = photonIndex[index],
                    .distance = distance
                };
                closest[occupied - 1] = newPhoton;
            }   
        }
    }
}



vec3 castRay(vec4 s, vec4 d, const vector<Triangle>& triangles, const vector<Circle>& circles, int bounces, const vector<Photon> &kdTree, const vector<unsigned long> &photonIndex, const vector<int> &dimensions, const int photonSize){
    if (bounces > 2) return vec3(0,0,0);
    Intersection i;
    vec3 hitColor(0,0,0);

    if (ClosestIntersection(s, d, triangles, circles, i, 1)){

        vec4 normal;
        vec3 objectColor;
        bool isMirror;
        if (i.isTriangle){
            objectColor = triangles[i.triangleIndex].color;  
            isMirror = triangles[i.triangleIndex].isMirror;
            normal = triangles[i.triangleIndex].normal;
        }
        else {
            objectColor = circles[i.circleIndex].color;
        }

        if (isMirror){
            vec4 reflectedDirection = specularReflection(normal, d);
            hitColor += 0.8f * castRay(i.position, reflectedDirection, triangles, circles, bounces+1, kdTree, photonIndex, dimensions, photonSize);
        } else {
            vec3 indirect;
            vec3 directLight = DirectLight(i, triangles, circles);
            vec3 acc (0,0,0);
            int index;
            vector<Photon> copy = kdTree;
            int neighbours = 1000;
            int occupied = 0;
            vector<WhichPhoton> closest;
            closest.resize(neighbours);

            searchkdTree(0, closest, neighbours, occupied, photonIndex, dimensions, kdTree, photonSize, i.position, 0.5);   

            for (int j = 0; j < occupied; j++) {
                float weight = (max(glm::dot((kdTree[closest[j].index].incomingDirection), normal), 0.f));
                acc += kdTree[closest[j].index].color * weight;
            }

            float radius = abs(glm::distance(i.position, normalize(kdTree[index].hitPos)));
            float area = float(PI) * pow(radius, 2);
            acc /= area;
            acc /= occupied;
            if (acc.x == 0 && acc.y == 0 && acc.z ==0)std::cout<<glm::to_string(acc)<<std::endl;
            
            hitColor = acc;
        }
    } else {
        return vec3(0,0,0);
    }
    return hitColor;
}

void swap(Photon *a, Photon *b){
    Photon *temp = a;
    *a = *b;
    *b = *temp;
}

int getMaximumDimension( vector<Photon>& photons){
    float difference = minFloat;
    int dim = -1;

    for (int dimension = 0; dimension < 3; dimension++){
        float max = minFloat; 
        float min = maxFloat; 

        for (uint i = 0; i < photons.size(); i++){
            if (photons[i].hitPos[dimension] > max) max = photons[i].hitPos[dimension];
            if (photons[i].hitPos[dimension] < min) min = photons[i].hitPos[dimension];
        }
        if (max - min > difference) {
            difference = max-min;
            dim = dimension;
        }
    }

    return dim;
}

int partition(vector <Photon> &photons, int dim, int low, int high){
    int pivot = photons[high].hitPos[dim];
    int i = (low - 1);

    for (int j = low; j <= high - 1; j++)
    {
        if (photons[j].hitPos[dim] <= pivot)
        {
            i++;
            swap(&photons[i], &photons[j]);
        }
    }
    swap(&photons[i+1], &photons[high]);
    return (i + 1);
}

void quickSort(vector<Photon> &photons, int dim, int low, int high){
    if (low < high) 
    {
        int pi = partition(photons, dim, low, high);
        quickSort(photons, dim, low, pi-1);
        quickSort(photons, dim, pi+1, high);
    }
}

int calculateMedian(vector<Photon> photons, int& dim){
    quickSort(photons, dim,  0, photons.size()-1);

    return floor(photons.size()/2);
}

void createkdTree(const unsigned long index, vector<unsigned long>& photonIndex, vector<int>& dimensions, vector<Photon>& photons, const int photonSize, bool &isSuccessful){
    if (photons.size() == 0) {
        isSuccessful = false;
        return;
    }
    if (photons.size() == 1) {
        isSuccessful = true;
        if (index < photonSize) {
            photonIndex[index] = photons[0].originalIndex;
            dimensions[index] = -1;
        }
        return;
    }
    int dim = getMaximumDimension(photons);

    int median = calculateMedian(photons, dim);

    photonIndex[index] = photons[median].originalIndex;
    dimensions[index] = dim;

    vector <Photon> left;
    vector <Photon> right;

    for (uint i = 0; i < photons.size(); i++){
        if (i < median){
            left.push_back(photons[i]);
        }
        if (i > median){
            right.push_back(photons[i]);
        }
    }


    createkdTree(2 * index + 1, photonIndex, dimensions, left, photonSize, isSuccessful);
    createkdTree(2 * index + 2, photonIndex, dimensions, right, photonSize, isSuccessful);
    return;
}

void firePhoton(vec4 s, vec4 d, vector<Photon>& photons, Photon& photon, const vector<Triangle>& triangles, const vector<Circle>& circles, int n_fire){
    Intersection i;
    if (photon.color.x < 0.1 && photon.color.y < 0.1 && photon.color.z < 0.1) {
        return;
    }
    if (ClosestIntersection(s, d, triangles, circles, i, n_fire)){
        vec4 objectNormal;
        vec3 objectColor;
        Material objectMaterial;
        bool isCircle = false;
        if (i.isTriangle) {
            objectNormal = triangles[i.triangleIndex].normal;
            objectColor = triangles[i.triangleIndex].color;
            objectMaterial = triangles[i.triangleIndex].material;
        }
        else {
            isCircle = true;
            objectNormal = i.circleNormal;
            objectColor = circles[i.circleIndex].color;
            objectMaterial = circles[i.circleIndex].material;
        }


        // HERE IS WHERE THE DIFFERENT BOUNCES WILL HAPPEN ONLY STORE WHEN PHOTON HITS A DIFFUSE
        switch (bounceType(photon, objectMaterial))
        {
            case DIFFUSE: 
            {
                // photon.color *= objectColor;

                int theta = int(acos(d[2])*(256.0/PI) );
                if (theta>255)
                    photon.theta = 255;
                else
                    photon.theta = (unsigned char) theta;

                int phi = int(atan2(d[1],d[0])*(256.0/(2.0*PI)) );
                if (phi>255)
                    photon.phi = 255;
                else if (phi<0)
                    photon.phi = (unsigned char)(phi+256);
                else
                    photon.phi = (unsigned char)phi;
                photon.incomingDirection = normalize(d);
                photon.hitPos = i.position;
                photon.originalIndex = photons.size();
                photons.push_back(photon);
                

                    
                Photon newPhoton;
                scalePower(DIFFUSE, photon, newPhoton, objectMaterial);
                newPhoton.color *= objectColor;
                vec4 outDir = diffuseDirection(objectNormal);
                firePhoton(i.position + toVec4(outDir) * 0.001f, outDir, photons, newPhoton, triangles, circles, n_fire + 1);

                /* code */
                break;
            }
            case SPECULAR:
            {
                photon.hitPos = i.position;

                Photon newPhoton;
                scalePower(SPECULAR, photon, newPhoton, objectMaterial);
                vec4 outDir = specularReflection(objectNormal, d);
                firePhoton(i.position, outDir, photons, newPhoton, triangles, circles, n_fire + 1);
                break;
            }
            case ABSORPTION:
                break;
            default:
                break;
        }

    } else {
        return;
        printf("NO HIT\n");
    }

}

vec3 DirectLight(const Intersection& i, const vector<Triangle>& triangles, const vector<Circle>& circles){
    float r = glm::distance(light.position, i.position);
    float area = 4 * PI * pow(r, 2);
    vec4 normal; 
    if (i.isTriangle) {
        normal = triangles[i.triangleIndex].normal;
    }
    else {
        normal = i.circleNormal;
    }

    vec4 direction = normalize(light.position - i.position);
    float r_n = dot(direction, normal);
    vec3 d = (light.color * max((r_n), 0.f))/area;
    
    Intersection intersect;
    if (ClosestIntersection(i.position, direction, triangles, circles, intersect, i.triangleIndex)){
        if (glm::distance(i.position, intersect.position) < r && glm::distance(i.position, intersect.position) >= 1e-4) {
        return vec3(0,0,0);
        }
    }

    return d;

    
}

vec4 calculateLightDirection(){
    float x = distributionx(generator);
    float y = distributiony(generator);
    float z = distributionz(generator);

    vec4 randomPoint (x,y,z, 1);

    vec4 direction = randomPoint - light.position;

    return normalize(direction);

}

void startPhotons(vector<Photon>& photons, const vector<Triangle>& triangles, const vector<Circle>& circles) {
    int numberOfPhotons = 300000;

    for (int i = 0; i < numberOfPhotons; i++){
        Photon photon;
        photon.color = light.color;

        vec4 direction = calculateLightDirection();
        firePhoton(light.position, direction, photons, photon, triangles, circles, 1);

    }

}

void scalePower(const int& reflectionType, const Photon& incomingPhoton, Photon& outwardsPhoton, const Material& material){
    switch (reflectionType)
    {
        case DIFFUSE:
            /* code */
            outwardsPhoton.color = incomingPhoton.color * material.diffuse;
            break;
        case SPECULAR:
            outwardsPhoton.color = incomingPhoton.color * material.specular;
        default:
            break;
    }
}




int bounceType(const Photon &photon, const Material& material) {
    float randomNumber = distribution(generator);
    if (randomNumber < material.specular) {
        return SPECULAR;
    } else if (randomNumber < material.specular + material.refraction){
        return TRANSMISSION;
    } else if (randomNumber < material.specular + material.refraction + material.absorption){
        return ABSORPTION;
    } else if (randomNumber <= material.specular + material.refraction + material.absorption + material.diffuse){
        return DIFFUSE;
    } else {
        printf("THERE WAS AN ERROR AT BOUNCE TYPE\n");
        return -1;
    }
    return 0;
}



bool ClosestIntersection(vec4 s, vec4 d, const vector<Triangle>& triangles, const vector<Circle>& circles, Intersection& closestIntersection, int index){

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
  
    circleIntersection(s, d, circles, closestIntersection, index);

  if (closestIntersection.distance == maxFloat) return false;
  return true;
}

bool circleIntersection(vec4 s, vec4 d, const vector<Circle>& circles, Intersection& closestIntersection, int within_index, bool isInCircle) {

  if (isInCircle) closestIntersection.distance = maxFloat;

  float current = closestIntersection.distance;
  for (uint i = 0; i < circles.size(); i++) {

    if (isInCircle && (i != within_index)) continue;
    // if (!isInCircle && within_index != -1 && (i == within_index)) continue;
   
    float t0 = maxFloat;
    float t1 = maxFloat;

    vec4 L4 = s - circles[i].center;

    vec3 L = vec3(L4);
    vec3 d3 = vec3(d);

    float a = dot(d3, d3);
    float b = 2 * dot(d3, L);
    float c = dot(L, L) - pow(circles[i].radius, 2);

    if ((pow(b,2) - 4 * a * c) < 0) continue;

    // if (within_index == 2 && (pow(b,2) - 4 * a * c) > 0) printf("I'm here");
    
    

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
      closestIntersection.circleNormal = normalize((s + t0 * d) - circles[i].center);
    }
  }

  if (closestIntersection.distance == current) return false;
  return true;
}


vec4 diffuseDirection(vec4& normal){
    vec3 hitNormal = toVec3(normal);
    vec3 Nt, Nb;
    createCoordinateSystem(hitNormal, Nt, Nb);
    float r1 = distribution(generator);
    float r2 = distribution(generator);
    
    vec3 sampledVector = sampleDirectionVector(r1, r2);
    vec3 sampleWorld(
      sampledVector.x * Nb.x + sampledVector.y * hitNormal.x + sampledVector.z * Nt.x,
      sampledVector.x * Nb.y + sampledVector.y * hitNormal.y + sampledVector.z * Nt.y,
      sampledVector.x * Nb.z + sampledVector.y * hitNormal.z + sampledVector.z * Nt.z);
    
    sampleWorld = (sampleWorld);
    return normalize(toVec4(sampleWorld));

}

vec3 toVec3(const vec4 x){
  return vec3(x[0], x[1], x[2]);
}

vec4 toVec4(const vec3 x) {
  return vec4(x.x, x.y, x.z, 1);
}

void createCoordinateSystem(const vec3 &N, vec3 &Nt, vec3 &Nb){
  if (std::fabs(N[0]) > std::fabs(N[1])){
    Nt = vec3(N[2], 0, -N[0]) / sqrt(N[0] * N[0] + N[2] * N[2]);
  } else {
    Nt = vec3(0, -N[2], N[1]) / sqrt(N[1] * N[1] + N[2] * N[2]);
  }
  Nb = (glm::cross(N, Nt));
}

vec3 sampleDirectionVector(const float &r1, const float &r2){
  float sinTheta = sqrtf(1 - r1 * r1);
  float phi = 2 * PI * r2;
  float x = sinTheta * cosf(phi);
  float z = sinTheta * sinf(phi);
  return vec3(x, r1, z);
}

vec4 specularReflection(vec4& normal, vec4& incoming) {
    vec4 n = normal;
    vec4 reflectedDirection = (incoming - 2.f * (glm::dot(incoming, n)) * n);
    return normalize(reflectedDirection);
}

mat4 generateRotation(vec3 a){
  vec4 uno = vec4(cos(a[1])*cos(a[2]), cos(a[1])*sin(a[2]), -sin(a[1]), 0);
  vec4 dos = vec4(-cos(a[0])*sin(a[2]) + sin(a[0]*sin(a[1])*cos(a[2])), cos(a[0])*cos(a[2])+sin(a[0])*sin(a[1])*sin(a[2]), sin(a[0])*cos(a[1]), 0);
  vec4 tres = vec4(sin(a[0])*sin(a[2])+cos(a[0])*sin(a[1])*cos(a[2]), -sin(a[0])*cos(a[2])+cos(a[0])*sin(a[1])*sin(a[2]), cos(a[0])*cos(a[1]), 0);
  vec4 cuatro = vec4(0,0,0,1);
  return (mat4(uno, dos, tres, cuatro));
}



void Update()
{
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;
  /*Good idea to remove this*/
//   std::cout << "Render time: " << dt << " ms." << std::endl;
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
                camera.basis[3][2] += 0.1;
                break;
            case SDLK_DOWN:
            // Move camera backward
                camera.basis[3][2] += -0.1;
                break;
            case SDLK_LEFT:
            // Move camera to the left
                camera.basis[3][0] += -0.1;
                break;
            case SDLK_RIGHT:
            // Move camera to the right
                camera.basis[3][0] += 0.1;
                break;
            case SDLK_n:
                camera.basis[3][1] += -0.1;
                break;
            case SDLK_m:
                camera.basis[3][1] += 0.1;
                break;
            case SDLK_d:
                //Rotate camera right;
                camera.basis =  generateRotation(vec3(0, yaw, 0)) * camera.basis;
                break;
            case SDLK_a:
                //Rotate camera left;
                camera.basis =  generateRotation(vec3(0, -yaw, 0)) * camera.basis;

                break;
            case SDLK_w:
                //Rotate camera top;
                camera.basis =  generateRotation(vec3(yaw, 0, 0)) * camera.basis;

                break;
            case SDLK_s:
                //Rotate camera down;
                camera.basis =  generateRotation(vec3(-yaw, 0, 0)) * camera.basis;
                break;
            case SDLK_q:
                // camera.basis =  generateRotation(vec3(0, 0, -yaw)) * camera.basis;
                break;
            case SDLK_e:
                // camera.basis =  generateRotation(vec3(0, 0, yaw)) * camera.basis;
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
   }
 }
}