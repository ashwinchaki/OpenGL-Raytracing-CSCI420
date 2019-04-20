/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Ashwin Chakicherla
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>
#include <iostream>
#include <glm/glm.hpp>
#include "ray.h"

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

#define FURTHEST -1e8

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 320
#define HEIGHT 240

// aspect ratio of image render
double aspectRatio = (double) WIDTH / (double) HEIGHT;
double angle;

//the field of view of the camera
#define fov 60.0
#define PI 3.14159265

unsigned char buffer[HEIGHT][WIDTH][3];

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];
glm::vec3 ambient;

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

void clamp(double& num, float low, float high){
  if (num > high){
    num = high;
  }
  else if (num < low){
    num = low;
  }
}

<<<<<<< HEAD
// std::vector<Light> getLightSamples(Light light){
//   std::vector<Light> samples;
//   for (int i = 0; i < SHADOW_SAMPLES; i++){
//     // get random phi, theta, and radius
//     double phi = ((double) (rand() % 360) * PI / 180);
//     double theta = ((double) (rand() % 360) * PI / 180);
//     double radius = ((double) rand() / (RAND_MAX)) / 10.0f;

//     glm::vec3 newPos(light.position[0] + radius * std::sin(theta) * std::cos(phi), /* x offset */
//                      light.position[1] + radius * std::sin(theta) * std::sin(phi), /* y offset */
//                      light.position[2] + radius * std::cos(theta));                /* z offset */

//     Light newLight;
//     newLight.position[0] = newPos[0]; newLight.position[1] = newPos[1]; newLight.position[2] = newPos[2];
//     newLight.color[0] = light.color[0] / SHADOW_SAMPLES; newLight.color[1] = light.color[1] / SHADOW_SAMPLES; newLight.color[2] = light.color[2] / SHADOW_SAMPLES;
//     samples.push_back(newLight);
//   }

//   return samples;
// }

=======
>>>>>>> parent of 3b2bb8a... working main, broken soft shadows
glm::vec3 calculateLighting(Sphere& sphere, Light& light, glm::vec3& intersect){
  glm::vec3 normal = glm::normalize(intersect - glm::vec3(sphere.position[0], sphere.position[1], sphere.position[2]));
  glm::vec3 light_direction = glm::normalize(glm::vec3(light.position[0], light.position[1], light.position[2]) - intersect);
  
  // calculate light magnitude --> light_direction dot normal vector
  double light_magnitude = glm::dot(light_direction, normal);
  clamp(light_magnitude, 0, 1);

  // calculate reflection magnitude --> reflection vector dot normalized intersection vector
  glm::vec3 norm_intersect = glm::normalize(intersect * -1.0f);
  glm::vec3 reflect(2 * light_magnitude * normal[0] - light_direction[0],
                    2 * light_magnitude * normal[1] - light_direction[1],
                    2 * light_magnitude * normal[2] - light_direction[2]);
  glm::normalize(reflect);
  double reflection_magnitude = glm::dot(reflect, norm_intersect);
  clamp(reflection_magnitude, 0.0f, 1.0f);

  glm::vec3 diffuse(sphere.color_diffuse[0], sphere.color_diffuse[1], sphere.color_diffuse[2]);
  glm::vec3 spec(sphere.color_specular[0], sphere.color_specular[1], sphere.color_specular[2]);

  // phong shading algorithm

  return glm::vec3(light.color[0] * (diffuse[0] * light_magnitude + (spec[0] * pow(reflection_magnitude, sphere.shininess))),
                   light.color[1] * (diffuse[1] * light_magnitude + (spec[1] * pow(reflection_magnitude, sphere.shininess))),
                   light.color[2] * (diffuse[2] * light_magnitude + (spec[2] * pow(reflection_magnitude, sphere.shininess))));
}

glm::vec3 calculateLighting(Triangle& triangle, Light& light, glm::vec3& intersect){

  glm::vec3 light_direction = glm::normalize(glm::vec3(light.position[0], light.position[1], light.position[2]) - intersect);

  glm::vec3 point_a = glm::vec3(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
  glm::vec3 point_b = glm::vec3(triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
  glm::vec3 point_c = glm::vec3(triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);

  glm::vec3 vector_ab = point_b - point_a;
  glm::vec3 vector_ac = point_c - point_a;
  glm::vec3 vector_cb = point_c - point_b;
  glm::vec3 vector_ca = point_a - point_c;

  glm::vec3 intersect_distance_b = intersect - point_b;
  glm::vec3 intersect_distance_c = intersect - point_c;

  glm::vec3 cross_CB = glm::cross(vector_cb, intersect_distance_b);
  glm::vec3 cross_AC = glm::cross(vector_ca, intersect_distance_c); 

  glm::vec3 planar_vec = glm::cross(vector_ab, vector_ac);
  float denom = glm::dot(planar_vec, planar_vec);

  double u = glm::dot(planar_vec, cross_CB) / denom;
  double v = glm::dot(planar_vec, cross_AC) / denom;
  double w = 1.0 - u - v;

  glm::vec3 normal(u * triangle.v[0].normal[0] + v * triangle.v[1].normal[0] + w * triangle.v[2].normal[0],
                   u * triangle.v[0].normal[1] + v * triangle.v[1].normal[1] + w * triangle.v[2].normal[1],
                   u * triangle.v[0].normal[2] + v * triangle.v[1].normal[2] + w * triangle.v[2].normal[2]);
  glm::normalize(normal);

  double light_magnitude = glm::dot(light_direction, normal);
  clamp(light_magnitude, 0, 1);

  // calculate reflection magnitude --> reflection vector dot normalized intersection vector
  glm::vec3 norm_intersect = glm::normalize(intersect * -1.0f);
  glm::vec3 reflect(2 * light_magnitude * normal[0] - light_direction[0],
                    2 * light_magnitude * normal[1] - light_direction[1],
                    2 * light_magnitude * normal[2] - light_direction[2]);
  glm::normalize(reflect);
  double reflection_magnitude = glm::dot(reflect, norm_intersect);
  clamp(reflection_magnitude, 0, 1);

  glm::vec3 diffuse(u * triangle.v[0].color_diffuse[0] + v * triangle.v[1].color_diffuse[0] + w * triangle.v[2].color_diffuse[0],
                    u * triangle.v[0].color_diffuse[1] + v * triangle.v[1].color_diffuse[1] + w * triangle.v[2].color_diffuse[1],
                    u * triangle.v[0].color_diffuse[2] + v * triangle.v[1].color_diffuse[2] + w * triangle.v[2].color_diffuse[2]);

  glm::vec3 spec(u * triangle.v[0].color_specular[0] + v * triangle.v[1].color_specular[0] + w * triangle.v[2].color_specular[0],
                 u * triangle.v[0].color_specular[1] + v * triangle.v[1].color_specular[1] + w * triangle.v[2].color_specular[1],
                 u * triangle.v[0].color_specular[2] + v * triangle.v[1].color_specular[2] + w * triangle.v[2].color_specular[2]);

  double shininess = u * triangle.v[0].shininess + v * triangle.v[1].shininess + w * triangle.v[2].shininess;

  return glm::vec3(light.color[0] * (diffuse[0] * light_magnitude + (spec[0] * pow(reflection_magnitude, shininess))),
                   light.color[1] * (diffuse[1] * light_magnitude + (spec[1] * pow(reflection_magnitude, shininess))),
                   light.color[2] * (diffuse[2] * light_magnitude + (spec[2] * pow(reflection_magnitude, shininess))));
}

<<<<<<< HEAD
// FIXME: this needs to be fixed, not 100% working now

// glm::vec3 calculateCollisions(Ray &ray, glm::vec3 &color, double &nearest, bool alpha)
// {
//   glm::vec3 returnColor = color;

//   // iterate through each sphere
//   for (int i = 0; i < num_spheres; i++)
//   {
//     glm::vec3 intersect(0, 0, FURTHEST);

//     bool checkIntersect = ray.checkIntersectSpheres(spheres[i], intersect);
//     // FIXME: COUT STATEMENT
//     // std::cout << std::boolalpha << checkIntersect << " && " << intersect[2] << " && " << nearest << std::endl;
//     if (checkIntersect && intersect[2] > nearest)
//     {
//       // FIXME: COUT STATEMENT
//       std::cout << "INTERSECTED SPHERE" << std::endl;
//       returnColor = glm::vec3(0.0f);
//       // now check  for each light
//       for (int j = 0; j < num_lights; j++)
//       {
//         bool light_on = true;

//         std::vector<Light> samples = getLightSamples(lights[j]);
        
//         for (int l = 0; l < samples.size(); l++){
//           glm::vec3 current_light_pos(samples[l].position[0], samples[l].position[1], samples[l].position[2]);
//           Ray shadow_ray(intersect, glm::normalize(current_light_pos - intersect));

//           // check shadow ray intersections with all other objects
//           for (int k = 0; k < num_spheres; k++)
//           {
//             // check for sphere intersections now, excluding current object
//             std::cout << num_spheres << std::endl;
//             glm::vec3 shadow_intersect;
//             if (shadow_ray.checkIntersectSpheres(spheres[k], shadow_intersect) && k != i)
//             {
//               if (glm::length(shadow_intersect - intersect) < glm::length(current_light_pos - intersect))
//               {
//                 light_on = false;
//                 break;
//               }
//             }
//           }

//           // checks for triangle intersections second
//           for (int k = 0; k < num_triangles; k++)
//           {
//             glm::vec3 shadow_intersect;
//             if (shadow_ray.checkIntersectTriangles(triangles[k], shadow_intersect))
//             {
//               if (glm::length(shadow_intersect - intersect) < glm::length(current_light_pos - intersect))
//               {
//                 light_on = false;
//                 break;
//               }
//             }
//           }

//           if (light_on)
//           {
//             returnColor = returnColor + calculateLighting(spheres[i], samples[l], intersect);
//           }
//         }
//       }

//       // new closest intersect
//       nearest = intersect[2];
//     }
//   }

//   // now check triangle collisions
//   for (int i = 0; i < num_triangles; i++)
//   {
//     glm::vec3 intersect(0, 0, FURTHEST);
//     bool checkIntersect = ray.checkIntersectTriangles(triangles[i], intersect);
//     // FIXME: COUT STATEMENT
//     // std::cout << std::boolalpha << checkIntersect << " && " << (intersect[2] > nearest) << std::endl;
//     if (checkIntersect && intersect[2] > nearest)
//     {
//       // FIXME: COUT STATEMENT
//       std::cout << "INTERSECTED TRIANGLE" << std::endl;
//       returnColor = glm::vec3(0.0f);

//       for (int j = 0; j < num_lights; j++)
//       {

//         std::vector<Light> samples = getLightSamples(lights[j]);

//         for (int l = 0; l < samples.size(); l++){
//           glm::vec3 current_light_pos(samples[l].position[0], samples[l].position[1], samples[l].position[2]);
//           Ray shadow_ray(intersect, glm::normalize(current_light_pos - intersect));

//           bool light_on = true;

//           for (int k = 0; k < num_spheres; k++)
//           {
//             glm::vec3 shadow_intersect(0, 0, FURTHEST);
//             if (shadow_ray.checkIntersectSpheres(spheres[k], shadow_intersect)){
//               if (glm::length(shadow_intersect - intersect) < glm::length(current_light_pos - intersect)){
//                 light_on = false;
//                 break;
//               }
//             }
//           }

//           for (int k = 0; k < num_triangles; k++)
//           {
//             glm::vec3 shadow_intersect(0, 0, FURTHEST);
//             if (shadow_ray.checkIntersectTriangles(triangles[k], shadow_intersect) && k != i)
//             {
//               if (glm::length(shadow_intersect - intersect) < glm::length(current_light_pos - intersect))
//               {
//                 light_on = false;
//                 break;
//               }
//             }
//           }
//           if (light_on)
//           {
//             returnColor = returnColor + calculateLighting(triangles[i], lights[j], intersect);
//           }
//         }
//         // new shadow ray starts from intersection point, towards the light source
        
//       }

//       //new closest intersect
//       nearest = intersect[2];
//     }
//   }

//   return returnColor;
// }

=======
>>>>>>> parent of 3b2bb8a... working main, broken soft shadows
glm::vec3 calculateCollisions(Ray& ray, glm::vec3& color, double& nearest){
  glm::vec3 returnColor = color;

  // iterate through each sphere
  for (int i = 0; i < num_spheres; i++){
    glm::vec3 intersect(0, 0, FURTHEST);

    bool checkIntersect = ray.checkIntersectSpheres(spheres[i], intersect);
    // FIXME: COUT STATEMENT
    std::cout << std::boolalpha << checkIntersect << " && " << intersect[2] << " && " << nearest << std::endl;
    if (checkIntersect && intersect[2] > nearest){
      // FIXME: COUT STATEMENT
      std::cout << "INTERSECTED SPHERE" << std::endl;
        returnColor = glm::vec3(0.0f);
        // now check  for each light
        for (int j = 0; j < num_lights; j++){
          bool light_on = true;
          std::cout << "test 1" << std::endl;
          glm::vec3 current_light_pos(lights[j].position[0], lights[j].position[1], lights[j].position[2]);
          Ray shadow_ray(intersect, glm::normalize(current_light_pos - intersect));

          // check shadow ray intersections with all other objects
          std::cout << "test 2" << std::endl;
          for (int k = 0; k < num_spheres; k++){
          // check for sphere intersections now, excluding current object
            std::cout << "test 3 -- entered loop" << std::endl;
            std::cout << num_spheres << std::endl;
            glm::vec3 shadow_intersect;
            if  (shadow_ray.checkIntersectSpheres(spheres[k], shadow_intersect) && k != i){
              std::cout << "test 4 -- checkedIntersection" << std::endl;
              if (glm::length(shadow_intersect - intersect) < glm::length(current_light_pos - intersect))
              {
                light_on = false;
                break;
              }
            }
            std::cout << "test 5 -- pastIntersection ifs" << std::endl;
          }

          // checks for triangle intersections second
          for (int k = 0; k < num_triangles; k++)
          {
            glm::vec3 shadow_intersect;
            if (shadow_ray.checkIntersectTriangles(triangles[k], shadow_intersect))
            {
              if (glm::length(shadow_intersect - intersect) < glm::length(current_light_pos - intersect))
              {
                light_on = false;
                break;
              }
            }
          }

          if (light_on){
            returnColor = returnColor + calculateLighting(spheres[i], lights[j], intersect);
          }
        }

      // new closest intersect
      nearest = intersect[2];
    }
  }

  // now check triangle collisions
  for (int i = 0; i < num_triangles; i++){
    glm::vec3 intersect(0, 0, FURTHEST);
    bool checkIntersect = ray.checkIntersectTriangles(triangles[i], intersect);
    // FIXME: COUT STATEMENT
    std::cout << std::boolalpha << checkIntersect  << " && " << (intersect[2] > nearest) << std::endl;
    if (checkIntersect && intersect[2] > nearest){
      // FIXME: COUT STATEMENT
        std::cout << "INTERSECTED TRIANGLE" << std::endl;
        returnColor = glm::vec3(0.0f);

        for (int j = 0; j < num_lights; j++){
          // new shadow ray starts from intersection point, towards the light source
          glm::vec3 direction = glm::vec3(lights[j].position[0], lights[j].position[1], lights[j].position[2]) - intersect;
          Ray shadow_ray(intersect, glm::normalize(direction));

          bool light_on = true;
          
          for (int k = 0; k < num_spheres; k++){
            glm::vec3 shadow_intersect(0, 0, FURTHEST);
            if (shadow_ray.checkIntersectSpheres(spheres[k], shadow_intersect)){
              if (glm::length(shadow_intersect - intersect) < glm::length(glm::vec3(lights[j].position[0], lights[j].position[1], lights[j].position[2]) - intersect)){
                light_on = false;
                break;
              }
            }
          }

          for (int k = 0; k < num_triangles; k++){
            glm::vec3 shadow_intersect(0, 0, FURTHEST);
            if (shadow_ray.checkIntersectTriangles(triangles[k], shadow_intersect) && k != i){
              if (glm::length(shadow_intersect - intersect) < glm::length(glm::vec3(lights[j].position[0], lights[j].position[1], lights[j].position[2]) - intersect))
              {
                light_on = false;
                break;
              }
            }
          }
          if (light_on)
          {
            returnColor = returnColor + calculateLighting(triangles[i], lights[j], intersect);
          }
        }

      //new closest intersect
      nearest = intersect[2];
    }
  }

  return returnColor;
}

glm::vec3 traceRay(Ray &ray){
  
  // set color to white, 
  glm::vec3 color(1.0f);
  double closest_intersect = FURTHEST; // set the nearest intersect

<<<<<<< HEAD
  // if (SOFT_SHADOWS){
  //   color = calculateCollisions(ray, color, closest_intersect, SOFT_SHADOWS);
  // }
  // else{
    color = calculateCollisions(ray, color, closest_intersect);
  // }
 
=======
  color = calculateCollisions(ray, color, closest_intersect);
  

>>>>>>> parent of 3b2bb8a... working main, broken soft shadows
  color += ambient;
  // std::cout << "2r " << color[0] << " g " << color[1] << " b " << color[2] << std::endl;

  return color;

}

std::vector<Ray> calculateCameraRays(double x, double y){
  // 4x supersampling
  // TODO: 8x Supersampling?

  std::vector<Ray> cameraRays;

  double xCamera = ((2 * (x + 0.25) / (double) WIDTH) - 1) * aspectRatio * angle;
  double yCamera = ((2 * (y + 0.25) / (double)HEIGHT) - 1) * angle;
  cameraRays.push_back(Ray(glm::vec3(0.0), glm::normalize(glm::vec3(xCamera, yCamera, -1))));

  xCamera = ((2 * (x + 0.25) / (double)WIDTH) - 1) * aspectRatio * angle;
  yCamera = ((2 * (y + 0.75) / (double)HEIGHT) - 1) * angle;
  cameraRays.push_back(Ray(glm::vec3(0.0), glm::normalize(glm::vec3(xCamera, yCamera, -1))));

  xCamera = ((2 * (x + 0.75) / (double)WIDTH) - 1) * aspectRatio * angle;
  yCamera = ((2 * (y + 0.75) / (double)HEIGHT) - 1) * angle;
  cameraRays.push_back(Ray(glm::vec3(0.0), glm::normalize(glm::vec3(xCamera, yCamera, -1))));

  xCamera = ((2 * (x + 0.75) / (double)WIDTH) - 1) * aspectRatio * angle;
  yCamera = ((2 * (y + 0.25) / (double)HEIGHT) - 1) * angle;
  cameraRays.push_back(Ray(glm::vec3(0.0), glm::normalize(glm::vec3(xCamera, yCamera, -1))));

  return cameraRays;
}

Ray calculateCameraRay(double x, double y){
  double xCamera = ((2 * (x + 0.5) / (double) WIDTH) - 1) * aspectRatio * angle;
  double yCamera = ((2 * (y + 0.5) / (double) HEIGHT) - 1) * angle;
  glm::vec3 origin(0.0);
  glm::vec3 direction(xCamera, yCamera, -1);
  Ray ray(origin, glm::normalize(direction));
  return ray;
}
//MODIFY THIS FUNCTION
void draw_scene()
{
  angle = tan((fov / 2) * PI / 180);
  // begin ray tracing
  for (unsigned int x = 0; x < WIDTH;  x++){

    glPointSize(2.0);
    glBegin(GL_POINTS);

    for (unsigned int y = 0; y < HEIGHT; y++){
      // FIXME: COUT STATEMENT
      std::cout << "x: " << x << " y: " << y << std::endl;
      Ray trace = calculateCameraRay(x, y);
      glm::vec3 color = traceRay(trace);
      color = glm::clamp(color, 0.0f, 1.0f);
      plot_pixel(x, y, color[0] * 255, color[1] * 255, color[2] * 255);
    }

    glEnd();
    glFlush();

  }
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);
  ambient = glm::vec3(ambient_light[0], ambient_light[1], ambient_light[2]);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
  
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0.0f,0.0f,0.0f,0.0f);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

