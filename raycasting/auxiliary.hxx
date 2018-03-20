#include <cmath>
#include <iostream>

#define HEIGHT 100
#define WIDTH 100
#define PI 3.14159265

template <typename T>
struct Vec3 {
  T x, y, z;

  Vec3(){};
  Vec3(T xVal, T yVal, T zVal) : x(xVal), y(yVal), z(zVal) {}

  Vec3 operator -(Vec3& vector)
  {
    return Vec3(x - vector.x, y - vector.y, z - vector.z);
  }

  Vec3 operator +(Vec3& vector)
  {
    return Vec3(x + vector.x, y + vector.y, z + vector.z);
  }

};

template <typename T>
void PrintVector(Vec3<T> &vector) {
  std::cout << "{" << vector.x << ", " << vector.y << ", " << vector.z << "}"
            << std::endl;
}

template <typename T>
Vec3<T> Multiply(Vec3<T> &vector, T val) {
  Vec3<T> output = Vec3<T>(vector.x * val,vector.y * val,vector.z * val);
  return output;
}

template <typename T>
T Magnitude(Vec3<T> &vector) {
  T det = vector.x * vector.x + vector.y * vector.y + vector.z * vector.z;
  T magnitude = sqrt(det);
  return magnitude;
}

template <typename T>
void Normalize(Vec3<T> &vector) {
  T magnitude = Magnitude(vector);
  vector.x /= magnitude;
  vector.y /= magnitude;
  vector.z /= magnitude;
}

template <typename T>
Vec3<T> Cross(Vec3<T> &input1, Vec3<T> &input2) {
  Vec3<T> output = Vec3<T>(input1.y * input2.z - input1.z * input2.y,
                           input1.x * input2.z - input1.z * input2.x,
                           input1.x * input2.y - input1.y * input2.x);
  return output;
}

template <typename T>
struct Camera {
  T near, far;
  T angle;
  Vec3<T> position;
  Vec3<T> focus;
  Vec3<T> up;
};

template <typename T>
Camera<T> SetupCamera(void) {
  Camera<T> camera;

  camera.focus = Vec3<T>(0, 0, 0);
  camera.up = Vec3<T>(0, 1, 0);
  camera.position = Vec3<T>(-8.25e+7, -3.45e+7, 3.35e+7);

  camera.angle = 30;
  camera.near = 7.5e+7;
  camera.far = 1.4e+8;

  return camera;
}

template <typename T>
void CalculateViewParameters(Camera<T> &camera, Vec3<T> &look, Vec3<T> &forX,
                        Vec3<T> &forY) {
  look = Vec3<T>(camera.focus - camera.position);
  Normalize(look);
  Vec3<T> u = Cross(look, camera.up);
  Normalize(u);
  Vec3<T> v = Cross(look, u);
  Normalize(v);
  T radians = camera.angle * PI / 180.0;
  T x_multiplier = 2*tan(radians)/ (T)WIDTH;
  T y_multiplier = 2*tan(radians)/ (T)HEIGHT;
  forX = Multiply(u, x_multiplier);
  forY = Multiply(v, y_multiplier);
}
