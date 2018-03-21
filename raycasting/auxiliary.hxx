#include <cmath>
#include <iostream>

#define HEIGHT 100
#define WIDTH 100
#define SAMPLES 256

#define PI 3.14159265

template <typename T, typename U>
T Interpolate(U fieldValue1, U fieldValue2, T proportion)
{
  //return (1.0f - proportion)*fieldvalue1 + proportion*fieldvalue2;
  return fieldValue1 + proportion*(fieldValue2 - fieldValue1);
}

template <typename T> struct Vec3 {
  T x, y, z;

  Vec3(){};
  Vec3(T xVal, T yVal, T zVal) : x(xVal), y(yVal), z(zVal) {}

  Vec3 operator-(Vec3 &vector) {
    return Vec3(x - vector.x, y - vector.y, z - vector.z);
  }

  Vec3 operator+(Vec3 &vector) {
    return Vec3(x + vector.x, y + vector.y, z + vector.z);
  }
};

template <typename T> void PrintVector(Vec3<T> &vector) {
  std::cout << "{" << vector.x << ", " << vector.y << ", " << vector.z << "}"
            << std::endl;
}

template <typename T> Vec3<T> Multiply(Vec3<T> &vector, T val) {
  Vec3<T> output = Vec3<T>(vector.x * val, vector.y * val, vector.z * val);
  return output;
}

template <typename T> T Magnitude(Vec3<T> &vector) {
  T det = vector.x * vector.x + vector.y * vector.y + vector.z * vector.z;
  T magnitude = sqrt(det);
  return magnitude;
}

template <typename T> void Normalize(Vec3<T> &vector) {
  T magnitude = Magnitude(vector);
  vector.x /= magnitude;
  vector.y /= magnitude;
  vector.z /= magnitude;
}

template <typename T> Vec3<T> Inverse(Vec3<T> &vector) {
  T det = vector.x * vector.x + vector.y * vector.y + vector.z * vector.z;
  Vec3<T> inverse = Vec3<T>(vector.x / det, vector.y / det, vector.z / det);
  return inverse;
}

template <typename T> Vec3<T> Cross(Vec3<T> &input1, Vec3<T> &input2) {
  Vec3<T> output = Vec3<T>(input1.y * input2.z - input1.z * input2.y,
                           input1.z * input2.x - input1.x * input2.z,
                           input1.x * input2.y - input1.y * input2.x);
  return output;
}

int GetPointIndex(Vec3<int>& index, const int *dims)
{
  int sliceSize = (dims[0])*(dims[1]);
  int width = (dims[0]);
  return index.z*sliceSize + index.y*width + index.x;
}

template <typename T> struct Camera {
  T near, far;
  T angle;
  Vec3<T> position;
  Vec3<T> focus;
  Vec3<T> up;
};

template <typename T> Camera<T> SetupCamera(void) {
  Camera<T> camera;

  camera.focus = Vec3<T>(0, 0, 0);
  camera.up = Vec3<T>(0, -1, 0);
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
  T x_multiplier = 2 * tan(radians / 2.0) / (T)WIDTH;
  T y_multiplier = 2 * tan(radians / 2.0) / (T)HEIGHT;
  forX = Multiply(u, x_multiplier);
  forY = Multiply(v, y_multiplier);
}

struct TransferFunction {
  double min;
  double max;
  int numBins;
  unsigned char *colors; // size is 3*numBins
  double *opacities;     // size is numBins
  double scalingFactor;
  // Take in a value and applies the transfer function.
  // Step #1: figure out which bin "value" lies in.
  // If "min" is 2 and "max" is 4, and there are 10 bins, then
  //   bin 0 = 2->2.2
  //   bin 1 = 2.2->2.4
  //   bin 2 = 2.4->2.6
  //   bin 3 = 2.6->2.8
  //   bin 4 = 2.8->3.0
  //   bin 5 = 3.0->3.2
  //   bin 6 = 3.2->3.4
  //   bin 7 = 3.4->3.6
  //   bin 8 = 3.6->3.8
  //   bin 9 = 3.8->4.0
  // and, for example, a "value" of 3.15 would return the color in bin 5
  // and the opacity at "opacities[5]".
  template<typename T, typename U>
  void ApplyTransferFunction(T value, Vec3<U>& RGB,
                             T &opacity) {
     if(value < min || value > max)
       return;
     int bin = (value - min) * scalingFactor;
     RGB.x = colors[3*bin+0];
     RGB.y = colors[3*bin+1];
     RGB.z = colors[3*bin+2];
     opacity = opacities[bin];
  }
};

TransferFunction SetupTransferFunction(void) {
  int i;

  TransferFunction rv;
  rv.min = 10;
  rv.max = 15;
  rv.numBins = 256;
  rv.colors = new unsigned char[3 * 256];
  rv.opacities = new double[256];
  unsigned char charOpacity[256] = {
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   1,   1,   1,   2,   2,   3,   3,   4,   5,   6,   7,   8,
      9,   10,  11,  12,  13,  13,  14,  14,  14,  14,  14,  14,  14,  13,  12,
      11,  10,  9,   8,   7,   6,   5,   5,   4,   3,   2,   3,   3,   4,   5,
      6,   7,   8,   9,   10,  11,  12,  13,  14,  15,  16,  17,  17,  17,  17,
      17,  17,  17,  16,  16,  15,  14,  13,  12,  11,  9,   8,   7,   6,   5,
      5,   4,   3,   3,   3,   4,   5,   6,   7,   8,   9,   11,  12,  14,  16,
      18,  20,  22,  24,  27,  29,  32,  35,  38,  41,  44,  47,  50,  52,  55,
      58,  60,  62,  64,  66,  67,  68,  69,  70,  70,  70,  69,  68,  67,  66,
      64,  62,  60,  58,  55,  52,  50,  47,  44,  41,  38,  35,  32,  29,  27,
      24,  22,  20,  20,  23,  28,  33,  38,  45,  51,  59,  67,  76,  85,  95,
      105, 116, 127, 138, 149, 160, 170, 180, 189, 198, 205, 212, 217, 221, 223,
      224, 224, 222, 219, 214, 208, 201, 193, 184, 174, 164, 153, 142, 131, 120,
      109, 99,  89,  79,  70,  62,  54,  47,  40,  35,  30,  25,  21,  17,  14,
      12,  10,  8,   6,   5,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0};

  for (i = 0; i < 256; i++)
    rv.opacities[i] = charOpacity[i] / 255.0;
  const int numControlPoints = 8;
  unsigned char controlPointColors[numControlPoints * 3] = {
      71,  71,  219, 0,   0,  91, 0,   255, 255, 0,   127, 0,
      255, 255, 0,   255, 96, 0,  107, 0,   0,   224, 76,  76};
  double controlPointPositions[numControlPoints] = {0,     0.143, 0.285, 0.429,
                                                    0.571, 0.714, 0.857, 1.0};
  for (i = 0; i < numControlPoints - 1; i++) {
    int start = controlPointPositions[i] * rv.numBins;
    int end = controlPointPositions[i + 1] * rv.numBins + 1;
    /*cerr << "Working on " << i << "/" << i + 1 << ", with range " << start
         << "/" << end << endl;*/
    if (end >= rv.numBins)
      end = rv.numBins - 1;
    for (int j = start; j <= end; j++) {
      double proportion =
          (j / (rv.numBins - 1.0) - controlPointPositions[i]) /
          (controlPointPositions[i + 1] - controlPointPositions[i]);
      if (proportion < 0 || proportion > 1.)
        continue;
      for (int k = 0; k < 3; k++)
        rv.colors[3 * j + k] =
            proportion * (controlPointColors[3 * (i + 1) + k] -
                          controlPointColors[3 * i + k]) +
            controlPointColors[3 * i + k];
    }
  }
  rv.scalingFactor = 256 / (rv.max - rv.min);
  return rv;
}
