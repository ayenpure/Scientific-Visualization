#include <algorithm>
#include <cmath>
#include <limits>

#include "auxiliary.hxx"
#include "vtkDataArray.h"
#include "vtkDataSetReader.h"
#include "vtkPointData.h"
#include "vtkRectilinearGrid.h"

template <typename T>
void GenerateRay(Vec3<T> &ray, int pixelX, int pixelY, Vec3<T> look,
                 Vec3<T> forX, Vec3<T> forY) {
  T x_multiplier = (2 * pixelX + 1 - WIDTH) / 2.0;
  T y_multiplier = (2 * pixelY + 1 - HEIGHT) / 2.0;
  Vec3<T> dX = Multiply(forX, x_multiplier);
  Vec3<T> dY = Multiply(forY, y_multiplier);
  ray = look + dX + dY;
}

template <typename T>
bool CheckRayVolumeIntersection(Vec3<T> &origin, Vec3<T> &direction,
                                T *bounds) {
  T txmin, txmax, tymin, tymax, tzmin, tzmax;
  Vec3<T> inverse = Inverse(direction);

  txmin = (bounds[0] - direction.x) * inverse.x;
  txmax = (bounds[1] - direction.x) * inverse.x;
  tymin = (bounds[2] - direction.y) * inverse.y;
  tymax = (bounds[3] - direction.y) * inverse.y;

  if ((txmin > tymax) || (tymin > txmax))
    return false;
  if (tymin > txmin)
    txmin = tymin;
  if (tymax < txmax)
    txmax = tymax;

  tzmin = (bounds[4] - direction.z) * inverse.z;
  tzmax = (bounds[5] - direction.z) * inverse.z;

  if ((txmin > tzmax) || (tzmin > txmax))
    return false;
  if (tzmin > txmin)
    txmin = tzmin;
  if (tzmax < txmax)
    txmax = tzmax;
  return true;
}

template <typename T>
T GetFieldValueForSample(Vec3<T> &currentSample, int *dims, T *bounds,
                         Vec3<T> &gridprop, T *fieldData) {

  // TODO : Check if point is in bounds.
  Vec3<T> translated = Vec3<T>((currentSample.x - bounds[0]) * gridprop.x,
                                 (currentSample.y - bounds[2]) * gridprop.y,
                                 (currentSample.z - bounds[4]) * gridprop.z);

  Vec3<int> vertex0 =
      Vec3<int>((int)translated.x, (int)translated.y, (int)translated.z);
  int index0 = GetPointIndex(vertex0, dims);
  Vec3<int> vectex1 = Vec3<int>(vertex0.x + 1,vertex0.y + 0,vertex0.z + 0);
  int index1 = GetPointIndex(vertex0, dims);
  Vec3<int> vectex2 = Vec3<int>(vertex0.x + 0,vertex0.y + 0,vertex0.z + 1);
  int index2 = GetPointIndex(vertex0, dims);
  Vec3<int> vectex3 = Vec3<int>(vertex0.x + 1,vertex0.y + 0,vertex0.z + 1);
  int index3 = GetPointIndex(vertex0, dims);
  Vec3<int> vectex4 = Vec3<int>(vertex0.x + 0,vertex0.y + 1,vertex0.z + 0);
  int index4 = GetPointIndex(vertex0, dims);
  Vec3<int> vectex5 = Vec3<int>(vertex0.x + 1,vertex0.y + 1,vertex0.z + 0);
  int index5 = GetPointIndex(vertex0, dims);
  Vec3<int> vectex6 = Vec3<int>(vertex0.x + 0,vertex0.y + 1,vertex0.z + 1);
  int index6 = GetPointIndex(vertex0, dims);
  Vec3<int> vectex7 = Vec3<int>(vertex0.x + 1,vertex0.y + 1,vertex0.z + 1);
  int index7 = GetPointIndex(vertex0, dims);

  T difference1 = translated.x - vertex0.x;
  T difference2 = translated.y - vertex0.y;
  T difference3 = translated.z - vertex0.z;

  //Interpolation in X (0,1), (2,3), (4,5), (6,7)
  T f01 = Interpolate(fieldData[index0], fieldData[index1], difference1);
  T f23 = Interpolate(fieldData[index2], fieldData[index3], difference1);
  T f45 = Interpolate(fieldData[index4], fieldData[index5], difference1);
  T f67 = Interpolate(fieldData[index6], fieldData[index7], difference1);

  //Interpolation in Y
  T f0123 = Interpolate(f01, f23, difference2);
  T f4567 = Interpolate(f45, f67, difference2);

  //Interpolaiton in Z
  T sampleValue = Interpolate(f0123, f4567, difference3);
  return sampleValue;
}

template <typename T>
void Sample(Camera<T> &camera, Vec3<T> &ray, const int samplerate,
            T samplingdiff, T *bounds, int *dims, Vec3<T> &gridprop,
            T *fieldData, Vec3<T> &pixelColor, TransferFunction &transfer) {

  double distance = camera.near;
  Vec3<T> origin = camera.position;

  PrintVector(origin);
  PrintVector(ray);
  std::cout << "Sampling step size : " << samplingdiff << std::endl;
  Vec3<T> currentSample;
  int count = 0;
  do {
    Vec3<T> offset = Multiply(ray, distance);
    currentSample = origin + offset;
    PrintVector(currentSample);
    T sampleValue = GetFieldValueForSample(currentSample, dims, bounds, gridprop, fieldData);
    std::cout << "field value" << sampleValue << std::endl;
    distance += samplingdiff;
    ++count;
  } while (count < 10);
}

int main(int argc, char **argv) {

  if (argc < 2) {
    std::cout << "usage : raycaster <filename> <variable>" << std::endl;
    exit(EXIT_FAILURE);
  }

  const int height = HEIGHT;
  const int width = WIDTH;
  const int samplerate = 256;
  const std::string filename(argv[1]);
  int dims[3];
  double bounds[6];

  // Read in the DataSet, which is a Rectilinear/Uniform Grid
  vtkDataSetReader *reader = vtkDataSetReader::New();
  reader->SetFileName(filename.c_str());
  reader->Update();
  vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *)reader->GetOutput();
  rgrid->GetDimensions(dims);
  rgrid->GetBounds(bounds);

  std::cout << "Dataset dimensions " << dims[0] << ", " << dims[1] << ", "
            << dims[2] << std::endl;

  double *xCoords = (double *)rgrid->GetXCoordinates()->GetVoidPointer(0);
  double *yCoords = (double *)rgrid->GetYCoordinates()->GetVoidPointer(0);
  double *zCoords = (double *)rgrid->GetZCoordinates()->GetVoidPointer(0);
  double *fieldData =
      (double *)rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

  Vec3<double> gridprop = Vec3<double>(dims[0] / (bounds[1] - bounds[0]),
                                       dims[1] / (bounds[3] - bounds[2]),
                                       dims[2] / (bounds[5] - bounds[4]));

  Camera<double> camera = SetupCamera<double>();
  Vec3<double> look, forX, forY;
  CalculateViewParameters(camera, look, forX, forY);
  std::cout << "Sampling rate : " << camera.far - camera.near << std::endl;
  double samplingdiff = (camera.far - camera.near) / (double)(samplerate - 1);

  TransferFunction transfer = SetupTransferFunction();

  long int intersections = 0;
  Vec3<double> ray;

  // RGB store for colors in the final image
  double *imagedata = new double[height * width * 3];

  for (size_t pixelX = 50; pixelX < 51 /*width*/; ++pixelX)
    for (size_t pixelY = 50; pixelY < 51 /*height*/; ++pixelY) {
      GenerateRay(ray, pixelX, pixelY, look, forX, forY);
      Normalize(ray);
      bool intersects =
          CheckRayVolumeIntersection(camera.position, ray, bounds);
      Vec3<double> pixelcolor(0, 0, 0);
      if (intersects) {
        std::cout << pixelX << ", " << pixelY << std::endl;
        Sample(camera, ray, samplerate, samplingdiff, bounds, dims, gridprop,
               fieldData, pixelcolor, transfer);
        ++intersections;
      }
      int pixelindex = (pixelY * width + pixelX) * 3;
      imagedata[pixelindex++] = pixelcolor.x;
      imagedata[pixelindex++] = pixelcolor.y;
      imagedata[pixelindex++] = pixelcolor.z;
    }
  std::cout << "Number of intersection : " << intersections << " out of "
            << height * width << " pixels" << std::endl;
}
