#include <algorithm>
#include <cmath>
#include <limits>

#include "auxiliary.hxx"
#include "vtkDataArray.h"
#include "vtkDataSetReader.h"
#include "vtkImageData.h"
#include "vtkPNGWriter.h"
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
                                double *bounds) {
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

template <typename T, typename U>
T GetFieldValueForSample(Vec3<T> &currentSample, int *dims, double *bounds,
                         Vec3<U> &gridprop, U *fieldData, U *xCoords,
                         U *yCoords, U *zCoords) {

  if ((currentSample.x < bounds[0] || currentSample.x > bounds[1]) ||
      (currentSample.y < bounds[2] || currentSample.y > bounds[3]) ||
      (currentSample.z < bounds[4] || currentSample.x > bounds[5]))
    return 0;

  Vec3<int> vertex0(0, 0, 0);
  for (int i = 0; i < dims[0] - 1; i++) {
    if (xCoords[i] <= currentSample.x && currentSample.x < xCoords[i + 1])
      vertex0.x = i;
  }
  for (int i = 0; i < dims[1] - 1; i++) {
    if (yCoords[i] <= currentSample.y && currentSample.y < yCoords[i + 1])
      vertex0.y = i;
  }
  for (int i = 0; i < dims[2] - 1; i++) {
    if (zCoords[i] <= currentSample.z && currentSample.z < zCoords[i + 1])
      vertex0.z = i;
  }

  //  std::cout << "Looking for point" << "\t";
  // PrintVector(currentSample);
  /*Vec3<int> vertex0 =
      Vec3<int>((int)translated.x, (int)translated.y, (int)translated.z);*/
  //  std::cout << "Found Cell" << "\t";
  // PrintVector(vertex0);

  int index0 = GetPointIndex(vertex0, dims);
  Vec3<int> vertex1 = Vec3<int>(vertex0.x + 1, vertex0.y + 0, vertex0.z + 0);
  int index1 = GetPointIndex(vertex1, dims);
  Vec3<int> vertex2 = Vec3<int>(vertex0.x + 0, vertex0.y + 0, vertex0.z + 1);
  int index2 = GetPointIndex(vertex2, dims);
  Vec3<int> vertex3 = Vec3<int>(vertex0.x + 1, vertex0.y + 0, vertex0.z + 1);
  int index3 = GetPointIndex(vertex3, dims);
  Vec3<int> vertex4 = Vec3<int>(vertex0.x + 0, vertex0.y + 1, vertex0.z + 0);
  int index4 = GetPointIndex(vertex4, dims);
  Vec3<int> vertex5 = Vec3<int>(vertex0.x + 1, vertex0.y + 1, vertex0.z + 0);
  int index5 = GetPointIndex(vertex5, dims);
  Vec3<int> vertex6 = Vec3<int>(vertex0.x + 0, vertex0.y + 1, vertex0.z + 1);
  int index6 = GetPointIndex(vertex6, dims);
  Vec3<int> vertex7 = Vec3<int>(vertex0.x + 1, vertex0.y + 1, vertex0.z + 1);
  int index7 = GetPointIndex(vertex7, dims);

  U proportionX = (currentSample.x - xCoords[vertex0.x]) /
                  (xCoords[vertex0.x + 1] - xCoords[vertex0.x]);
  U proportionY = (currentSample.y - yCoords[vertex0.y]) /
                  (yCoords[vertex0.y + 1] - yCoords[vertex0.y]);
  U proportionZ = (currentSample.z - zCoords[vertex0.z]) /
                  (zCoords[vertex0.z + 1] - zCoords[vertex0.z]);

  // Interpolation in X (0,1), (2,3), (4,5), (6,7)
  U f01 = Interpolate(fieldData[index0], fieldData[index1], proportionX);
  U f23 = Interpolate(fieldData[index2], fieldData[index3], proportionX);
  U f45 = Interpolate(fieldData[index4], fieldData[index5], proportionX);
  U f67 = Interpolate(fieldData[index6], fieldData[index7], proportionX);

  // Interpolation in Z
  U f0123 = Interpolate(f01, f23, proportionZ);
  U f4567 = Interpolate(f45, f67, proportionZ);

  // Interpolaiton in Y
  U sampleValue = Interpolate(f0123, f4567, proportionY);
  return (T)sampleValue;
}

template <typename T, typename U>
void Sample(Camera<T> &camera, Vec3<T> &ray, const int samples, T samplingdiff,
            double *bounds, int *dims, Vec3<U> &gridprop, U *fieldData,
            Vec3<unsigned char> &pixelColor, TransferFunction &transfer,
            U *xCoords, U *yCoords, U *zCoords) {

  double distance = camera.near;
  Vec3<T> origin = camera.position;

  // PrintVector(origin);
  // PrintVector(ray);
  // std::cout << "Sampling step size : " << samplingdiff << std::endl;
  Vec3<T> currentSample;
  int count = 0;
  Vec3<double> color(0, 0, 0);
  T opacity = 0;
  do {
    //std::cout << "********************* sample :" << count << "*********************" << std::endl;
    Vec3<T> offset = Multiply(ray, distance);
    currentSample = origin + offset;
    //PrintVector(currentSample);
    T sampleValue =
        GetFieldValueForSample(currentSample, dims, bounds, gridprop, fieldData,
                               xCoords, yCoords, zCoords);
    //std::cout << "field value : " << sampleValue << std::endl;
    T backopacity = 0;
    Vec3<unsigned char> backcolor(0, 0, 0);
    transfer.ApplyTransferFunction(sampleValue, backcolor, backopacity);
/*
    std::cout << "mapped color : " << (int)backcolor.x << ", "
              << (int)backcolor.y << ", " << (int)backcolor.z << std::endl;
    std::cout << "mapped opacity : " << backopacity << std::endl;
*/
    backopacity = 1 - pow((1 - backopacity), 500 / (T)samples);
    //std::cout << "Corrected opacity : " << backopacity << std::endl;

    // Blending values
    color.x =
        color.x + (1 - opacity) * backopacity * (backcolor.x / 255.0);
    color.y =
        color.y + (1 - opacity) * backopacity * (backcolor.y / 255.0);
    color.z =
        color.z + (1 - opacity) * backopacity * (backcolor.z / 255.0);
    opacity = opacity + (1 - opacity) * backopacity;

    /*color.x =
        backopacity * backcolor.x / 255.0 + (1 - backopacity) * opacity * (color.x);
    color.x =
        backopacity * backcolor.y / 255.0 + (1 - backopacity) * opacity * (color.y);
    color.x =
        backopacity * backcolor.z / 255.0 + (1 - backopacity) * opacity * (color.z);
    opacity = backopacity + (1 - backopacity) * opacity;
*/

/*    std::cout << "Running : " << color.x << ", " << color.y << ", " << color.z
              << std::endl;
    std::cout << "Opacity : " << opacity << std::endl;
*/
    distance += samplingdiff;
    ++count;
  } while (count < samples);
  pixelColor.x = (unsigned char)(color.x * 255);
  pixelColor.y = (unsigned char)(color.y * 255);
  pixelColor.z = (unsigned char)(color.z * 255);
/*  std::cout << "Final : " << (int)pixelColor.x << ", " << (int)pixelColor.y << ", " << (int)pixelColor.z
              << std::endl; */
}

int main(int argc, char **argv) {

  if (argc < 2) {
    std::cout << "usage : raycaster <filename> <variable>" << std::endl;
    exit(EXIT_FAILURE);
  }

  const int height = HEIGHT;
  const int width = WIDTH;
  const int samples = SAMPLES;
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

  float *xCoords = (float *)rgrid->GetXCoordinates()->GetVoidPointer(0);
  float *yCoords = (float *)rgrid->GetYCoordinates()->GetVoidPointer(0);
  float *zCoords = (float *)rgrid->GetZCoordinates()->GetVoidPointer(0);
  float *fieldData =
      (float *)rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

  Vec3<float> gridprop =
      Vec3<float>(dims[0] / (xCoords[dims[0] - 1] - xCoords[0]),
                  dims[1] / (yCoords[dims[1] - 1] - yCoords[0]),
                  dims[2] / (zCoords[dims[2] - 1] - xCoords[0]));

  Camera<double> camera = SetupCamera<double>();
  Vec3<double> look, forX, forY;
  CalculateViewParameters(camera, look, forX, forY);
  // std::cout << "Sampling rate : " << camera.far - camera.near << std::endl;
  double samplingdiff = (camera.far - camera.near) / (double)(samples - 1);

  TransferFunction transfer = SetupTransferFunction();

  long int intersections = 0;
  Vec3<double> ray;

  vtkImageData *image = vtkImageData::New();
  image->SetDimensions(width, height, 1);
  image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
  unsigned char *imagedata = (unsigned char *)image->GetScalarPointer(0, 0, 0);

  for (size_t pixelX = 0; pixelX < width; ++pixelX)
    for (size_t pixelY = 0; pixelY < height; ++pixelY) {
      GenerateRay(ray, pixelX, pixelY, look, forX, forY);
      Normalize(ray);
      /*bool intersects =
          CheckRayVolumeIntersection(camera.position, ray, bounds);*/
      Vec3<unsigned char> pixelcolor(0, 0, 0);
      // if (intersects) {
      // std::cout << pixelX << ", " << pixelY << std::endl;
      Sample(camera, ray, samples, samplingdiff, bounds, dims, gridprop,
             fieldData, pixelcolor, transfer, xCoords, yCoords, zCoords);
      //++intersections;
      //}
      int pixelindex = (pixelY * width + pixelX) * 3;
      imagedata[pixelindex++] = pixelcolor.x;
      imagedata[pixelindex++] = pixelcolor.y;
      imagedata[pixelindex++] = pixelcolor.z;
    }
  /*std::cout << "Number of intersection : " << intersections << " out of "
            << height * width << " pixels" << std::endl;*/

  const std::string outfilename("output.png");
  vtkPNGWriter *writer = vtkPNGWriter::New();
  writer->SetInputData(image);
  writer->SetFileName(outfilename.c_str());
  writer->Write();
  writer->Delete();
}
