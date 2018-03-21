#include <algorithm>
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
bool CheckRayVolumeIntersection(Vec3<T> &origin, Vec3<T> &direction, T *bounds,
                                Vec3<T> near, Vec3<T> far) {
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

int main(int argc, char **argv) {

  if (argc < 2) {
    std::cout << "usage : raycaster <filename> <variable>" << std::endl;
    exit(EXIT_FAILURE);
  }

  const int height = HEIGHT;
  const int width = WIDTH;
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

  Camera<double> camera = SetupCamera<double>();

  Vec3<double> look, forX, forY;

  CalculateViewParameters(camera, look, forX, forY);
  long int intersections = 0;

  Vec3<double> ray;
  for (size_t pixelY = 0; pixelY < height; ++pixelY)
    for (size_t pixelX = 0; pixelX < width; ++pixelX) {
      GenerateRay(ray, pixelX, pixelY, look, forX, forY);
      Normalize(ray);
      bool intersects =
          CheckRayVolumeIntersection(camera.position, ray, bounds);
      if (intersects) {
        Sample(Camera, ray, bounds);
      }
    }
  std::cout << "Number of intersection : " << intersections << " out of "
            << height * width << " pixels" << std::endl;
}
