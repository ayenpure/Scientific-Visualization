#include "auxiliary.hxx"
#include "vtkDataSetReader.h"

template <typename T>
void GenerateRay(Vec3<T> &ray, unsigned int pixelX, unsigned int pixelY,
                 Vec3<T> look, Vec3<T> forX, Vec3<T> forY)
{
  T x_multiplier = (2*pixelX + 1 - WIDTH) / (T)2.0;
  T y_multiplier = (2*pixelY + 1 - HEIGHT) / (T)2.0;
  Vec3<T> dX = Multiply(forX, x_multiplier);
  Vec3<T> dY = Multiply(forY, y_multiplier);
  ray = look + dX + dY;
}

int main(int argc, char **argv) {

  if(argc < 3)
  {
    std::cout << "usage : raycaster <filename> <variable>" << std::endl;
    exit(EXIT_FAILURE);
  }

  const int height = HEIGHT;
  const int width = WIDTH;
  const std::string filename(argv[1]);
  int dims[3];

  // Read in the DataSet, which is a Rectilinear/Uniform Grid
  vtkDataSetReader *reader = vtkDataSetReader::New();
  reader->SetFileName(filename.c_str());
  reader->Update();
  vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
  rgrid->GetDimensions(dims);

  double *xCoords = (double*) rgrid->GetXCoordinates()->GetVoidPointer(0);
  double *yCoords = (double*) rgrid->GetYCoordinates()->GetVoidPointer(0);
  double *zCoords = (double*) rgrid->GetZCoordinates()->GetVoidPointer(0);
  double *fieldData = (double*) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

  Camera<double> camera = SetupCamera<double>();

  Vec3<double> look, forX, forY;
  CalculateViewParameters(camera, look, forX, forY);

  for (size_t pixelX = 0; pixelX < width; ++pixelX)
    for (size_t pixelY = 0; pixelY < height; ++pixelY) {
      Vec3<double> ray;
      GenerateRay(ray, pixelX, pixelY, look, forX, forY);
      Normalize(ray);
    }
}
