#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  std::ofstream dataset;
  std::string variable("velvec");
  std::string filename("analytical.vtk");

  dataset.open(filename);
  dataset << "# vtk DataFile Version 3.0" << std::endl;
  dataset << "vtk output" << std::endl;
  dataset << "ASCII" << std::endl;
  dataset << "DATASET RECTILINEAR_GRID" << std::endl;
  dataset << "DIMENSIONS 11 11 1" << std::endl;

 float currentPoint;
  currentPoint = -5.0f;
  dataset << "X_COORDINATES 11 float" << std::endl;
  while (currentPoint <= 5.0f) {
    dataset << currentPoint << " ";
    currentPoint += 1.0f;
  }
  dataset << std::endl;
  dataset << "Y_COORDINATES 11 float" << std::endl;
  currentPoint = -5.0f;
  while (currentPoint <= 5.0f) {
    dataset << currentPoint << " ";
    currentPoint += 1.0f;
  }
  dataset << std::endl;
  dataset << "Z_COORDINATES 1 float" << std::endl;
  dataset << 0 << std::endl;
  dataset << "CELL_DATA " << 10*10 << std::endl;
  dataset << "POINT_DATA " << 11*11 << std::endl;
  dataset << "VECTORS " << variable << " float" << std::endl;
  int count = 0;
  for (float x = -5.0f; x <= 5.0f; x+= 1.0f)
    for (float y = -5.0f; y <= 5.0f; y+= 1.0f) {
      float rVec[2] = {x, -y};
      float vel[2] = { -rVec[0], -rVec[1]};
      float magnitude = sqrt((vel[0]*vel[0]) + (vel[1]*vel[1]));
      float normal[2] = {0.0f, 0.0f};
      if(magnitude != 0.0f)
      {
        normal[0] = vel[0]/magnitude;
        normal[1] = vel[1]/magnitude;
      }
      std::cout << ++count << "(" << x << "," << y << ") :" << normal[0] << " " << normal[1] << " " << 0 << std::endl;
      dataset << normal[0] << " " << normal[1] << " " << 0 << std::endl;
    }
  dataset << std::endl;
  dataset.close();
  return 0;
}
