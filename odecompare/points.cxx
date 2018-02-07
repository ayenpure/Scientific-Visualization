#include <fstream>
#include <string>

int main(int argc, char** argv)
{
  std::string filename("points.txt");
  std::ofstream pointsFile;
  pointsFile.open(filename);
  float xVal = 4.0f;
  float diff = (4.0f - 2.0f)/49.0f;
  for(int i = 0; i < 50; i++)
  {
    pointsFile << xVal << " " << 0 << " " << 0 << std::endl;
    xVal -= diff;
  }
  pointsFile.close();
}
