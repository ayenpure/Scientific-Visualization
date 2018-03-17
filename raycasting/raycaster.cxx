#include "auxiliary.hxx"

#define HEIGHT 1000;
#define WIDTH  1000;

template<typename T>
void GenerateRay(Vec3<T>& ray, unsigned int pixelX, unsigned int pixelY)
{

}

template<typename T>
CalculateViewParameters(Camera<T>& camera, Vec3<T>& look, Vec3<T> forX, Vec3<T> forY)
{
  look 
}

int main(int argc, char** argv)
{
  const int height = HEIGHT;
  const int width  = WIDTH;

  Camera<double> camera = SetupCamera<double>();

  Vec3<T> look, forX, forY;
  CalculateViewParameters(camera, look, forX, forY);

  for(size_t pixelX; pixelX < width; ++pixelX)
    for(size_t pixelY; pixelY < height; ++pixelY)
    {
      Vec3<double> ray;
      GenerateRay(ray, pixelX, pixelY);
    }
}
