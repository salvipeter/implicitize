#include <fstream>
#include <iostream>

#include <geometry.hh>

double Px[4], Py[4];
extern "C" {
  double cubic(double x, double y);
  double cubic_x(double x, double y);
  double cubic_y(double x, double y);
}

using namespace Geometry;

int main(int argc, char **argv) {
  Px[0] = 0; Py[0] = 0;
  Px[1] = 1; Py[1] = 2;
  Px[2] = 2; Py[2] = 3;
  Px[3] = 4; Py[3] = 0;

  // Write curve
  {
    BSCurve curve({{Px[0],Py[0],0},{Px[1],Py[1],0},{Px[2],Py[2],0},{Px[3],Py[3],0}});
    std::ofstream f("/tmp/curve.obj");
    for (size_t i = 0; i < 100; ++i)
      f << "v " << curve.eval(i / 99.0) << std::endl;
    f << 'l';
    for (size_t i = 1; i <= 100; ++i)
      f << ' ' << i;
    f << std::endl;
  }

  // Write SDF
  {
    size_t res = 200;
    double min_x = -2, min_y = -2, max_x = 7, max_y = 7;
    std::ofstream f("/tmp/sdf.vtk");
    f << "# vtk DataFile Version 2.0" << std::endl;
    f << "Implicitized Bezier curve" << std::endl;
    f << "ASCII" << std::endl;
    f << "DATASET POLYDATA" << std::endl;
    f << "POINTS " << res * res << " float" << std::endl;
    for (size_t i = 0; i < res; ++i) {
      double x = min_x + (max_x - min_x) * i / (res - 1);
      for (size_t j = 0; j < res; ++j) {
        double y = min_y + (max_y - min_y) * j / (res - 1);
        f << x << ' ' << y << ' ' << 0 << std::endl;
      }
    }
    f << "POLYGONS " << (res - 1) * (res - 1) << ' ' << (res - 1) * (res - 1) * 5 << std::endl;
    for (size_t i = 1; i < res; ++i)
      for (size_t j = 1; j < res; ++j)
        f << "4 "
          << (i - 1) * res + j - 1 << ' '
          << (i - 1) * res + j << ' '
          << i * res + j << ' '
          << i * res + j - 1 << std::endl;
    f << "POINT_DATA " << res * res << std::endl;
    f << "SCALARS distance float 1" << std::endl;
    f << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < res; ++i) {
      double x = min_x + (max_x - min_x) * i / (res - 1);
      for (size_t j = 0; j < res; ++j) {
        double y = min_y + (max_y - min_y) * j / (res - 1);
        auto dx = cubic_x(x, y);
        auto dy = cubic_y(x, y);
        auto gradnorm = Vector2D(dx, dy).norm();
        f << cubic(x, y) / gradnorm << std::endl;
      }
    }
    f << "NORMALS gradient float" << std::endl;
    for (size_t i = 0; i < res; ++i) {
      double x = min_x + (max_x - min_x) * i / (res - 1);
      for (size_t j = 0; j < res; ++j) {
        double y = min_y + (max_y - min_y) * j / (res - 1);
        auto dx = cubic_x(x, y);
        auto dy = cubic_y(x, y);
        f << Vector3D(dx, dy, 0).normalize() << std::endl;
      }
    }
  }
}
