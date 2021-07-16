// Test program for quadratic-by-linear and cubic-by-linear patches

#include <cassert>
#include <fstream>
#include <iostream>

#include <dc.hh>

size_t n, m;                    // degrees : 2 <= n <= 3, m = 1
double Px[4][2], Py[4][2], Pz[4][2];

double f21(double x, double y, double z);
double f31(double x, double y, double z);
double (*fnm)(double, double, double);

void readSurface(std::string filename) {
  try {
    std::ifstream f(filename.c_str());
    f.exceptions(std::ios::failbit | std::ios::badbit);
    f >> n >> m;
    for (size_t i = 0; i <= n; ++i)
      for (size_t j = 0; j <= m; ++j)
        f >> Px[i][j] >> Py[i][j] >> Pz[i][j];
  } catch(std::ifstream::failure &) {
    std::cerr << "Error reading file" << std::endl;
    return;
  }
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <surface.bzr> [resolution] [scaling]" << std::endl;
    return 1;
  }

  readSurface(argv[1]);

  assert((n == 2 || n == 3) && m == 1);
  if (n == 2)
    fnm = f21;
  else
    fnm = f31;


  size_t resolution = 50;
  if (argc >= 3)
    resolution = std::atoi(argv[2]);

  double scaling = 1.0;
  if (argc >= 4)
    scaling = std::strtod(argv[3], nullptr);

  DualContouring::Point3D min = { Px[0][0], Py[0][0], Pz[0][0] }, max = min;
  for (size_t i = 0; i <= n; ++i)
    for (size_t j = 0; j <= m; ++j) {
      double c = Px[i][j];
      if (c < min[0])
        min[0] = c;
      if (c > max[0])
        max[0] = c;
    }
  for (size_t i = 0; i <= n; ++i)
    for (size_t j = 0; j <= m; ++j) {
      double c = Py[i][j];
      if (c < min[1])
        min[1] = c;
      if (c > max[1])
        max[1] = c;
    }
  for (size_t i = 0; i <= n; ++i)
    for (size_t j = 0; j <= m; ++j) {
      double c = Pz[i][j];
      if (c < min[2])
        min[2] = c;
      if (c > max[2])
        max[2] = c;
    }

  for (size_t i = 0; i < 3; ++i) {
    double med = (min[i] + max[i]) / 2;
    double len = (max[i] - min[i]) / 2;
    min[i] = med - len * scaling;
    max[i] = med + len * scaling;
  }

  auto mesh =
    DualContouring::isosurface([&](const DualContouring::Point3D &p) {
      return fnm(p[0], p[1], p[2]);
    }, 0, { min, max }, { resolution, resolution, resolution });

  mesh.writeOBJ("/tmp/implicitize-output.obj");
}
