#include <math.h>

extern double Px[4], Py[4];

double cubic(double x, double y) {
  return Px[3] *
             (Px[2] *
                  (Px[1] *
                       (54 * pow(y, 3) +
                        (27 * Py[3] - 81 * Py[2] - 81 * Py[1] - 27 * Py[0]) *
                            pow(y, 2) +
                        (-(27 * pow(Py[0], 2)) - 54 * Py[0] * Py[3] +
                         (81 * Py[1] + 81 * Py[0]) * Py[2] +
                         81 * Py[0] * Py[1]) *
                            y +
                        27 * Py[3] * pow(Py[0], 2) -
                        81 * Py[0] * Py[1] * Py[2]) +
                   Px[0] *
                       (-(18 * pow(y, 3)) +
                        (-(9 * Py[3]) - 27 * Py[2] + 153 * Py[1] - 63 * Py[0]) *
                            pow(y, 2) +
                        (54 * pow(Py[2], 2) - 81 * pow(Py[1], 2) +
                         (9 * Py[1] + 9 * Py[0]) * Py[3] +
                         (108 * Py[0] - 162 * Py[1]) * Py[2] +
                         9 * Py[0] * Py[1]) *
                            y -
                        54 * Py[0] * pow(Py[2], 2) +
                        81 * Py[2] * pow(Py[1], 2) -
                        9 * Py[0] * Py[1] * Py[3])) +
              Px[0] * Px[1] *
                  (18 * pow(y, 3) +
                   (63 * Py[3] - 153 * Py[2] + 27 * Py[1] + 9 * Py[0]) *
                       pow(y, 2) +
                   (81 * pow(Py[2], 2) - 54 * pow(Py[1], 2) +
                    (-(9 * Py[2]) - 108 * Py[1] - 9 * Py[0]) * Py[3] +
                    (162 * Py[1] - 9 * Py[0]) * Py[2]) *
                       y -
                   81 * Py[1] * pow(Py[2], 2) +
                   Py[3] * (54 * pow(Py[1], 2) + 9 * Py[0] * Py[2])) +
              pow(Px[0], 2) *
                  (-(3 * pow(y, 3)) +
                   (-(21 * Py[3]) + 54 * Py[2] - 27 * Py[1] + 3 * Py[0]) *
                       pow(y, 2) +
                   (-(3 * pow(Py[3], 2)) - 81 * pow(Py[2], 2) +
                    (27 * Py[2] + 27 * Py[1] - 6 * Py[0]) * Py[3] +
                    27 * Py[1] * Py[2]) *
                       y +
                   3 * Py[0] * pow(Py[3], 2) + 27 * pow(Py[2], 3) -
                   27 * Py[1] * Py[2] * Py[3]) +
              pow(Px[1], 2) *
                  (-(27 * pow(y, 3)) +
                   (-(54 * Py[3]) + 162 * Py[2] - 54 * Py[1] + 27 * Py[0]) *
                       pow(y, 2) +
                   (-(81 * pow(Py[2], 2)) + (54 * Py[1] + 54 * Py[0]) * Py[3] -
                    162 * Py[0] * Py[2] + 54 * Py[0] * Py[1]) *
                       y +
                   81 * Py[0] * pow(Py[2], 2) - 54 * Py[0] * Py[1] * Py[3]) +
              pow(Px[2], 2) *
                  (-(27 * pow(y, 3)) + (27 * Py[2] + 54 * Py[0]) * pow(y, 2) +
                   (-(27 * pow(Py[0], 2)) - 54 * Py[0] * Py[2]) * y +
                   27 * Py[2] * pow(Py[0], 2))) +
         Px[2] *
             (pow(Px[1], 2) *
                  (81 * pow(y, 3) +
                   (-(81 * Py[3]) - 81 * Py[2] - 81 * Py[0]) * pow(y, 2) +
                   ((81 * Py[2] + 81 * Py[0]) * Py[3] + 81 * Py[0] * Py[2]) *
                       y -
                   81 * Py[0] * Py[2] * Py[3]) +
              pow(Px[0], 2) *
                  (9 * pow(y, 3) +
                   (9 * Py[3] - 54 * Py[2] + 18 * Py[1]) * pow(y, 2) +
                   (-(18 * pow(Py[3], 2)) + 27 * pow(Py[2], 2) +
                    (54 * Py[2] - 36 * Py[1]) * Py[3]) *
                       y +
                   18 * Py[1] * pow(Py[3], 2) - 27 * Py[3] * pow(Py[2], 2)) +
              Px[0] * Px[1] *
                  (-(54 * pow(y, 3)) +
                   (27 * Py[3] + 81 * Py[2] + 81 * Py[1] - 27 * Py[0]) *
                       pow(y, 2) +
                   (27 * pow(Py[3], 2) +
                    (-(81 * Py[2]) - 81 * Py[1] + 54 * Py[0]) * Py[3] -
                    81 * Py[1] * Py[2]) *
                       y -
                   27 * Py[0] * pow(Py[3], 2) + 81 * Py[1] * Py[2] * Py[3])) +
         pow(Px[2], 2) *
             (Px[0] * (27 * pow(y, 3) +
                       (-(27 * Py[3]) + 54 * Py[2] - 162 * Py[1] + 54 * Py[0]) *
                           pow(y, 2) +
                       (81 * pow(Py[1], 2) +
                        (-(54 * Py[2]) + 162 * Py[1] - 54 * Py[0]) * Py[3] -
                        54 * Py[0] * Py[2]) *
                           y +
                       Py[3] * (54 * Py[0] * Py[2] - 81 * pow(Py[1], 2))) +
              Px[1] *
                  (-(81 * pow(y, 3)) +
                   (81 * Py[3] + 81 * Py[1] + 81 * Py[0]) * pow(y, 2) +
                   ((-(81 * Py[1]) - 81 * Py[0]) * Py[3] - 81 * Py[0] * Py[1]) *
                       y +
                   81 * Py[0] * Py[1] * Py[3])) +
         pow(Px[3], 2) *
             (Px[2] * (9 * pow(y, 3) + (-(9 * Py[1]) - 18 * Py[0]) * pow(y, 2) +
                       (9 * pow(Py[0], 2) + 18 * Py[0] * Py[1]) * y -
                       9 * Py[1] * pow(Py[0], 2)) +
              Px[0] *
                  (3 * pow(y, 3) +
                   (-(3 * Py[3]) + 27 * Py[2] - 54 * Py[1] + 21 * Py[0]) *
                       pow(y, 2) +
                   (81 * pow(Py[1], 2) + 3 * pow(Py[0], 2) + 6 * Py[0] * Py[3] +
                    (-(27 * Py[1]) - 27 * Py[0]) * Py[2] - 27 * Py[0] * Py[1]) *
                       y -
                   27 * pow(Py[1], 3) - 3 * Py[3] * pow(Py[0], 2) +
                   27 * Py[0] * Py[1] * Py[2]) +
              Px[1] *
                  (-(9 * pow(y, 3)) +
                   (-(18 * Py[2]) + 54 * Py[1] - 9 * Py[0]) * pow(y, 2) +
                   (-(27 * pow(Py[1], 2)) + 18 * pow(Py[0], 2) +
                    36 * Py[0] * Py[2] - 54 * Py[0] * Py[1]) *
                       y +
                   27 * Py[0] * pow(Py[1], 2) - 18 * Py[2] * pow(Py[0], 2))) +
         pow(Px[2], 3) *
             (27 * pow(y, 3) + (-(27 * Py[3]) - 54 * Py[0]) * pow(y, 2) +
              (27 * pow(Py[0], 2) + 54 * Py[0] * Py[3]) * y -
              27 * Py[3] * pow(Py[0], 2)) +
         Px[0] * pow(Px[1], 2) *
             (27 * pow(y, 3) + (-(54 * Py[3]) - 27 * Py[1]) * pow(y, 2) +
              (27 * pow(Py[3], 2) + 54 * Py[1] * Py[3]) * y -
              27 * Py[1] * pow(Py[3], 2)) +
         pow(Px[0], 3) * (pow(y, 3) - 3 * Py[3] * pow(y, 2) +
                          3 * pow(Py[3], 2) * y - pow(Py[3], 3)) +
         pow(Px[3], 3) * (-pow(y, 3) + 3 * Py[0] * pow(y, 2) -
                          3 * pow(Py[0], 2) * y + pow(Py[0], 3)) +
         Px[1] * pow(Px[0], 2) *
             (-(9 * pow(y, 3)) + (18 * Py[3] + 9 * Py[2]) * pow(y, 2) +
              (-(9 * pow(Py[3], 2)) - 18 * Py[2] * Py[3]) * y +
              9 * Py[2] * pow(Py[3], 2)) +
         pow(Px[1], 3) *
             (-(27 * pow(y, 3)) + (54 * Py[3] + 27 * Py[0]) * pow(y, 2) +
              (-(27 * pow(Py[3], 2)) - 54 * Py[0] * Py[3]) * y +
              27 * Py[0] * pow(Py[3], 2)) +
         x * (Px[2] *
                  (Px[0] *
                       ((18 * Py[3] - 54 * Py[2] + 54 * Py[1] - 18 * Py[0]) *
                            pow(y, 2) +
                        (9 * pow(Py[3], 2) - 108 * pow(Py[2], 2) -
                         81 * pow(Py[1], 2) +
                         (81 * Py[2] - 180 * Py[1] + 45 * Py[0]) * Py[3] +
                         243 * Py[1] * Py[2] - 9 * Py[0] * Py[1]) *
                            y +
                        (27 * Py[0] - 36 * Py[1]) * pow(Py[3], 2) +
                        Py[3] * (54 * pow(Py[2], 2) + 162 * pow(Py[1], 2) +
                                 (-(81 * Py[1]) - 108 * Py[0]) * Py[2] +
                                 9 * Py[0] * Py[1]) +
                        54 * Py[0] * pow(Py[2], 2) -
                        81 * Py[2] * pow(Py[1], 2)) +
                   Px[1] * ((-(54 * Py[3]) + 162 * Py[2] - 162 * Py[1] +
                             54 * Py[0]) *
                                pow(y, 2) +
                            (-(27 * pow(Py[3], 2)) + 27 * pow(Py[0], 2) +
                             (243 * Py[1] - 81 * Py[2]) * Py[3] -
                             243 * Py[0] * Py[2] + 81 * Py[0] * Py[1]) *
                                y +
                            27 * Py[0] * pow(Py[3], 2) +
                            Py[3] * (-(27 * pow(Py[0], 2)) +
                                     (162 * Py[0] - 81 * Py[1]) * Py[2] -
                                     162 * Py[0] * Py[1]) +
                            81 * Py[0] * Py[1] * Py[2])) +
              Px[3] *
                  (Px[1] *
                       ((18 * Py[3] - 54 * Py[2] + 54 * Py[1] - 18 * Py[0]) *
                            pow(y, 2) +
                        (81 * pow(Py[2], 2) + 108 * pow(Py[1], 2) -
                         9 * pow(Py[0], 2) + (9 * Py[2] - 45 * Py[0]) * Py[3] +
                         (180 * Py[0] - 243 * Py[1]) * Py[2] -
                         81 * Py[0] * Py[1]) *
                            y +
                        (81 * Py[1] - 162 * Py[0]) * pow(Py[2], 2) -
                        54 * Py[0] * pow(Py[1], 2) +
                        Py[3] * (-(54 * pow(Py[1], 2)) - 27 * pow(Py[0], 2) -
                                 9 * Py[0] * Py[2] + 108 * Py[0] * Py[1]) +
                        Py[2] * (36 * pow(Py[0], 2) + 81 * Py[0] * Py[1])) +
                   Px[0] *
                       ((-(6 * Py[3]) + 18 * Py[2] - 18 * Py[1] + 6 * Py[0]) *
                            pow(y, 2) +
                        (6 * pow(Py[3], 2) + 27 * pow(Py[2], 2) -
                         27 * pow(Py[1], 2) - 6 * pow(Py[0], 2) +
                         (45 * Py[1] - 45 * Py[2]) * Py[3] -
                         45 * Py[0] * Py[2] + 45 * Py[0] * Py[1]) *
                            y -
                        6 * Py[0] * pow(Py[3], 2) - 54 * pow(Py[2], 3) +
                        (81 * Py[1] + 54 * Py[0]) * pow(Py[2], 2) +
                        54 * pow(Py[1], 3) +
                        Py[3] * (-(54 * pow(Py[1], 2)) + 6 * pow(Py[0], 2) +
                                 (54 * Py[1] - 9 * Py[0]) * Py[2] +
                                 9 * Py[0] * Py[1]) +
                        Py[2] * (-(81 * pow(Py[1], 2)) - 54 * Py[0] * Py[1])) +
                   Px[2] *
                       ((-(18 * Py[3]) + 54 * Py[2] - 54 * Py[1] + 18 * Py[0]) *
                            pow(y, 2) +
                        (-(54 * pow(Py[2], 2)) + 81 * pow(Py[1], 2) +
                         63 * pow(Py[0], 2) + (45 * Py[0] - 9 * Py[1]) * Py[3] +
                         (81 * Py[1] - 81 * Py[0]) * Py[2] -
                         126 * Py[0] * Py[1]) *
                            y +
                        54 * Py[0] * pow(Py[2], 2) +
                        Py[2] * (-(81 * pow(Py[1], 2)) - 54 * pow(Py[0], 2) +
                                 81 * Py[0] * Py[1]) +
                        18 * Py[1] * pow(Py[0], 2) +
                        Py[3] * (9 * Py[0] * Py[1] - 27 * pow(Py[0], 2)))) +
              pow(Px[1], 2) *
                  ((27 * Py[3] - 81 * Py[2] + 81 * Py[1] - 27 * Py[0]) *
                       pow(y, 2) +
                   (54 * pow(Py[3], 2) + 81 * pow(Py[2], 2) +
                    (-(81 * Py[2]) - 108 * Py[1] + 27 * Py[0]) * Py[3] +
                    81 * Py[0] * Py[2] - 54 * Py[0] * Py[1]) *
                       y +
                   (27 * Py[1] - 81 * Py[0]) * pow(Py[3], 2) -
                   81 * Py[0] * pow(Py[2], 2) +
                   (81 * Py[0] * Py[2] + 54 * Py[0] * Py[1]) * Py[3]) +
              pow(Px[2], 2) *
                  ((27 * Py[3] - 81 * Py[2] + 81 * Py[1] - 27 * Py[0]) *
                       pow(y, 2) +
                   (-(81 * pow(Py[1], 2)) - 54 * pow(Py[0], 2) +
                    (54 * Py[2] - 81 * Py[1] - 27 * Py[0]) * Py[3] +
                    108 * Py[0] * Py[2] + 81 * Py[0] * Py[1]) *
                       y +
                   Py[3] * (81 * pow(Py[1], 2) + 81 * pow(Py[0], 2) -
                            54 * Py[0] * Py[2] - 81 * Py[0] * Py[1]) -
                   27 * Py[2] * pow(Py[0], 2)) +
              pow(Px[0], 2) *
                  ((3 * Py[3] - 9 * Py[2] + 9 * Py[1] - 3 * Py[0]) * pow(y, 2) +
                   (21 * pow(Py[3], 2) + 54 * pow(Py[2], 2) +
                    (-(63 * Py[2]) + 9 * Py[1] + 6 * Py[0]) * Py[3] -
                    27 * Py[1] * Py[2]) *
                       y +
                   3 * pow(Py[3], 3) +
                   (-(9 * Py[2]) - 18 * Py[1] - 3 * Py[0]) * pow(Py[3], 2) -
                   27 * pow(Py[2], 3) +
                   Py[3] * (27 * pow(Py[2], 2) + 27 * Py[1] * Py[2])) +
              pow(Px[3], 2) *
                  ((3 * Py[3] - 9 * Py[2] + 9 * Py[1] - 3 * Py[0]) * pow(y, 2) +
                   (-(54 * pow(Py[1], 2)) - 21 * pow(Py[0], 2) -
                    6 * Py[0] * Py[3] + (27 * Py[1] - 9 * Py[0]) * Py[2] +
                    63 * Py[0] * Py[1]) *
                       y +
                   27 * pow(Py[1], 3) - 27 * Py[0] * pow(Py[1], 2) -
                   3 * pow(Py[0], 3) +
                   Py[2] * (18 * pow(Py[0], 2) - 27 * Py[0] * Py[1]) +
                   3 * Py[3] * pow(Py[0], 2) + 9 * Py[1] * pow(Py[0], 2)) +
              Px[0] * Px[1] *
                  ((-(18 * Py[3]) + 54 * Py[2] - 54 * Py[1] + 18 * Py[0]) *
                       pow(y, 2) +
                   (-(63 * pow(Py[3], 2)) - 81 * pow(Py[2], 2) +
                    54 * pow(Py[1], 2) +
                    (126 * Py[2] + 81 * Py[1] - 45 * Py[0]) * Py[3] +
                    (9 * Py[0] - 81 * Py[1]) * Py[2]) *
                       y +
                   (-(18 * Py[2]) + 54 * Py[1] + 27 * Py[0]) * pow(Py[3], 2) +
                   81 * Py[1] * pow(Py[2], 2) +
                   Py[3] * ((-(81 * Py[1]) - 9 * Py[0]) * Py[2] -
                            54 * pow(Py[1], 2)))) +
         pow(x, 2) *
             (Px[2] *
                  ((9 * pow(Py[3], 2) + 81 * pow(Py[2], 2) +
                    81 * pow(Py[1], 2) + 9 * pow(Py[0], 2) +
                    (-(54 * Py[2]) + 54 * Py[1] - 18 * Py[0]) * Py[3] +
                    (54 * Py[0] - 162 * Py[1]) * Py[2] - 54 * Py[0] * Py[1]) *
                       y +
                   (18 * Py[1] - 27 * Py[0]) * pow(Py[3], 2) -
                   54 * Py[0] * pow(Py[2], 2) +
                   Py[3] *
                       (-(27 * pow(Py[2], 2)) - 162 * pow(Py[1], 2) -
                        54 * pow(Py[0], 2) + (81 * Py[1] + 27 * Py[0]) * Py[2] +
                        153 * Py[0] * Py[1]) +
                   Py[2] * (81 * pow(Py[1], 2) + 54 * pow(Py[0], 2) -
                            81 * Py[0] * Py[1]) -
                   9 * Py[1] * pow(Py[0], 2)) +
              Px[0] *
                  ((3 * pow(Py[3], 2) + 27 * pow(Py[2], 2) +
                    27 * pow(Py[1], 2) + 3 * pow(Py[0], 2) +
                    (-(18 * Py[2]) + 18 * Py[1] - 6 * Py[0]) * Py[3] +
                    (18 * Py[0] - 54 * Py[1]) * Py[2] - 18 * Py[0] * Py[1]) *
                       y -
                   3 * pow(Py[3], 3) +
                   (18 * Py[2] + 9 * Py[1] - 21 * Py[0]) * pow(Py[3], 2) +
                   54 * pow(Py[2], 3) +
                   (-(81 * Py[1]) - 54 * Py[0]) * pow(Py[2], 2) +
                   Py[3] *
                       (-(54 * pow(Py[2], 2)) - 27 * pow(Py[1], 2) -
                        3 * pow(Py[0], 2) + (27 * Py[1] + 63 * Py[0]) * Py[2] -
                        9 * Py[0] * Py[1]) -
                   27 * pow(Py[1], 3) +
                   Py[2] * (81 * pow(Py[1], 2) + 27 * Py[0] * Py[1])) +
              Px[3] *
                  ((-(3 * pow(Py[3], 2)) - 27 * pow(Py[2], 2) -
                    27 * pow(Py[1], 2) - 3 * pow(Py[0], 2) +
                    (18 * Py[2] - 18 * Py[1] + 6 * Py[0]) * Py[3] +
                    (54 * Py[1] - 18 * Py[0]) * Py[2] + 18 * Py[0] * Py[1]) *
                       y +
                   3 * Py[0] * pow(Py[3], 2) + 27 * pow(Py[2], 3) +
                   (27 * Py[0] - 81 * Py[1]) * pow(Py[2], 2) -
                   54 * pow(Py[1], 3) +
                   Py[2] * (81 * pow(Py[1], 2) - 9 * pow(Py[0], 2) -
                            27 * Py[0] * Py[1]) +
                   Py[3] *
                       (54 * pow(Py[1], 2) + 21 * pow(Py[0], 2) +
                        (9 * Py[0] - 27 * Py[1]) * Py[2] - 63 * Py[0] * Py[1]) +
                   54 * Py[0] * pow(Py[1], 2) + 3 * pow(Py[0], 3) -
                   18 * Py[1] * pow(Py[0], 2)) +
              Px[1] *
                  ((-(9 * pow(Py[3], 2)) - 81 * pow(Py[2], 2) -
                    81 * pow(Py[1], 2) - 9 * pow(Py[0], 2) +
                    (54 * Py[2] - 54 * Py[1] + 18 * Py[0]) * Py[3] +
                    (162 * Py[1] - 54 * Py[0]) * Py[2] + 54 * Py[0] * Py[1]) *
                       y +
                   (9 * Py[2] - 54 * Py[1] + 54 * Py[0]) * pow(Py[3], 2) +
                   (162 * Py[0] - 81 * Py[1]) * pow(Py[2], 2) +
                   Py[3] * (54 * pow(Py[1], 2) + 27 * pow(Py[0], 2) +
                            (81 * Py[1] - 153 * Py[0]) * Py[2] -
                            27 * Py[0] * Py[1]) +
                   27 * Py[0] * pow(Py[1], 2) +
                   Py[2] * (-(18 * pow(Py[0], 2)) - 81 * Py[0] * Py[1]))) +
         (pow(Py[3], 3) +
          (-(9 * Py[2]) + 9 * Py[1] - 3 * Py[0]) * pow(Py[3], 2) -
          27 * pow(Py[2], 3) +
          Py[3] * (27 * pow(Py[2], 2) + 27 * pow(Py[1], 2) + 3 * pow(Py[0], 2) +
                   (18 * Py[0] - 54 * Py[1]) * Py[2] - 18 * Py[0] * Py[1]) +
          (81 * Py[1] - 27 * Py[0]) * pow(Py[2], 2) + 27 * pow(Py[1], 3) -
          27 * Py[0] * pow(Py[1], 2) +
          Py[2] *
              (-(81 * pow(Py[1], 2)) - 9 * pow(Py[0], 2) + 54 * Py[0] * Py[1]) -
          pow(Py[0], 3) + 9 * Py[1] * pow(Py[0], 2)) *
             pow(x, 3);
}

double cubic_x(double x, double y) {
  return Px[2] * (Px[0] * ((18 * Py[3] - 54 * Py[2] + 54 * Py[1] - 18 * Py[0]) *
                               pow(y, 2) +
                           (9 * pow(Py[3], 2) - 108 * pow(Py[2], 2) -
                            81 * pow(Py[1], 2) +
                            (81 * Py[2] - 180 * Py[1] + 45 * Py[0]) * Py[3] +
                            243 * Py[1] * Py[2] - 9 * Py[0] * Py[1]) *
                               y +
                           (27 * Py[0] - 36 * Py[1]) * pow(Py[3], 2) +
                           Py[3] * (54 * pow(Py[2], 2) + 162 * pow(Py[1], 2) +
                                    (-(81 * Py[1]) - 108 * Py[0]) * Py[2] +
                                    9 * Py[0] * Py[1]) +
                           54 * Py[0] * pow(Py[2], 2) -
                           81 * Py[2] * pow(Py[1], 2)) +
                  Px[1] * ((-(54 * Py[3]) + 162 * Py[2] - 162 * Py[1] +
                            54 * Py[0]) *
                               pow(y, 2) +
                           (-(27 * pow(Py[3], 2)) + 27 * pow(Py[0], 2) +
                            (243 * Py[1] - 81 * Py[2]) * Py[3] -
                            243 * Py[0] * Py[2] + 81 * Py[0] * Py[1]) *
                               y +
                           27 * Py[0] * pow(Py[3], 2) +
                           Py[3] * (-(27 * pow(Py[0], 2)) +
                                    (162 * Py[0] - 81 * Py[1]) * Py[2] -
                                    162 * Py[0] * Py[1]) +
                           81 * Py[0] * Py[1] * Py[2])) +
         Px[3] *
             (Px[1] *
                  ((18 * Py[3] - 54 * Py[2] + 54 * Py[1] - 18 * Py[0]) *
                       pow(y, 2) +
                   (81 * pow(Py[2], 2) + 108 * pow(Py[1], 2) -
                    9 * pow(Py[0], 2) + (9 * Py[2] - 45 * Py[0]) * Py[3] +
                    (180 * Py[0] - 243 * Py[1]) * Py[2] - 81 * Py[0] * Py[1]) *
                       y +
                   (81 * Py[1] - 162 * Py[0]) * pow(Py[2], 2) -
                   54 * Py[0] * pow(Py[1], 2) +
                   Py[3] * (-(54 * pow(Py[1], 2)) - 27 * pow(Py[0], 2) -
                            9 * Py[0] * Py[2] + 108 * Py[0] * Py[1]) +
                   Py[2] * (36 * pow(Py[0], 2) + 81 * Py[0] * Py[1])) +
              Px[0] * ((-(6 * Py[3]) + 18 * Py[2] - 18 * Py[1] + 6 * Py[0]) *
                           pow(y, 2) +
                       (6 * pow(Py[3], 2) + 27 * pow(Py[2], 2) -
                        27 * pow(Py[1], 2) - 6 * pow(Py[0], 2) +
                        (45 * Py[1] - 45 * Py[2]) * Py[3] - 45 * Py[0] * Py[2] +
                        45 * Py[0] * Py[1]) *
                           y -
                       6 * Py[0] * pow(Py[3], 2) - 54 * pow(Py[2], 3) +
                       (81 * Py[1] + 54 * Py[0]) * pow(Py[2], 2) +
                       54 * pow(Py[1], 3) +
                       Py[3] * (-(54 * pow(Py[1], 2)) + 6 * pow(Py[0], 2) +
                                (54 * Py[1] - 9 * Py[0]) * Py[2] +
                                9 * Py[0] * Py[1]) +
                       Py[2] * (-(81 * pow(Py[1], 2)) - 54 * Py[0] * Py[1])) +
              Px[2] *
                  ((-(18 * Py[3]) + 54 * Py[2] - 54 * Py[1] + 18 * Py[0]) *
                       pow(y, 2) +
                   (-(54 * pow(Py[2], 2)) + 81 * pow(Py[1], 2) +
                    63 * pow(Py[0], 2) + (45 * Py[0] - 9 * Py[1]) * Py[3] +
                    (81 * Py[1] - 81 * Py[0]) * Py[2] - 126 * Py[0] * Py[1]) *
                       y +
                   54 * Py[0] * pow(Py[2], 2) +
                   Py[2] * (-(81 * pow(Py[1], 2)) - 54 * pow(Py[0], 2) +
                            81 * Py[0] * Py[1]) +
                   18 * Py[1] * pow(Py[0], 2) +
                   Py[3] * (9 * Py[0] * Py[1] - 27 * pow(Py[0], 2)))) +
         pow(Px[1], 2) *
             ((27 * Py[3] - 81 * Py[2] + 81 * Py[1] - 27 * Py[0]) * pow(y, 2) +
              (54 * pow(Py[3], 2) + 81 * pow(Py[2], 2) +
               (-(81 * Py[2]) - 108 * Py[1] + 27 * Py[0]) * Py[3] +
               81 * Py[0] * Py[2] - 54 * Py[0] * Py[1]) *
                  y +
              (27 * Py[1] - 81 * Py[0]) * pow(Py[3], 2) -
              81 * Py[0] * pow(Py[2], 2) +
              (81 * Py[0] * Py[2] + 54 * Py[0] * Py[1]) * Py[3]) +
         pow(Px[2], 2) *
             ((27 * Py[3] - 81 * Py[2] + 81 * Py[1] - 27 * Py[0]) * pow(y, 2) +
              (-(81 * pow(Py[1], 2)) - 54 * pow(Py[0], 2) +
               (54 * Py[2] - 81 * Py[1] - 27 * Py[0]) * Py[3] +
               108 * Py[0] * Py[2] + 81 * Py[0] * Py[1]) *
                  y +
              Py[3] * (81 * pow(Py[1], 2) + 81 * pow(Py[0], 2) -
                       54 * Py[0] * Py[2] - 81 * Py[0] * Py[1]) -
              27 * Py[2] * pow(Py[0], 2)) +
         pow(Px[0], 2) *
             ((3 * Py[3] - 9 * Py[2] + 9 * Py[1] - 3 * Py[0]) * pow(y, 2) +
              (21 * pow(Py[3], 2) + 54 * pow(Py[2], 2) +
               (-(63 * Py[2]) + 9 * Py[1] + 6 * Py[0]) * Py[3] -
               27 * Py[1] * Py[2]) *
                  y +
              3 * pow(Py[3], 3) +
              (-(9 * Py[2]) - 18 * Py[1] - 3 * Py[0]) * pow(Py[3], 2) -
              27 * pow(Py[2], 3) +
              Py[3] * (27 * pow(Py[2], 2) + 27 * Py[1] * Py[2])) +
         pow(Px[3], 2) *
             ((3 * Py[3] - 9 * Py[2] + 9 * Py[1] - 3 * Py[0]) * pow(y, 2) +
              (-(54 * pow(Py[1], 2)) - 21 * pow(Py[0], 2) - 6 * Py[0] * Py[3] +
               (27 * Py[1] - 9 * Py[0]) * Py[2] + 63 * Py[0] * Py[1]) *
                  y +
              27 * pow(Py[1], 3) - 27 * Py[0] * pow(Py[1], 2) -
              3 * pow(Py[0], 3) +
              Py[2] * (18 * pow(Py[0], 2) - 27 * Py[0] * Py[1]) +
              3 * Py[3] * pow(Py[0], 2) + 9 * Py[1] * pow(Py[0], 2)) +
         Px[0] * Px[1] *
             ((-(18 * Py[3]) + 54 * Py[2] - 54 * Py[1] + 18 * Py[0]) *
                  pow(y, 2) +
              (-(63 * pow(Py[3], 2)) - 81 * pow(Py[2], 2) + 54 * pow(Py[1], 2) +
               (126 * Py[2] + 81 * Py[1] - 45 * Py[0]) * Py[3] +
               (9 * Py[0] - 81 * Py[1]) * Py[2]) *
                  y +
              (-(18 * Py[2]) + 54 * Py[1] + 27 * Py[0]) * pow(Py[3], 2) +
              81 * Py[1] * pow(Py[2], 2) +
              Py[3] *
                  ((-(81 * Py[1]) - 9 * Py[0]) * Py[2] - 54 * pow(Py[1], 2))) +
         2 * x *
             (Px[2] *
                  ((9 * pow(Py[3], 2) + 81 * pow(Py[2], 2) +
                    81 * pow(Py[1], 2) + 9 * pow(Py[0], 2) +
                    (-(54 * Py[2]) + 54 * Py[1] - 18 * Py[0]) * Py[3] +
                    (54 * Py[0] - 162 * Py[1]) * Py[2] - 54 * Py[0] * Py[1]) *
                       y +
                   (18 * Py[1] - 27 * Py[0]) * pow(Py[3], 2) -
                   54 * Py[0] * pow(Py[2], 2) +
                   Py[3] *
                       (-(27 * pow(Py[2], 2)) - 162 * pow(Py[1], 2) -
                        54 * pow(Py[0], 2) + (81 * Py[1] + 27 * Py[0]) * Py[2] +
                        153 * Py[0] * Py[1]) +
                   Py[2] * (81 * pow(Py[1], 2) + 54 * pow(Py[0], 2) -
                            81 * Py[0] * Py[1]) -
                   9 * Py[1] * pow(Py[0], 2)) +
              Px[0] *
                  ((3 * pow(Py[3], 2) + 27 * pow(Py[2], 2) +
                    27 * pow(Py[1], 2) + 3 * pow(Py[0], 2) +
                    (-(18 * Py[2]) + 18 * Py[1] - 6 * Py[0]) * Py[3] +
                    (18 * Py[0] - 54 * Py[1]) * Py[2] - 18 * Py[0] * Py[1]) *
                       y -
                   3 * pow(Py[3], 3) +
                   (18 * Py[2] + 9 * Py[1] - 21 * Py[0]) * pow(Py[3], 2) +
                   54 * pow(Py[2], 3) +
                   (-(81 * Py[1]) - 54 * Py[0]) * pow(Py[2], 2) +
                   Py[3] *
                       (-(54 * pow(Py[2], 2)) - 27 * pow(Py[1], 2) -
                        3 * pow(Py[0], 2) + (27 * Py[1] + 63 * Py[0]) * Py[2] -
                        9 * Py[0] * Py[1]) -
                   27 * pow(Py[1], 3) +
                   Py[2] * (81 * pow(Py[1], 2) + 27 * Py[0] * Py[1])) +
              Px[3] *
                  ((-(3 * pow(Py[3], 2)) - 27 * pow(Py[2], 2) -
                    27 * pow(Py[1], 2) - 3 * pow(Py[0], 2) +
                    (18 * Py[2] - 18 * Py[1] + 6 * Py[0]) * Py[3] +
                    (54 * Py[1] - 18 * Py[0]) * Py[2] + 18 * Py[0] * Py[1]) *
                       y +
                   3 * Py[0] * pow(Py[3], 2) + 27 * pow(Py[2], 3) +
                   (27 * Py[0] - 81 * Py[1]) * pow(Py[2], 2) -
                   54 * pow(Py[1], 3) +
                   Py[2] * (81 * pow(Py[1], 2) - 9 * pow(Py[0], 2) -
                            27 * Py[0] * Py[1]) +
                   Py[3] *
                       (54 * pow(Py[1], 2) + 21 * pow(Py[0], 2) +
                        (9 * Py[0] - 27 * Py[1]) * Py[2] - 63 * Py[0] * Py[1]) +
                   54 * Py[0] * pow(Py[1], 2) + 3 * pow(Py[0], 3) -
                   18 * Py[1] * pow(Py[0], 2)) +
              Px[1] *
                  ((-(9 * pow(Py[3], 2)) - 81 * pow(Py[2], 2) -
                    81 * pow(Py[1], 2) - 9 * pow(Py[0], 2) +
                    (54 * Py[2] - 54 * Py[1] + 18 * Py[0]) * Py[3] +
                    (162 * Py[1] - 54 * Py[0]) * Py[2] + 54 * Py[0] * Py[1]) *
                       y +
                   (9 * Py[2] - 54 * Py[1] + 54 * Py[0]) * pow(Py[3], 2) +
                   (162 * Py[0] - 81 * Py[1]) * pow(Py[2], 2) +
                   Py[3] * (54 * pow(Py[1], 2) + 27 * pow(Py[0], 2) +
                            (81 * Py[1] - 153 * Py[0]) * Py[2] -
                            27 * Py[0] * Py[1]) +
                   27 * Py[0] * pow(Py[1], 2) +
                   Py[2] * (-(18 * pow(Py[0], 2)) - 81 * Py[0] * Py[1]))) +
         3 * (pow(Py[3], 3) + (-(9 * Py[2]) + 9 * Py[1] - 3 * Py[0]) * pow(Py[3], 2) - 27 * pow(Py[2], 3) + Py[3] * (27 * pow(Py[2], 2) + 27 * pow(Py[1], 2) + 3 * pow(Py[0], 2) + (18 * Py[0] - 54 * Py[1]) * Py[2] - 18 * Py[0] * Py[1]) + (81 * Py[1] - 27 * Py[0]) * pow(Py[2], 2) + 27 * pow(Py[1], 3) - 27 * Py[0] * pow(Py[1], 2) + Py[2] * (-(81 * pow(Py[1], 2)) - 9 * pow(Py[0], 2) + 54 * Py[0] * Py[1]) - pow(Py[0], 3) + 9 * Py[1] * pow(Py[0], 2)) *
             pow(x, 2);
}

double cubic_y(double x, double y) {
  return Px[3] *
             (Px[2] * (Px[1] * (162 * pow(y, 2) +
                                2 *
                                    (27 * Py[3] - 81 * Py[2] - 81 * Py[1] -
                                     27 * Py[0]) *
                                    y -
                                27 * pow(Py[0], 2) - 54 * Py[0] * Py[3] +
                                (81 * Py[1] + 81 * Py[0]) * Py[2] +
                                81 * Py[0] * Py[1]) +
                       Px[0] * (-(54 * pow(y, 2)) +
                                2 *
                                    (-(9 * Py[3]) - 27 * Py[2] + 153 * Py[1] -
                                     63 * Py[0]) *
                                    y +
                                54 * pow(Py[2], 2) - 81 * pow(Py[1], 2) +
                                (9 * Py[1] + 9 * Py[0]) * Py[3] +
                                (108 * Py[0] - 162 * Py[1]) * Py[2] +
                                9 * Py[0] * Py[1])) +
              Px[0] * Px[1] *
                  (54 * pow(y, 2) +
                   2 * (63 * Py[3] - 153 * Py[2] + 27 * Py[1] + 9 * Py[0]) * y +
                   81 * pow(Py[2], 2) - 54 * pow(Py[1], 2) +
                   (-(9 * Py[2]) - 108 * Py[1] - 9 * Py[0]) * Py[3] +
                   (162 * Py[1] - 9 * Py[0]) * Py[2]) +
              pow(Px[0], 2) *
                  (-(9 * pow(y, 2)) +
                   2 * (-(21 * Py[3]) + 54 * Py[2] - 27 * Py[1] + 3 * Py[0]) *
                       y -
                   3 * pow(Py[3], 2) - 81 * pow(Py[2], 2) +
                   (27 * Py[2] + 27 * Py[1] - 6 * Py[0]) * Py[3] +
                   27 * Py[1] * Py[2]) +
              pow(Px[1], 2) *
                  (-(81 * pow(y, 2)) +
                   2 * (-(54 * Py[3]) + 162 * Py[2] - 54 * Py[1] + 27 * Py[0]) *
                       y -
                   81 * pow(Py[2], 2) + (54 * Py[1] + 54 * Py[0]) * Py[3] -
                   162 * Py[0] * Py[2] + 54 * Py[0] * Py[1]) +
              pow(Px[2], 2) *
                  (-(81 * pow(y, 2)) + 2 * (27 * Py[2] + 54 * Py[0]) * y -
                   27 * pow(Py[0], 2) - 54 * Py[0] * Py[2])) +
         Px[2] *
             (pow(Px[1], 2) *
                  (243 * pow(y, 2) +
                   2 * (-(81 * Py[3]) - 81 * Py[2] - 81 * Py[0]) * y +
                   (81 * Py[2] + 81 * Py[0]) * Py[3] + 81 * Py[0] * Py[2]) +
              pow(Px[0], 2) * (27 * pow(y, 2) +
                               2 * (9 * Py[3] - 54 * Py[2] + 18 * Py[1]) * y -
                               18 * pow(Py[3], 2) + 27 * pow(Py[2], 2) +
                               (54 * Py[2] - 36 * Py[1]) * Py[3]) +
              Px[0] * Px[1] *
                  (-(162 * pow(y, 2)) +
                   2 * (27 * Py[3] + 81 * Py[2] + 81 * Py[1] - 27 * Py[0]) * y +
                   27 * pow(Py[3], 2) +
                   (-(81 * Py[2]) - 81 * Py[1] + 54 * Py[0]) * Py[3] -
                   81 * Py[1] * Py[2])) +
         pow(Px[2], 2) *
             (Px[0] *
                  (81 * pow(y, 2) +
                   2 * (-(27 * Py[3]) + 54 * Py[2] - 162 * Py[1] + 54 * Py[0]) *
                       y +
                   81 * pow(Py[1], 2) +
                   (-(54 * Py[2]) + 162 * Py[1] - 54 * Py[0]) * Py[3] -
                   54 * Py[0] * Py[2]) +
              Px[1] *
                  (-(243 * pow(y, 2)) +
                   2 * (81 * Py[3] + 81 * Py[1] + 81 * Py[0]) * y +
                   (-(81 * Py[1]) - 81 * Py[0]) * Py[3] - 81 * Py[0] * Py[1])) +
         pow(Px[3], 2) *
             (Px[2] * (27 * pow(y, 2) + 2 * (-(9 * Py[1]) - 18 * Py[0]) * y +
                       9 * pow(Py[0], 2) + 18 * Py[0] * Py[1]) +
              Px[0] *
                  (9 * pow(y, 2) +
                   2 * (-(3 * Py[3]) + 27 * Py[2] - 54 * Py[1] + 21 * Py[0]) *
                       y +
                   81 * pow(Py[1], 2) + 3 * pow(Py[0], 2) + 6 * Py[0] * Py[3] +
                   (-(27 * Py[1]) - 27 * Py[0]) * Py[2] - 27 * Py[0] * Py[1]) +
              Px[1] * (-(27 * pow(y, 2)) +
                       2 * (-(18 * Py[2]) + 54 * Py[1] - 9 * Py[0]) * y -
                       27 * pow(Py[1], 2) + 18 * pow(Py[0], 2) +
                       36 * Py[0] * Py[2] - 54 * Py[0] * Py[1])) +
         pow(Px[2], 3) *
             (81 * pow(y, 2) + 2 * (-(27 * Py[3]) - 54 * Py[0]) * y +
              27 * pow(Py[0], 2) + 54 * Py[0] * Py[3]) +
         Px[0] * pow(Px[1], 2) *
             (81 * pow(y, 2) + 2 * (-(54 * Py[3]) - 27 * Py[1]) * y +
              27 * pow(Py[3], 2) + 54 * Py[1] * Py[3]) +
         pow(Px[0], 3) * (3 * pow(y, 2) - 6 * Py[3] * y + 3 * pow(Py[3], 2)) +
         pow(Px[3], 3) *
             (-(3 * pow(y, 2)) + 6 * Py[0] * y - 3 * pow(Py[0], 2)) +
         Px[1] * pow(Px[0], 2) *
             (-(27 * pow(y, 2)) + 2 * (18 * Py[3] + 9 * Py[2]) * y -
              9 * pow(Py[3], 2) - 18 * Py[2] * Py[3]) +
         pow(Px[1], 3) *
             (-(81 * pow(y, 2)) + 2 * (54 * Py[3] + 27 * Py[0]) * y -
              27 * pow(Py[3], 2) - 54 * Py[0] * Py[3]) +
         x * (Px[2] *
                  (Px[0] * (2 *
                                (18 * Py[3] - 54 * Py[2] + 54 * Py[1] -
                                 18 * Py[0]) *
                                y +
                            9 * pow(Py[3], 2) - 108 * pow(Py[2], 2) -
                            81 * pow(Py[1], 2) +
                            (81 * Py[2] - 180 * Py[1] + 45 * Py[0]) * Py[3] +
                            243 * Py[1] * Py[2] - 9 * Py[0] * Py[1]) +
                   Px[1] * (2 *
                                (-(54 * Py[3]) + 162 * Py[2] - 162 * Py[1] +
                                 54 * Py[0]) *
                                y -
                            27 * pow(Py[3], 2) + 27 * pow(Py[0], 2) +
                            (243 * Py[1] - 81 * Py[2]) * Py[3] -
                            243 * Py[0] * Py[2] + 81 * Py[0] * Py[1])) +
              Px[3] * (Px[1] * (2 *
                                    (18 * Py[3] - 54 * Py[2] + 54 * Py[1] -
                                     18 * Py[0]) *
                                    y +
                                81 * pow(Py[2], 2) + 108 * pow(Py[1], 2) -
                                9 * pow(Py[0], 2) +
                                (9 * Py[2] - 45 * Py[0]) * Py[3] +
                                (180 * Py[0] - 243 * Py[1]) * Py[2] -
                                81 * Py[0] * Py[1]) +
                       Px[0] * (2 *
                                    (-(6 * Py[3]) + 18 * Py[2] - 18 * Py[1] +
                                     6 * Py[0]) *
                                    y +
                                6 * pow(Py[3], 2) + 27 * pow(Py[2], 2) -
                                27 * pow(Py[1], 2) - 6 * pow(Py[0], 2) +
                                (45 * Py[1] - 45 * Py[2]) * Py[3] -
                                45 * Py[0] * Py[2] + 45 * Py[0] * Py[1]) +
                       Px[2] * (2 *
                                    (-(18 * Py[3]) + 54 * Py[2] - 54 * Py[1] +
                                     18 * Py[0]) *
                                    y -
                                54 * pow(Py[2], 2) + 81 * pow(Py[1], 2) +
                                63 * pow(Py[0], 2) +
                                (45 * Py[0] - 9 * Py[1]) * Py[3] +
                                (81 * Py[1] - 81 * Py[0]) * Py[2] -
                                126 * Py[0] * Py[1])) +
              pow(Px[1], 2) *
                  (2 * (27 * Py[3] - 81 * Py[2] + 81 * Py[1] - 27 * Py[0]) * y +
                   54 * pow(Py[3], 2) + 81 * pow(Py[2], 2) +
                   (-(81 * Py[2]) - 108 * Py[1] + 27 * Py[0]) * Py[3] +
                   81 * Py[0] * Py[2] - 54 * Py[0] * Py[1]) +
              pow(Px[2], 2) *
                  (2 * (27 * Py[3] - 81 * Py[2] + 81 * Py[1] - 27 * Py[0]) * y -
                   81 * pow(Py[1], 2) - 54 * pow(Py[0], 2) +
                   (54 * Py[2] - 81 * Py[1] - 27 * Py[0]) * Py[3] +
                   108 * Py[0] * Py[2] + 81 * Py[0] * Py[1]) +
              pow(Px[0], 2) *
                  (2 * (3 * Py[3] - 9 * Py[2] + 9 * Py[1] - 3 * Py[0]) * y +
                   21 * pow(Py[3], 2) + 54 * pow(Py[2], 2) +
                   (-(63 * Py[2]) + 9 * Py[1] + 6 * Py[0]) * Py[3] -
                   27 * Py[1] * Py[2]) +
              pow(Px[3], 2) *
                  (2 * (3 * Py[3] - 9 * Py[2] + 9 * Py[1] - 3 * Py[0]) * y -
                   54 * pow(Py[1], 2) - 21 * pow(Py[0], 2) - 6 * Py[0] * Py[3] +
                   (27 * Py[1] - 9 * Py[0]) * Py[2] + 63 * Py[0] * Py[1]) +
              Px[0] * Px[1] *
                  (2 * (-(18 * Py[3]) + 54 * Py[2] - 54 * Py[1] + 18 * Py[0]) *
                       y -
                   63 * pow(Py[3], 2) - 81 * pow(Py[2], 2) +
                   54 * pow(Py[1], 2) +
                   (126 * Py[2] + 81 * Py[1] - 45 * Py[0]) * Py[3] +
                   (9 * Py[0] - 81 * Py[1]) * Py[2])) +
         (Px[2] * (9 * pow(Py[3], 2) + 81 * pow(Py[2], 2) + 81 * pow(Py[1], 2) +
                   9 * pow(Py[0], 2) +
                   (-(54 * Py[2]) + 54 * Py[1] - 18 * Py[0]) * Py[3] +
                   (54 * Py[0] - 162 * Py[1]) * Py[2] - 54 * Py[0] * Py[1]) +
          Px[0] * (3 * pow(Py[3], 2) + 27 * pow(Py[2], 2) + 27 * pow(Py[1], 2) +
                   3 * pow(Py[0], 2) +
                   (-(18 * Py[2]) + 18 * Py[1] - 6 * Py[0]) * Py[3] +
                   (18 * Py[0] - 54 * Py[1]) * Py[2] - 18 * Py[0] * Py[1]) +
          Px[3] * (-(3 * pow(Py[3], 2)) - 27 * pow(Py[2], 2) -
                   27 * pow(Py[1], 2) - 3 * pow(Py[0], 2) +
                   (18 * Py[2] - 18 * Py[1] + 6 * Py[0]) * Py[3] +
                   (54 * Py[1] - 18 * Py[0]) * Py[2] + 18 * Py[0] * Py[1]) +
          Px[1] * (-(9 * pow(Py[3], 2)) - 81 * pow(Py[2], 2) -
                   81 * pow(Py[1], 2) - 9 * pow(Py[0], 2) +
                   (54 * Py[2] - 54 * Py[1] + 18 * Py[0]) * Py[3] +
                   (162 * Py[1] - 54 * Py[0]) * Py[2] + 54 * Py[0] * Py[1])) *
             pow(x, 2);
}
