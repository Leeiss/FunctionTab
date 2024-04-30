#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

namespace constans {
const int STEPS = 1024;
const float LEFT_BORDER = 0;
const float EPSILON = 1e-6;
const float E = 2.71828182846;
}  // namespace constans

float Tfunc(float x) {
  int n = 0;
  float node_0 = x;
  float ans = x;
  while (fabs(node_0) > 1e-6) {
    float q = ((-1) * x * x * (2 * n + 1)) / ((2 * n + 2) * (2 * n + 3) * (2 * n + 3));
    node_0 *= q;
    ans += node_0;
    n++;
  }
  return ans;
}

float func(float t) {
    if (t != 0) {
        return sin(t) / t; 
    }
    return 1;
}
     
float leftRectangles(float (*func)(float), const float &a, float b,
                      int steps) {
  float result = 0.0;
  float stepSize = (b - a) / steps;
  for (int i = 0; i < steps; i++) {
    float x_i = a + stepSize * i;
    result += stepSize * func(x_i);
  }
  return result;
}

float rightRectangles(float (*func)(float), const float &a, float b,
                       int steps) {
  float result = 0.0;
  float stepSize = (b - a) / steps;
  for (int i = 1; i <= steps; i++) {
    float x_i_1 = a + stepSize * i;
    result += stepSize * func(x_i_1);
  }
  return result;
}

float middleRectangles(float (*func)(float), const float &a, float b,
                        int steps) {
  float result = 0.0;
  float stepSize = (b - a) / steps;
  for (int i = 0; i < steps; i++) {
    float x_i = a + stepSize * i;
    float x_i_1 = a + stepSize * (i + 1);
    result += stepSize * func((x_i + x_i_1) / 2);
  }
  return result;
}

float trapezeFormula(float (*func)(float), const float &a, float b,
                      int steps) {
  float result = func(a) + func(b);
  float stepSize = (b - a) / steps;
  for (int i = 1; i < steps; i++) {
    float x_i_1 = a + stepSize * i;
    result += 2 * func(x_i_1);
  }
  result *= stepSize / 2;
  return result;
}

float SypmsonsFormula(float (*func)(float), const float &a, float b,
                       int steps) {
  float stepSize = (b - a) / steps;
  float result = 0;
  float x = a;
  for (int i = 0; i < steps; i++)
  {
      result += (func(x) + 4 * func(x + stepSize / 2) + func(x + stepSize)) * stepSize / 6;
      x += stepSize;
  }
  return result;
}

float GaussFormula(float (*func)(float), const float &a, float b, int steps) {
  float stepSize = (b - a) / steps;
  float ad1 = (1 - 1.0 / sqrt(3)) * stepSize / 2;
  float ad2 = (1 + 1.0 / sqrt(3)) * stepSize / 2;
  float result = 0;
  float x = a;
  for (int i = 0; i < steps; i++)
  {
      result += (func(x + ad1) + func(x + ad2)) * stepSize / 2;
      x += stepSize;
  }
  return result;
}

void CalculateFunc(vector<float> points,
                   float (*function)(float (*)(float), const float &,
                                      float, int)) {
  for (auto point : points) {
    int steps = 1;
    float lastResult = 0.0;
    float currentResult = 0.0;
    do {
      steps *= 2;
      lastResult = currentResult;
      currentResult = function(func, constans::LEFT_BORDER, point, steps);
    } while (abs(lastResult - currentResult) > constans::EPSILON && steps < constans::STEPS);

    float difference = abs(Tfunc(point) - currentResult);

    printf(
        "x_i = %.1lf | J_o = %.10lf | J_n = %.10lf | |J_o - J_n| = %.10lf | N "
        "= %d\n",
        point, Tfunc(point), currentResult, difference, steps);
  }
}

int main() {
  vector<float> points = {0.0, 0.4, 0.8, 1.2, 1.6,
                           2.0, 2.4, 2.8, 3.2, 3.6, 4.0};
  cout << "Правые прямоугольники\n";
  CalculateFunc(points, rightRectangles);
  cout << "Левые прямоугольники\n";
  CalculateFunc(points, leftRectangles);
  cout << "Центральные прямоугольники\n";
  CalculateFunc(points, middleRectangles);
  cout << "Трапеции\n";
  CalculateFunc(points, trapezeFormula);
  cout << "Симпсон\n";
  CalculateFunc(points, SypmsonsFormula);
  cout << "Гаус\n";
  CalculateFunc(points, GaussFormula);
  return 0;
}
