#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>

void Sum(double x, double &s1) {
    double eps = 0.000001;
    double s = x;
    double dr = x;
    int i = 1;
    s1 = x;

    do {
        s = s1;
        dr = -1 * dr * x * x / (2 * i * (2 * i + 1));
        s1 = s + dr / (2 * i + 1);
        i++;
    } while (std::abs(s1 - s) > eps);
}

void Interpolation(double x1, std::vector<double> &x, std::vector<double> &si, int n, double &l) {
    std::vector<double> x_(n);
    std::vector<double> si_(n);

    for (int i = 0; i < n; i++) {
        x_[i] = x[i * 2];
        si_[i] = si[i * 2];
    }

    l = 0;

    for (int i = 0; i < n; i++) {
        double p = 1;

        for (int j = 0; j < n; j++) {
            if (j != i) {
                p = p * ((x1 - x_[j]) / (x_[i] - x_[j]));
            }
        }

        l += si_[i] * p;
    }
}

void Error(std::vector<double> &a, std::vector<double> &b, int n, std::vector<double> &e) {
    for (int i = 0; i < n; i++) {
        e[i] = std::abs(a[i] - b[i]);
    }
}

void CalculateInterpolationAndError(int m, std::ofstream &file, std::vector<double> &max_errors, std::vector<double> &max_errors_cb) {
    int a = 0;
    int b = 4;
    int n = 2 * m + 1;
    std::vector<double> x(n);
    std::vector<double> si(n);
    std::vector<double> pol(n);
    x[0] = 0;
    si[0] = 0;
    pol[0] = 0;
    std::vector<double> err(n);

    std::vector<double> x_cb(n);
    std::vector<double> si_cb(n);
    std::vector<double> pol_cb(n);
    x_cb[0] = 0;
    si_cb[0] = 0;

    for (int i = 1; i < n; i++) {
        double s, s_cb;
        x[i] = i * (b - a) / static_cast<double>(2 * m);
        x_cb[i] = (b - a) / 2 * std::cos(std::acos(-1.0) * (2 * i + 1) / (2 * n + 2)) + (a + b) / 2;
        Sum(x[i], s);
        Sum(x_cb[i], s_cb);
        file << "Si(" << x[i] << ") = " << s << std::endl;
        file << "Si_cb(" << x_cb[i] << ") = " << s_cb << std::endl;
        si[i] = s;
        si_cb[i] = s_cb;
    }

    file << std::endl;

    for (int i = 0; i < n; i++) {
        double l, l_cb;
        Interpolation(x[i], x, si, m + 1, l);
        Interpolation(x_cb[i], x_cb, si_cb, m + 1, l_cb);
        file << "L(" << x[i] << ") = " << l << std::endl;
        file << "L_cb(" << x_cb[i] << ") = " << l_cb << std::endl;
        pol[i] = l;
        pol_cb[i] = l_cb;
    }

    file << std::endl;

    Error(si, pol, n, err);
    double max = err[0];

    for (int i = 0; i < n; i++) {
        file << "S(" << x[i] << ") - L(" << x[i] << ") = " << err[i] << std::endl;

        if (err[i] > max) {
            max = err[i];
        }
    }

    file << std::fixed << std::setprecision(10) << max << "\t" << std::endl;

    Error(si_cb, pol_cb, n, err);
    double max_cb = err[0];

    for (int i = 0; i < n; i++) {
        file << "S_cb(" << x_cb[i] << ") - L_cb(" << x_cb[i] << ") = " << err[i] << std::endl;

        if (err[i] > max_cb) {
            max_cb = err[i];
        }
    }

    file << std::fixed << std::setprecision(10) << max_cb << "\t" << std::endl;

    max_errors.push_back(max);
    max_errors_cb.push_back(max_cb);

    if (m % 5 == 0) {
        file << m << "\t" << max << "\t" << n << "\t" << max_cb << std::endl;
    }
}

void OutputResults(std::vector<double> &max_errors, std::vector<double> &max_errors_cb) {
    std::cout << "Max Errors:" << std::endl;
    for (size_t i = 0; i < max_errors.size(); i++) {
        std::cout << "Regular " << i + 6 << ": " << std::fixed << std::setprecision(10) << max_errors[i]
                  << "\tChebyshev " << i + 6 << ": " << std::fixed << std::setprecision(10) << max_errors_cb[i] << std::endl;
    }
}

int main() {
    std::ofstream f("Points.txt", std::ios::out);

    std::vector<double> max_errors;
    std::vector<double> max_errors_cb;

    for (int m = 5; m < 50; m++) {
        CalculateInterpolationAndError(m, f, max_errors, max_errors_cb);
    }

    OutputResults(max_errors, max_errors_cb);

    f.close();

    return 0;
}
