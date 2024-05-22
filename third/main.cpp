
using namespace std;

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

double f_x(double x) { 
    return -pow(x, 6) + 26 * pow(x, 4) + 4 * pow(x, 3) - 12 * pow(x, 2); 
} 
 
double u_x(double x) { 
    return pow(x, 4) * (1 - x); 
} 
 
double p_x(double x) { 
   return 1 + x;
}

double q_x(double x) {
    return 1 + x;
}

vector<double> func(int n) {
    double h = 1. / n;
    vector<double> b(n,0);
    for (int i = 1; i <= n; ++i) {
        b[i - 1] = f_x(i * h) * pow(h,2);
    }
    return b;
}

double calculate_error(vector<double> &x, vector<vector<double>> &A, vector<double> &b) {

    vector<double> r = calculate_r(A, b, x);
    double max_err = 0.0;
    for (int i = 0; i < r.size(); ++i) {
        if (abs(r[i]) > max_err)
            max_err = abs(r[i]);
    }

    return max_err;
}

void progonka_method(vector<vector<double>> &A, vector<double> &x, int n, vector<double> &b) {
    
    vector<double> alpha(n + 1), betta(n + 1);
    double h = 1.0 / n;

    // прямой ход
    alpha[0] = A[0][1] / A[0][0];
    betta[0] = (b[0]) / A[0][0];

    for (int i = 1; i < n; ++i) {
        double del = 1.0 / (A[i][i] - alpha[i - 1] * A[i][i - 1]);
        alpha[i] = A[i][i + 1] * del;
        betta[i] = (-A[i][i - 1] * betta[i - 1] + b[i]) * del;
    }

    // обратный ход
    x[n - 1] = betta[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = -alpha[i] * x[i + 1] + betta[i];
    }
    for (int i = 1; i <= n; ++i) {
        printf("ih = %4.2lf | y_i = %9.6lf | u(ih) = %8.6lf | |y_i - u(ih)| = %8.6lf\n", i * h, x[i - 1], u_x(i * h), abs(x[i - 1] - u_x(i * h)));
    }
}

int Jacobi_method(int n, vector<double> &x, vector<vector<double>> &A, vector<double> &b) {
    double h = 1.0 / n;
    int k = 0;
    vector<double> new_x(n, 0.);
    double error = 1.0;

    while (error > 1.0 / pow(n, 3)) {
        vector<double> curr_x = x;
        for (int i = 0; i < n; ++i) {
            new_x[i] = Jacobi_calculate_new_x(i, new_x, n, A, b);
        }
        x = new_x;
        error = calculate_error(new_x, A, b);
        k++;
    }
    return k;
}

double Jacobi_calculate_new_x(int i, vector<double>& x, int n, vector<vector<double>> &A, vector<double> &b) {
    double h = 1.0 / n;
    double sum = 0.0;


    if (i > 0) {
        for (int j = 0; j <= i - 1; ++j) {
            sum += A[i][j] * x[j];
        }
    }
    for (int j = i + 1; j < n + 1; ++j) {
        sum += A[i][j] * x[j];
    }

    return (b[i] - sum) * (1.0 / A[i][i]);
}


int relax_bottom(int n, vector<double>& x, vector<vector<double>>& A, vector<double>& b) {
    double h = 1.0 / n;
    double omega = 0.8;
    double error = 1.0;
    int k = 0;
    vector<double> new_x(n, 0.);

    while (error > 1.0 / pow(n, 3)) {
        for (int i = 0; i < n; ++i) {
            new_x[i] = Relax_calculate_new_x(i, new_x, n, A, omega, b);
        }
        x = new_x;
        error = calculate_error(new_x, A, b);
        k++;
    }
    return k;
}

double Relax_calculate_new_x(int i, vector<double>& x, int n, vector<vector<double>> &A, double omega, vector<double> &b) {
    
    double sum = 0.0;
    if (i > 0) {
        for (int j = 0; j <= i - 1; ++j) {
            sum += A[i][j] * x[j];
        }
    }
    for (int j = i + 1; j < n + 1; ++j) {
        sum += A[i][j] * x[j];
    }
    double new_x = (b[i] - sum) * (1.0 / A[i][i]);
    return x[i] + omega * (new_x - x[i]);
}


int spusk(int n, vector<double>& x, vector<vector<double>>& A,  vector<double>& b) {
    double h = 1.0 / n;
    int k = 1;
    vector<double> new_x(n, 0.);
    double error = 1.0;

    while (error > 1.0 / pow(n, 3)) {
        for (int i = 0; i < n; ++i) {
            new_x[i] = spusk_calculate_new_x(i, new_x, n, A, b);
        }
        x = new_x;
        error = calculate_error(new_x, A, b);
        k++;
    }
}

double spusk_calculate_new_x(int i,  vector<double>& x, int n,  vector<vector<double>>& A,  vector<double>& b) {
    vector<double> r = calculate_r(A, b, x);
    vector<double> Ar = spusk_calculate_matrix_vector_multiplication(A, r);

    return x[i] -  spusk_calculate_tau(r, Ar) * r[i];
}

vector<double> calculate_r( vector<vector<double>>& A,  vector<double>& b,  vector<double> &x) {
    vector<double> Ax = spusk_calculate_matrix_vector_multiplication(A, x);
    vector<double> r(Ax.size());

    for (size_t i = 0; i < r.size(); ++i) {
        r[i] = Ax[i] - b[i];
    }

    return r;
}

vector<double> spusk_calculate_matrix_vector_multiplication(vector<vector<double>>& A,  vector<double> &x) {
    int n = A.size();
    int m = x.size();
    vector<double> result(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            result[i] += A[i][j] * x[j];
        }
    }
    return result;
}

double spusk_calculate_tau( vector<double>& r,  vector<double>& Ar) {
    double a = 0.;
    double b = 0.;
    for (size_t i = 0; i < r.size(); ++i) {
        a += r[i] * r[i];
        b += Ar[i] * r[i];
    }
    if (abs(b) < 1e-10) {
        return 0.0;
    }
    return a * (1.0 / b);
}


vector<vector<double>> create_matrix(int n){
    double h = 1. / n;

    vector<vector<double>> matrix_res(n + 1, vector<double>(n + 1,0.0));

    for (int i = 0; i < n; ++i) {
        if (i == 0) {
            double b = (p_x((i + 1) * h) + p_x((i + 2) * h) + (pow(h,2) * q_x((i + 1) * h)));
            double c = -p_x((i + 2) * h);
            matrix_res[i][i] = b;
            matrix_res[i][i + 1] = c;
            continue;
        }
        if (i == (n - 1)) {
            double b = (p_x((i + 1) * h) + p_x((i + 2) * h) + (pow(h,2) * q_x((i + 1) * h)));
            double a = -p_x((i + 1) * h);
            matrix_res[i][i] = b;
            matrix_res[i][i - 1] = a;
            continue;
        }
        double a = -p_x((i + 1) * h);
        double b = (p_x((i + 1) * h) + p_x((i + 2) * h) + (pow(h,2) * q_x((i + 1) * h)));
        double c = -p_x((i + 2) * h);

        matrix_res[i][i] = b;
        matrix_res[i][i + 1] = c;
        matrix_res[i][i - 1] = a;
    }

    return matrix_res;
}

void m_print(int k, vector<double> &y_i, vector<double> &y_ik) {
    for (int i = 0; i < y_i.size() - 1; ++i) {
        printf("ih = %3d | y_i = %9.6lf | y_ik = %9.6lf | |y_i - y_ik| = %9.6lf | k = %d\n", i, y_i[i], y_ik[i], abs(y_i[i] - y_ik[i]), k);
    }
}


int main() {

    int n = 10.0;

    std::vector<double> y_i_p(n + 1), y_i_s(n + 1), y_i_r(n + 1), y_i_sp(n + 1), b(n);
    vector<vector<double>> A = create_matrix(n);
    b = func(n);
    printf("Прогонка\n");
    progonka_method(A, y_i_p, n, b);
    int k_s = Jacobi_method(n, y_i_s, A, b);
    int k_r = relax_bottom(n, y_i_r, A, b);
    int k_dec = spusk(n, y_i_sp, A, b);
    printf("Прогонка - Якоби\n");
    m_print(k_s, y_i_p, y_i_s);
    printf("Прогонка - Нижняя релаксация\n");
    m_print(k_r, y_i_p, y_i_r);
    printf("Прогонка - Наискорейший спуск\n");
    m_print(k_dec, y_i_p, y_i_sp);
    return 0;
}