#include <iostream>
#include <cmath>
#include <iomanip>

double Summ(double x) {
    const double eps = 0.000001;
    int i = 0;
    double s = 0;
    double ai = x;
    double qn = 0;

    while (std::fabs(ai) >= eps) {
        s = s + ai;
        qn = (-x * x * (i + 1)) / ((2 * i + 3) * (2 * i + 3) * (2 * i + 2));
        ai = ai * qn;
        i++;
    }

    return s;
}

// Определение функции Si(x) с использованием заданного ряда
    double Si(double x, double eps = 1e-6) {
    int n = 0;
    double term = x;
    double result = term;

    while (std::fabs(term) >= eps) {
        n++;
        term = (-term * x * x) / ((2 * n + 1) * (2 * n + 1));
        result += term;
    }

    return result;
}

// Функция для вычисления полинома Лагранжа
double LagrangeInterpolation(double x, int n, double xi[]) {
    double result = 0;
    for (int i = 0; i <= n; i++) {
        double term = Si(xi[i]);
        for (int j = 0; j <= n; j++) {
            if (j != i) {
                term *= (x - xi[j]) / (xi[i] - xi[j]);
            }
        }
        result += term;
    }
    return result;
}

void LagrangePolynomial() {
    /*
    Задание 2
    По полученной таблице значений построить интерполяционный полином Лагранжа, приближающий Si(x):  
    L_n = (сумма S(x_i) от i=0 до n)(произведение ((x-x_j)/(x_i - x_j) от j=0, j!=i до n))
    Затем вычислить погрешность интерполирования
    
    */
    double a = 0;
    double b = 4;
    double h = 0.3;
    int n = static_cast<int>((b - a) / h);

    double xi[100];
    for (int i = 0; i <= n; i++) {
        xi[i] = a + i * h;
    }

    std::cout << "x_i\t\tSi(x)\t\tL_n(x)\t\t e(x)\n";

    for (int i = 0; i <= n; i++) {
        double x = xi[i];
        double Si_x = Si(x);
        double L_n_x = LagrangeInterpolation(x, i, xi);
        double error = std::abs(Si_x - L_n_x);

        std::cout << std::fixed << std::setprecision(6)
                  << x << "\t\t" << Si_x << "\t\t" << L_n_x << "\t\t" << error << "\n";
    }
}

void Tabulation() {
/*  
    Задание 1
    Протабулировать Si(x) на отрезке [a,b] с шагом h с точностью e, основываясь на ряде Тейлора, 
    предварительно вычислив его
    Si(x)=сумма от n=0 до бесконечности ((-1)^n*x^(2n+1))/((2n+1)*(2n+1)!),
    где a=0, b=4, h=0,3, e=10^-6 и получить, таким образом, таблицу.
    x_0 x_1 x_2 ... x_n
    f_0 f_1 f_2 ... f_n
    f_1=Si(x), x_i=a+i*h, i=0,…,n.        
*/
    double a = 0; 
    double b = 4; 
    double h = 0.3; 
    double e = 1e-6; 

    std::cout << std::left << std::setw(10) << "x_i";
    for (double x = a; x <= b; x += h) {
        std::cout << std::left << std::setw(10) << x;
    }
    std::cout << std::endl;

    std::cout << std::left << std::setw(10) << "f_i";
    for (double x = a; x <= b; x += h) {
        double result = Summ(x);
        std::cout << std::left << std::setw(10) << std::fixed << std::setprecision(6) << result;
    }
    std::cout << std::endl;
}

int main() {
    Tabulation();
    LagrangePolynomial(); 
    return 0;
}

