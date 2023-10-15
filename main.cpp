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

void Tabulation() {
/*  
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

void ConstructionOfNewtonInterpolationPolynomial(){
/*
    По полученной таблице значений посторить интерполяционный полином
    Ньютона, приближающий Si(x)
    L_n(x) = f(x_0)+(x-x_o)f(x_0,x_1)+(x-x_0)(x-x_1)f(x_0,x_1,x_2)+...
    +(x-x_0)(x-x_1)...(x-x_n)f(x_0,x_1,...,x_n)
    и вычислить погрешность интерполирования
    e_n = max(e(x)) (x∈(a,b)), e(x) = |Si(x)-L_n(x)|

*/

}



int main() {
    Tabulation();
    return 0;
}
