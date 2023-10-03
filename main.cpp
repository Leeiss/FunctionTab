//Протабулировать Si(x) на отрезке [a,b] с шагом h с точностью e, основываясь на ряде Еейлора, 
//предварительно вычислив его
//Si(x)=сумма от n=0 до бесконечности ((-1)^n*x^(2n+1))/((2n+1)*(2n+1)!),
//где a=0, b=4, h=0,3, e=10^-6 и получить, таким образом, таблицу.
// x_0 x_1 x_2 ... x_n
// f_0 f_1 f_2 ... f_n
//f_1=Si(x), x_i=a+i*h, i=0,…,n.


#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

void Printing(vector<double> arr_f, vector<double> arr_x) {
    for (int i = 0; i < 5; i++) {
        cout << arr_x[i] << ":   " << arr_f[i] << endl;
    }
}

void FillX(vector<double>* arr_x, double a, double h) {
    for (double i = 0; i < 5; i++) {
        arr_x->push_back(a + i * h);
    }
}

double GetQ(double x, int n) {
    return 0 - ((pow(x, 2) * (2.0 * n + 1.0)) / ((n + 1.0) * (2.0 * n + 3.0)));
}

double GetRes(double x)
{
    double e = pow(10, -6);

    double q_current = 0 - (pow(x, 2) / 3.0);

    double a_current = x;

    double summ = 0;
    
    summ += a_current ; // summ = a0

    a_current = a_current * q_current; // a1

    int i = 2; // start to count a2 = (a1 * q1)
    while (abs(a_current) >= e)
    {
        summ += a_current;  // int the first iteration we add a1
        a_current = GetQ(x, i) * a_current;  // int the first iteration we count a2
    }

    return summ;
};

int main()
{
    int a = 0;
    int b = 4;

    const int n = 5;

    const double h = 0.3;

    vector<double> arr_x;
    vector<double> arr_f;

    FillX(&arr_x, a, h);

    for (int x = 0; x < n; x++) {
        arr_f.push_back(GetRes(arr_x[x]));
    }

    for (int i = 0; i < n; i++) {
        arr_f[i] *= (2.0 / pow(3.14, 0.5));
    }

    Printing(arr_f, arr_x);

}
