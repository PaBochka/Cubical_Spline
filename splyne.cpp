#include "splyne.h"
#include <cmath>
#include <vector>
#include <iostream>
Splyne::Splyne()
{
}
double Splyne::spline_s(const std::vector<double> &a, const std::vector<double> &d, const std::vector<double> &c, const std::vector<double> &b, double x, const std::vector<double> &X, int n)
{
    for (int i = 1; i <= n; i++) {
        if(x >= X[i - 1] && x <= X[i])
            return (a[i] + b[i] * (x - X[i]) + (c[i] / 2) * pow((x - X[i]), 2) + (d[i] / 6) * pow((x - X[i]), 3));
    }
    return 0.0;
}
double Splyne::func_fi(double x)
{
    double res;
    if((x <= 0) && (x >= -1))
    {
        res = pow(x, 3) + 3 * pow(x, 2);
    }
    else
    {
        res = -pow(x, 3) + 3 * pow(x, 2);
    }
    return res;
}

double Splyne::first_dev_fi(double x)
{
    double res;
    if((x <= 0) && (x >= -1))
    {
        res = 3 * pow(x, 2) + 6 * x;
    }
    else
    {
        res = -3 * pow(x, 2) + 6 * x;
    }
    return res;
}
double Splyne::second_dev_fi(double x)
{
    double res;
    if((x <= 0) && (x >= -1))
    {
        res = 6 * x + 6;
    }
    else
    {
        res = -6 * x + 6;
    }
    return res;
}
double Splyne::first_dev_spline_s(const std::vector<double> &d, const std::vector<double> &c, const std::vector<double> &b, double x, const std::vector<double> &X, int n)
{
    for (int i =0; i <= n; i++) {
        if(x >= X[i - 1] && x <= X[i])
            return (b[i] + (c[i] / 2.0) * 2.0 * (x - X[i]) + (d[i] / 2.0) * pow((x - X[i]), 2));
    }
    return 0.0;
}
double Splyne::second_dev_spline_s(const std::vector<double> &d, const std::vector<double> &c, double x, const std::vector<double> &X, int n)
{
    for (int i = 0; i <= n; i++) {
        if(x >= X[i - 1] && x <= X[i])
            return (c[i] + d[i] * (x - X[i]));
    }
    return 0.0;
}
std::vector<double> Splyne::TDMASolve(double myu1, double myu2, const std::vector<double> &A, std::vector<double> &C, int n, double h, double a, double b)
{
    std::vector<double> Alpha(n + 1);
    std::vector<double> Beta(n + 1);
    Alpha[1]=0;
    Beta[1]=myu1;
    for (int i = 1; i <=n-1; i++)
    {
    Alpha[i + 1] = h / (-4 * h - Alpha[i] * h);
    Beta[i + 1] = ((-6.0 / h) * (A[i + 1] - 2 * A[i] + A[i - 1]) + Beta[i] * h) / (-4 * h - Alpha[i] * h);
    }
    C[n]=myu2;
    for (int i = n ; i >= 1; i--)
    {
    C[i - 1] = Alpha[i] * C[i] + Beta[i];
    }
    return C;
}

