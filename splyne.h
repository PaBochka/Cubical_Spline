#ifndef SPLYNE_H
#define SPLYNE_H
#include <vector>

class Splyne
{
public:
    double h;
    Splyne();
    std::vector<double> TDMASolve(double myu1, double myu2, const std::vector<double> &A, std::vector<double> &C, int n, double h, double a, double b);
    double spline_s(double a, double d, double c, double b, double x, double ih);
    double func_fi(double x);

    double first_dev_fi(double x);
    double second_dev_fi(double x);

    double first_dev_spline_s(double d, double c, double b, double x, double ih);
    double second_dev_spline_s(double d, double c, double x, double ih);
};

#endif // SPLYNE_H
