#ifndef SPLYNE_H
#define SPLYNE_H
#include <vector>

class Splyne
{
public:
    double h;
    Splyne();
    std::vector<double> TDMASolve(double myu1, double myu2, const std::vector<double> &A, std::vector<double> &C, int n, double h, double a, double b);
    double spline_s(const std::vector<double> &a, const std::vector<double> &d, const std::vector<double> &c, const std::vector<double> &b, double x, const std::vector<double> &X, int n);
    double func_fi(double x);

    double first_dev_fi(double x);
    double second_dev_fi(double x);

    double first_dev_spline_s(const std::vector<double> &d, const std::vector<double> &c, const std::vector<double> &b, double x, const std::vector<double> &X, int n);
    double second_dev_spline_s(const std::vector<double> &d, const std::vector<double> &c, double x, const std::vector<double> &X, int n);
};

#endif // SPLYNE_H
