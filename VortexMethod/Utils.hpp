#pragma once
#ifndef UTILS
#define UTILS
#include <vector>

using std::vector;

double round_(const double& num, const int n_places)
{
    return round(num * pow(10, n_places)) / pow(10, n_places);
}

double arctan(double& angle)
{
    return std::atan(angle);
}

double square(double& num)
{
    return num * num;
}

double trapz(const std::vector<double>& y, double dx)
{
    double integral = 0.0;
    for (size_t i = 1; i < y.size(); ++i)
    {
        integral += 0.5 * dx * (y[i] + y[i - 1]);
    }
    return integral;
}

double trapz(const vector<double>& y, const vector<double> x)
{
    double integral = 0.0f;
    for (size_t i = 1; i < y.size(); ++i) {
        //std::cout << "x[i]: " << x[i] << ", x[i-1]: " << x[i - 1] << ", y[i]: " << y[i] << ", y[i-1]: " << y[i - 1] << std::endl;
        integral += (x[i] - x[i - 1]) * (y[i] + y[i - 1]);
    }
    return integral * 0.5;
}

double trapz(const vector<vector<double>>& y, int second_index, vector<double>& x)
{
    double integral = 0.0f;
    for (size_t i = 1; i < y.size(); ++i) {
        //std::cout << "x[i]: " << x[i] << ", x[i-1]: " << x[i - 1] << ", y[i]: " << y[i] << ", y[i-1]: " << y[i - 1] << std::endl;
        integral += (x[i] - x[i - 1]) * (y[i][second_index] + y[i - 1][second_index]);
    }
    return integral * 0.5;
}

double h1_ig(const double L_ig, const double LHS, double mu, double alpha)
{
    return LHS - 2 * L_ig * sqrt(pow(mu, 2) + pow(mu * tan(abs(alpha)) + L_ig, 2));
}

double h2_ig(const double& L_ig, const double& LHS, const double& mu, const double& alpha)
{
    const double term2 = sqrt(pow(mu, 2) + pow(mu * tan(abs(alpha)) + L_ig, 2));
    const double diff = - (2 * term2 + (2 * L_ig * (mu * tan(abs(alpha)) + L_ig) / term2));
    return diff;
}

double h1_i(const double& L_i, const double& L_ig_calc, const double& mu, const double& angle, const double& radius, const double& R_tip)
{
    const double LHS = (L_i) / (L_ig_calc);
    const double RHS = 1 + ((1.333 * mu) / (1.2 * (L_i + L_ig_calc) + mu)) * ((radius * cos(angle)) / R_tip);
    return LHS - RHS;
}

double h2_i(const double& L_i, const double& L_ig_calc, const double& mu, const double& angle, const double& radius, const double& R_tip)
{
    return (1 / L_ig_calc) + ((1.333 * 1.2 * mu / pow(1.2 * (L_i + L_ig_calc) + mu, 2)) * radius * cos(angle) / R_tip);
}

double sumVecElem(const vector<double>& vec)
{
    double sum = 0.0f;
    for (double el : vec) sum += el;
    return sum;
}

vector<double> elementWiseMul(const vector<double>& a, const vector<double>& b)
{
    if (a.size() != b.size())
    {
        return {};
    }
    vector<double> prod;
    for (int i = 0; i < a.size(); i++)
    {
        prod.push_back(a[i] * b[i]);
    }
    return prod;
}

vector<double> linspace(double start, double end, int num) {
    vector<double> result(num);
    double step = (end - start) / (num - 1);

    for (int i = 0; i < num; ++i) {
        result[i] = start + i * step;
    }

    return result;
}

vector<double> vectorAdd(const vector<double>& a, const vector<double>& b)
{
    if (a.size() != b.size())
    {
        return {};
    }
    vector<double> sum;
    for (int i = 0; i < a.size(); i++)
    {
        sum.push_back(a[i] - b[i]);
    }
    return sum;
}

vector<double> vectorSub(vector<double>& a, vector<double>& b)
{
    if (a.size() != b.size())
    {
        return {};
    }
    vector<double> diff(a.size(), 0);
    for (int i = 0; i < a.size(); i++)
    {
        diff[i] = a[i] - b[i];
    }
    return diff;
}

vector<vector<double>> vectorSub(vector<vector<double>>& a, vector<vector<double>>& b)
{
    if (a.size() != b.size()) return {};
    vector<vector<double>> diff;
    for (int i = 0; i < a.size(); i++)
    {
        diff.push_back(vectorSub(a[i], b[i]));
    }
    return diff;
}

//vector<vector<double>> scalarMul(const double&k, vector<double)

//vector<double> scalarMul(const double& k, const vector<double>& vec)
//{
//    vector<double> res;
//    for (auto& elem : vec)
//    {
//        res.push_back(k * elem);
//    }
//    return res;
//}

vector<double> vectorDiv(vector<double>& dividend, vector<double>& divisor)
{
    if (dividend.size() != divisor.size())
    {
        return {};
    }
    vector<double> res(dividend.size(), 0);
    for (int i = 0; i < dividend.size(); i++)
    {
        res[i] = dividend[i] / divisor[i];
    }
    return res;
}

template<typename T>
vector<T> scalarMul(const double& k, const vector<T>& vec)
{
    std::vector<T> res;
    for (const auto& elem : vec)
    {
        if constexpr (std::is_arithmetic_v<T>) {
            res.push_back(k * elem);
        }
        else {
            res.push_back(scalarMul(k, elem));
        }
    }
    return res;
}

vector<double> scalarAdd(const double& k, const vector<double>& vec)
{
    vector<double> res;
    for (auto& elem : vec)
    {
        res.push_back(k + elem);
    }
    return res;
}

std::pair<double, double> rotate(double x, double y, double theta)
{
    const double x_ = x * std::cos(theta) - y * std::sin(theta);
    const double y_ = x * std::sin(theta) + y * std::cos(theta);
    return { x_, y_ };
}

#endif