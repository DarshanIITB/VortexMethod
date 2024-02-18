#pragma once
#ifndef UTILS
#define UTILS
#include <vector>

using std::vector;

float arctan(float& angle)
{
    return std::atan(angle);
}

float square(float& num)
{
    return num * num;
}

float trapz(const std::vector<float>& y, float dx)
{
    float integral = 0.0;
    for (size_t i = 1; i < y.size(); ++i)
    {
        integral += 0.5 * dx * (y[i] + y[i - 1]);
    }
    return integral;
}

float trapz(const std::vector<float>& y, const std::vector<float> x)
{
    float integral = 0.0f;
    for (size_t i = 1; i < y.size(); ++i) {
        //std::cout << "x[i]: " << x[i] << ", x[i-1]: " << x[i - 1] << ", y[i]: " << y[i] << ", y[i-1]: " << y[i - 1] << std::endl;
        integral += (x[i] - x[i - 1]) * (y[i] + y[i - 1]);
    }
    return integral * 0.5;
}

float h1_ig(const float L_ig, const float LHS, float mu, float alpha)
{
    return LHS - 2 * L_ig * sqrt(pow(mu, 2) + pow(mu * tan(abs(alpha)) + L_ig, 2));
}

float h2_ig(const float& L_ig, const float& LHS, const float& mu, const float& alpha)
{
    const float term2 = sqrt(pow(mu, 2) + pow(mu * tan(abs(alpha)) + L_ig, 2));
    const float diff = - (2 * term2 + (2 * L_ig * (mu * tan(abs(alpha)) + L_ig) / term2));
    return diff;
}

float h1_i(const float& L_i, const float& L_ig_calc, const float& mu, const float& angle, const float& radius, const float& R_tip)
{
    const float LHS = (L_i) / (L_ig_calc);
    const float RHS = 1 + ((1.333 * mu) / (1.2 * (L_i + L_ig_calc) + mu)) * ((radius * cos(angle)) / R_tip);
    return LHS - RHS;
}

float h2_i(const float& L_i, const float& L_ig_calc, const float& mu, const float& angle, const float& radius, const float& R_tip)
{
    return (1 / L_ig_calc) + ((1.333 * 1.2 * mu / pow(1.2 * (L_i + L_ig_calc) + mu, 2)) * radius * cos(angle) / R_tip);
}

vector<float> elementWiseMul(const vector<float>& a, const vector<float>& b)
{
    if (a.size() != b.size())
    {
        return {};
    }
    vector<float> prod;
    for (int i = 0; i < a.size(); i++)
    {
        prod.push_back(a[i] * b[i]);
    }
    return prod;
}

vector<float> linspace(float start, float end, int num) {
    vector<float> result(num);
    float step = (end - start) / (num - 1);

    for (int i = 0; i < num; ++i) {
        result[i] = start + i * step;
    }

    return result;
}

vector<float> vectorAdd(const vector<float>& a, const vector<float>& b)
{
    if (a.size() != b.size())
    {
        return {};
    }
    vector<float> sum;
    for (int i = 0; i < a.size(); i++)
    {
        sum.push_back(a[i] - b[i]);
    }
    return sum;
}

vector<float> vectorSub(vector<float>& a, vector<float>& b)
{
    if (a.size() != b.size())
    {
        return {};
    }
    vector<float> diff(a.size(), 0);
    for (int i = 0; i < a.size(); i++)
    {
        diff[i] = a[i] - b[i];
    }
    return diff;
}

//vector<vector<float>> scalarMul(const float&k, vector<float)

//vector<float> scalarMul(const float& k, const vector<float>& vec)
//{
//    vector<float> res;
//    for (auto& elem : vec)
//    {
//        res.push_back(k * elem);
//    }
//    return res;
//}

vector<float> vectorDiv(vector<float>& dividend, vector<float>& divisor)
{
    if (dividend.size() != divisor.size())
    {
        return {};
    }
    vector<float> res(dividend.size(), 0);
    for (int i = 0; i < dividend.size(); i++)
    {
        res[i] = dividend[i] / divisor[i];
    }
    return res;
}

template<typename T>
vector<T> scalarMul(const float& k, const vector<T>& vec)
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

vector<float> scalarAdd(const float& k, const vector<float>& vec)
{
    vector<float> res;
    for (auto& elem : vec)
    {
        res.push_back(k + elem);
    }
    return res;
}

std::pair<float, float> rotate(float x, float y, float theta)
{
    const float x_ = x * std::cos(theta) - y * std::sin(theta);
    const float y_ = x * std::sin(theta) + y * std::cos(theta);
    return { x_, y_ };
}

#endif