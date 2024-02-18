#include "Coordinates.h"
#include <iostream>
#include <vector>
#include <cmath>

using std::sqrt, std::pow, std::sin, std::cos;

Coordinates::Coordinates(double xCoord, double yCoord, double zCoord) : x(xCoord), y(yCoord), z(zCoord) {}

double Coordinates::distanceTo(const Coordinates& other) const
{
    return sqrt(pow(x - other.x, 2) + pow(y - other.y, 2) + pow(z - other.z, 2));
}

void Coordinates::rotateIn2D(double angle)
{
    double xNew = x * cos(angle) - y * sin(angle);
    double yNew = x * sin(angle) + y * cos(angle);
    x = xNew;
    y = yNew;
}
