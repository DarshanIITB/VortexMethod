#include "Coordinates.h"
#include <iostream>
#include <vector>
#include <cmath>

using std::sqrt;
using std::pow;
using std::sin;
using std::cos;

Coordinates::Coordinates(float xCoord, float yCoord, float zCoord) : x(xCoord), y(yCoord), z(zCoord) {}

float Coordinates::distanceTo(const Coordinates& other) const
{
    return sqrt(pow(x - other.x, 2) + pow(y - other.y, 2) + pow(z - other.z, 2));
}

void Coordinates::rotateIn2D(float angle)
{
    float xNew = x * cos(angle) - y * sin(angle);
    float yNew = x * sin(angle) + y * cos(angle);
    x = xNew;
    y = yNew;
}
