#pragma once
#ifndef COORDINATES_H
#define COORDINATES_H
#include <vector>

class Coordinates {
public:
    double x;
    double y;
    double z;
    // Constructor
    Coordinates(double xCoord, double yCoord, double zCoord);

    double distanceTo(const Coordinates& other) const;

    void rotateIn2D(double angle);

};

#endif // COORDINATES_H
