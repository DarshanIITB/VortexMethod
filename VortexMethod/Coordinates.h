#pragma once
#ifndef COORDINATES_H
#define COORDINATES_H
#include <vector>

class Coordinates {
public:
    float x;
    float y;
    float z;
    // Constructor
    Coordinates(float xCoord, float yCoord, float zCoord);

    float distanceTo(const Coordinates& other) const;

    void rotateIn2D(float angle);

};

#endif // COORDINATES_H
