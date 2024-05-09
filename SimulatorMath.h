#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#define PI 3.14159265358979323846

struct SimulatorMathRect
{
    double x1; double x2; double y1; double y2;
};


typedef double (*Function2D)(double, double);

typedef double (*NumericalIntegral)(Function2D, struct SimulatorMathRect*, int);



