#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#define PI 3.14159265358979323846

/*
* A struct containing the rectnagular bounds of a 2D integration region
*
* x1: Smaller 1st dimension bound
*
* x2: Larger 1st dimension bound
*
* y1: Smaller 2nd dimension bound
*
* y2: Larger 2nd dimension bound 
*/
struct SimulatorMathRect
{
    double x1; double x2; double y1; double y2;
};

// Represents a mathemtatical function of two variables
typedef double (*Function2D)(double, double);

/*
* Represents a numerical integration procedure.
*
* Inputs:
*
* A Function2D to be integrated
*
* A SimulatorMathRect containing integration bounds
*
* An integer representing sample count for the numerical intergation procedure
*/
typedef double (*NumericalIntegral)(Function2D, struct SimulatorMathRect*, int);



