#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#define PI 3.14159265358979323846

/*
* A struct containing the rectangular bounds of a 2D integration region. 
* This rectangle has corners (x1,y1), (x1,y2), (x2,y1), and (x2,y2). 
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
* Performs 2 dimensional numerical integration on the given function using a Midpoint Riemann Sum with n rectangles of equal area spanning the integration bounds.
* 
* Time Complexity: O(n^2)
* 
*func: The function of 2 variables to be integrated
*
* bounds: The bounds of the integration region in the form of a pointer to a SimulatorMathRect
*
* n: The amount of rectnagles to be used
*
* returns: The approximated value of the integral
*/
__declspec( dllexport ) double MidpointSumIntegral2D(Function2D func, const struct SimulatorMathRect* bounds, int n);

/*
* Performs 2 dimensional numerical integration on the given function using an extension of the Trapezoidal Rule with n shapes spanning the integration bounds.
* 
* Time Complexity: O(n^2)
*
* func: The function of 2 variables to be integrated
*
* bounds: The bounds of the integration region in the form of a pointer to a SimulatorMathRect
*
* n: The amount of rectnagles to be used
*
* returns: The approximated value of the integral
*/
__declspec( dllexport ) double TrapezoidalSumIntegral2D(Function2D func, const struct SimulatorMathRect* bounds, int n);

/*
* Performs 2 dimensional numerical integration on the given function using an extension of Simpson's 3/8 Rule with n shapes spanning the integration bounds.
* 
* Time Complexity: O(n^2)
*
* func: The function of 2 variables to be integrated
*
* bounds: The bounds of the integration region in the form of a pointer to a SimulatorMathRect
*
* n: The amount of rectnagles to be used
*
* returns: The approximated value of the integral
*/
__declspec( dllexport ) double SimpsonsIntegral2D(Function2D func, const struct SimulatorMathRect* bounds, int n);

/*
* Performs 2 dimensional numerical integration on the given function using the Monte Carlo Method with a modified gaussian probability distribution.
* Uses the Central Limit Theorem with n samples of size 50 to approximate true mean result.
* 
*
* Time Complexity: O(n)
*
* func: The function of 2 variables to be integrated
*
* bounds: The bounds of the integration region in the form of a pointer to a SimulatorMathRect
*
* n: The amount of rectangles to be used
*
* returns: The approximated value of the integral
*/
__declspec( dllexport ) double GaussianMonteCarloIntegral2D(Function2D func, const struct SimulatorMathRect* bounds, int n);

/*
* Performs 2 dimensional numerical integration on the given function using the Monte Carlo Method with a pseudo-random probability distribution.
* Uses the Central Limit Theorem with n samples of size 100 to approximate true mean result.
* 
*
* Time Complexity: O(n)
*
* func: The function of 2 variables to be integrated
*
* bounds: The bounds of the integration region in the form of a pointer to a SimulatorMathRect
*
* n: The amount of rectnagles to be used
*
* returns: The approximated value of the integral
*/
__declspec( dllexport ) double PsuedoRandMonteCarloIntegral2D(Function2D func, const struct SimulatorMathRect* bounds, int n);

/*
* Pre-computes halton sequence values for bases 2 and 3 for use in the Quasi Monte Carlo Integral.
*
* Must be called before use the Quasi Monte Carlo Integral for max effeciency on repeated use.
*
* Must be freed with FreeHaltons after use.
*/
__declspec( dllexport ) void FillHaltons(int n);

/*
* Frees the halton sequences generated by FillHaltons
*/
__declspec( dllexport ) void FreeHaltons();

/*
* Performs 2 dimensional numerical integration on the given function using a Quasi Monte Carlo method 
* selecting n values from the Hatlon sequence. Must ensure that the size of the Halton sequence arrays matches n if pre-filled. 
* Uses the Central Limit Theorem with samples of size 50, thus also handling large values of n.
* It is highly recommended to use fillHaltons before the execution of this algorithm and freeHaltons after for max effeciency.
* 
* Time Complexity: O(n) if fillHaltons is properly called, O(n^2 log(n)) otherwise
*
* func: The function of 2 variables to be integrated
*
* bounds: The bounds of the integration region in the form of a pointer to a SimulatorMathRect
*
* n: The amount of rectnagles to be used
*
* returns: The approximated value of the integral
*/
__declspec( dllexport ) double QuasiMonteCarloIntegral2D(Function2D func, const struct SimulatorMathRect* bounds, int n);