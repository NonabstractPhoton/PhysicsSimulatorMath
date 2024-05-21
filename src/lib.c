#include "SimulatorMath.h"

// Monte Carlo vars
int precision = 1000;

// Quasi Monte Carlo Halton Sequence vars
double* halton2 = NULL;
double* halton3 = NULL;
int haltonSize = 0;

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
double MidpointSumIntegral2D(Function2D func, struct SimulatorMathRect* bounds, int n)
{    
    double x1 = bounds->x1, x2 = bounds->x2, y1 = bounds->y1, y2 = bounds->y2;

    double xStep = (x2-x1)/n;
    double yStep = (y2-y1)/n;

    double xStart = x1 + xStep/2;
    double yStart = y1 + yStep/2;

    double area = 0;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            area += func(xStart + i*xStep, yStart + j*yStep);
        }
    }

    return area * xStep * yStep;
}

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
double TrapezoidalSumIntegral2D(Function2D func, struct SimulatorMathRect* bounds, int n)
{
    double x1 = bounds->x1, x2 = bounds->x2, y1 = bounds->y1, y2 = bounds->y2;

    double xStep = (x2-x1)/n;
    double yStep = (y2-y1)/n;

    double area = 0;

    // Corners
    area += func(x1, y1)
    + func(x2, y1)
    + func(x1, y2)
    + func(x2, y2);

    for (int i = 1; i < n; i++)
    {
        // Sides
        area += 2 * (func(x1 + i*xStep, y1)
        + func(x1 + i*xStep, y2)
        + func(x1, y1 + i*yStep)
        + func(x2, y1 + i*yStep));

        // Inner 
        for (int j = 1; j < n; j++)
        {
            area += 4*func(x1 + i*xStep, y1 + j*yStep);
        }
    }

    return area * xStep * yStep / 4;      
}

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
double SimpsonsIntegral2D(Function2D func, struct SimulatorMathRect* bounds, int n)
{
    double x1 = bounds->x1, x2 = bounds->x2, y1 = bounds->y1, y2 = bounds->y2;
    double xStep = (x2-x1)/n;
    double yStep = (y2-y1)/n;

    double area = 0;

    // Corners
    area += func(x1, y1)
    + func(x2, y1)
    + func(x1, y2)
    + func(x2, y2);

    // Odd x Terms
    for (int i = 1; i < n; i+=2)
    {
        // Sides
        area += 4 * (func(x1+i*xStep, y1)
        + func(x1+i*xStep, y2)
        + func(x1, y1+i*yStep)
        + func(x2, y1+i*yStep));

        // Odd y terms
        for (int j = 1; j < n; j+=2)
        {
            area += 16*func(x1+i*xStep, y1+j*yStep);
        }

        // Even y terms
        for (int j = 2; j < n; j+=2)
        {
            area += 8*func(x1+i*xStep, y1+j*yStep);
        }

    }

    // Even x terms
    for (int i = 2; i < n; i+=2)
    {
        // Sides
        area += 2 * (func(x1+i*xStep, y1)
        + func(x1+i*xStep, y2)
        + func(x1, y1+i*yStep)
        + func(x2, y1+i*yStep));

        // Odd y terms
        for (int j = 1; j < n; j+=2)
        {
            area += 8*func(x1+i*xStep, y1+j*yStep);
        }

        // Even y terms
        for (int j = 2; j < n; j+=2)
        {
            area += 4*func(x1+i*xStep, y1+j*yStep);
        }
    }

    return area * xStep * yStep / 9;
}
/*
* Generates a random number from a Gaussian distribution that can vary at most by a given amount of standard deviations from the mean.
* More precisely, this function returns mean +- fittedRandomValue * range, where fittedRandomValue is given by 1-exp(-1*stdDev^2), with stdDev randomly generated according to the given parameters.
*
* lower: The lower bound of the distribution
*
* upper: The upper bound of the distribution
*
* precision: A value representing the precision of the random values used in the fitting process. The step between potential standard deviations from the center is given by (double) 2/precision.
*
* 
* returns: The randomly selected value based on the given parameters and the Gaussian distribution 
*/
double _genGaussRand(double lower, double upper, int precision, int samples)
{   
    // Defaults at ~2.87
    double maxStdDev =  (1.25 * (log(upper-lower)/log(10) + log(samples)/log(2))) - 8.75;

    // Ex precision of 100 gives rand [0, 100] gives [-50, 50] gives [-1.0, 1.0] with a step of .02 given by 2/(precision). This value is the amount of stdDevs we can go from the center  
    double stdDev = (rand() % (precision+1) - precision/2) / (double) (precision/2) * maxStdDev;
    
    double fittedValue = 1-exp(-1*stdDev*stdDev/2);

    double center = (upper+lower)/2;
    double range = (upper-lower);

    // center +- fittedValue * distance
    return center + ((stdDev > 0) - (stdDev < 0)) * fittedValue * range/2;


}

/*
* Sets the precision for the Monte Carlo Numerical Integration. 
*
* p: The inverse of the noramlized difference between adjacent possible results in random value generation.
*/
void SetMonteCarloPrecision(int p)
{
    precision = p;
}

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
double GaussianMonteCarloIntegral2D(Function2D func, struct SimulatorMathRect* bounds, int n)
{
    double x1 = bounds->x1, x2 = bounds->x2, y1 = bounds->y1, y2 = bounds->y2;
    double workingAverage = 0;
    long double trueAverage = 0;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < 50; j++)
        {
            workingAverage += func(_genGaussRand(x1, x2, precision, n), _genGaussRand(y1, y2, precision, n));
        }

        trueAverage += workingAverage / 100;
        workingAverage = 0;
    }

    return trueAverage / n * (x2-x1) * (y2-y1);
}


/*
/* Generates a pseudo random number for use in the pseudorandom Monte-Carlo Numerical Integration
*
* lower: Lower bound
*
* upper: Upper bound
*
* precision: The inverse of the noramlized difference between adjacent possible results. More precisely, each generated value differs by (upper-lower)/precision from its adjacent possible value. 
*
* returns: The pseudorandom value according to the specifications
*/
double _genPseudoRand(double lower, double upper, int precision)
{
    return lower + (rand() % (int)(precision * (upper - lower + 1)))/precision;
}

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
double PsuedoRandMonteCarloIntegral2D(Function2D func, struct SimulatorMathRect* bounds, int n)
{
    double x1 = bounds->x1, x2 = bounds->x2, y1 = bounds->y1, y2 = bounds->y2;
    double workingAverage = 0;
    long double trueAverage = 0;


    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < 100; j++)
        {
            workingAverage += func(_genPseudoRand(x1, x2, precision), _genPseudoRand(y1, y2, precision));
        }

        trueAverage += workingAverage / 100;
        workingAverage = 0;
    }

    return trueAverage / n * (x2-x1) * (y2-y1);
}

/*
* Computes the value at a given index for the halton sequence of the provided base
*/
double _haltonSeq(int index, int base){
    
    double x = 1, y = 0;
    
    while (index > 0)
    {
        x /= base;
        y += x * (index % base);
        index /= base;
    }
    return y;
}
/*
* Pre-computes halton sequence values for bases 2 and 3 for use in the Quasi Monte Carlo Integral.
*
* Must be called before use the Quasi Monte Carlo Integral for max effeciency on repeated use.
*
* Must be freed with freeHaltons after use.
*/
void FillHaltons(int n)
{
    if (haltonSize == n)
        return;

    haltonSize = n;

    halton2 = (double*)malloc(sizeof(double) * n);
    halton3 = (double*)malloc(sizeof(double) * n);

    for (int i = 0; i < haltonSize; i++)
    {
        halton2[i] = _haltonSeq(i+1, 2);
        halton3[i] = _haltonSeq(i+1, 3);
    }
}

void FreeHaltons()
{
    free(halton2);
    free(halton3);
    haltonSize = 0;
}
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
double QuasiMonteCarloIntegral2D(Function2D func, struct SimulatorMathRect* bounds, int n)
{
    double x1 = bounds->x1, x2 = bounds->x2, y1 = bounds->y1, y2 = bounds->y2;
    double xRange = x2-x1, yRange = y2-y1;
    double workingAverage = 0;
    double trueAverage = 0;

    bool autoFillHaltons = haltonSize == 0;
    
    if (autoFillHaltons)
        FillHaltons(n);


    int cycles = n / 50 + 1;
    
    int index = 0;

    if (xRange > yRange)
    {
        for (int i = 1; i <= cycles; i++)
        {
            for (int j = (i-1)*n/cycles; j < i*n/cycles; j++)
            {
                workingAverage += func(x1+halton2[index]*xRange, y1+halton3[index]*yRange);
                index++;
            }

            trueAverage += workingAverage /n;
            workingAverage = 0;
        }
    }
    else
    {
        for (int i = 1; i <= cycles; i++)
        { 
            for (int j = (i-1)*n/cycles; j < i*n/cycles; j++)
            {
                workingAverage += func(x1+halton3[index]*xRange, y1+halton2[index]*yRange);
                index++;
            }

            trueAverage += workingAverage /n;
            workingAverage = 0;
        }
    }

    if (autoFillHaltons)
        FreeHaltons();

    return trueAverage * xRange * yRange;
}


/*
* Times a NumericalIntegral in milliseconds, optionally logging the result.
*
* integralFunc: The NumericalIntegral to be timed
*
* func: The Function2D which the NumericalIntegral uses
*
* bounds: The SimulatorMathRect that descirbes integration region bounds
*
* n: The number of samples the integralFunc will use
*
* rresult: The double ptr to store the result of the integral
*
* returns: The time elapsed in milliseconds from before and after the integralFunc is called.
*/
double _timeFunc(NumericalIntegral integralFunc, Function2D func, struct SimulatorMathRect* bounds, int n, double* result)
{
    clock_t start = clock();
    *result = integralFunc(func, bounds, n);
    clock_t end = clock();

    return ((double) (end - start)) / (double)CLOCKS_PER_SEC * pow(10,3);
}