#include "SimulatorMath.h"

double h = 10;
char buffer[250];


// Monte Carlo vars

int precision = 100;
double maxStdDev = 1.5;

// Quasi Monte Carlo Halton Sequence vars
int haltonSize = 100;
double halton2[100];
double halton3[100];
bool haltonsFilled = false;


/*
* Logs a double value with a given label. Uses a predeclared char buffer called buffer
*
* label: cstr of label to be printed
*
* d: double to be printed
*/
void _logDouble(const char* label, double d)
{
    sprintf(buffer, "%s: %f", label, d);
    puts(buffer);
}

//  Priamry Function 
double f(double x, double y)
{
    return h/pow((x*x+y*y+h*h), 1.5);
}

// Linear Test Function
double g(double x, double y)
{
    return 5*x+4*y;
}

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
double MidpointSumIntegral(Function2D func, struct SimulatorMathRect* bounds, int n)
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
            area += func(xStart + i*xStep, yStart + i*yStep);
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
double TrapezoidalSumIntegral(Function2D func, struct SimulatorMathRect* bounds, int n)
{
    double x1 = bounds->x1, x2 = bounds->x2, y1 = bounds->y1, y2 = bounds->y2;

    double xStep = (x2-x1)/n;
    double yStep = (y2-y1)/n;

    double area = 0;

    area += func(x1, y1);
    area += func(x2, y1);
    area += func(x1, y2);
    area += func(x2, y2);

    for (int i = 1; i < n; i++)
    {
        area += 2*func(x1 + i*xStep, y1);
        area += 2*func(x1 + i*xStep, y2);
        area += 2*func(x1, y1 + i*yStep);
        area += 2*func(x2, y1 + i*yStep);

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
double SimpsonsIntegral(Function2D func, struct SimulatorMathRect* bounds, int n)
{
    double x1 = bounds->x1, x2 = bounds->x2, y1 = bounds->y1, y2 = bounds->y2;
    double xStep = (x2-x1)/n;
    double yStep = (y2-y1)/n;

    double area = 0;

    area += func(x1, y1);
    area += func(x2, y1);
    area += func(x1, y2);
    area += func(x2, y2);

    int c = 2;

    for (int i = 1; i < n; i+=2)
    {
        area += 4*func(x1+i*xStep, y1);
        area += 4*func(x1+i*xStep, y2);
        area += 4*func(x1, y1+i*yStep);
        area += 4*func(x2, y1+i*yStep);

        for (int j = 1; j < n; j+=2)
        {
            area += 16*func(x1+i*xStep, y1+i*yStep);
        }

        for (int j = 2; j < n; j+=2)
        {
            area += 8*func(x1+i*xStep, y1+i*yStep);
        }

    }

    for (int i = 2; i < n; i+=2)
    {
        area += 2*func(x1+i*xStep, y1);
        area += 2*func(x1+i*xStep, y2);
        area += 2*func(x1, y1+i*yStep);
        area += 2*func(x2, y1+i*yStep);

        for (int j = 1; j < n; j+=2)
        {
            area += 8*func(x1+i*xStep, y1+i*yStep);
        }

        for (int j = 2; j < n; j+=2)
        {
            area += 4*func(x1+i*xStep, y1+i*yStep);
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
* maxStdDev: The maximum amount of standard deviations a value generated by this function will be allowed to fall from the center
* 
* returns: The randomly selected value based on the given parameters and the Gaussian distribution 
*/
double _genGaussRand(double lower, double upper, int precision, double maxStdDev)
{   
    // Ex precision of 100 gives rand [0, 100] gives [-50, 50] gives [-1.0, 1.0] with a step of .02 given by 2/(precision). This value is the amount of stdDevs we can go from the center  
    double stdDev = (rand() % (precision+1) - precision/2) / (double) (precision/2) * maxStdDev;
    
    double fittedValue = 1-exp(-1*stdDev*stdDev/2);

    double center = (upper+lower)/2;
    double range = (upper-lower);

    // center +- fittedValue * distance
    return center + ((stdDev > 0) - (stdDev < 0)) * fittedValue * range/2;


}

/*
* Performs 2 dimensional numerical integration on the given function using the Monte Carlo Method with n samples. Splits the computation into multiple cycles in order to handle large values of n.
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
double MonteCarloIntegral(Function2D func, struct SimulatorMathRect* bounds, int n)
{
    double x1 = bounds->x1, x2 = bounds->x2, y1 = bounds->y1, y2 = bounds->y2;
    long double workingAverage = 0;
    long double trueAverage = 0;

    int cycles = n % 100 + 1;

    for (int i = 1; i <= cycles; i++)
    {
        for (int j = (i-1)*n/cycles; j < i*n/cycles; j++)
        {
            workingAverage += func(_genGaussRand(x1, x2, precision, maxStdDev), _genGaussRand(y1, y2, precision, maxStdDev));
        }

        trueAverage += workingAverage /n * (x2-x1) * (y2-y1);
        workingAverage = 0;
    }

    return trueAverage;
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
* Pre-computes halton sequence values for bases 2 and 3 for use in the Quasi Monte Carlo Integral
*/
void fillHaltons()
{
    haltonsFilled = true;
    for (int i = 0; i < haltonSize; i++)
    {
        halton2[i] = _haltonSeq(i+1, 2);
        halton3[i] = _haltonSeq(i+1, 3);
    }
}
/*
* Performs 2 dimensional numerical integration on the given function using a Quasi Monte Carlo method selecting n values from the Hatlon sequence. Must ensure that the size of the Halton sequence arrays matches n. 
* 
* Time Complexity: O(n) after first call, O(nlog(n)) on first
*
* func: The function of 2 variables to be integrated
*
* bounds: The bounds of the integration region in the form of a pointer to a SimulatorMathRect
*
* n: The amount of rectnagles to be used
*
* returns: The approximated value of the integral
*/
double QuasiMonteCarloIntegral(Function2D func, struct SimulatorMathRect* bounds, int n)
{
    double x1 = bounds->x1, x2 = bounds->x2, y1 = bounds->y1, y2 = bounds->y2;
    double xRange = x2-x1, yRange = y2-y1;
    long double workingAverage = 0;
    long double trueAverage = 0;

    if (!haltonsFilled)
    {
        fillHaltons();
    }


    int cycles = n % 100 + 1;
    

    if (xRange > yRange)
    {
        for (int i = 0; i < n; i++)
        {
            for (int i = 1; i <= cycles; i++)
            {
                for (int j = (i-1)*n/cycles; j < i*n/cycles; j++)
                {
                    workingAverage += func(x1+halton2[i]*xRange, y1+halton3[i]*yRange);
                }

                trueAverage += workingAverage /n * (x2-x1) * (y2-y1);
                workingAverage = 0;
            }
        }
    }
    else
    {
        for (int i = 1; i <= cycles; i++)
        {
            for (int j = (i-1)*n/cycles; j < i*n/cycles; j++)
            {
                workingAverage += func(x1+halton2[i]*xRange, y1+halton3[i]*yRange);
            }

            trueAverage += workingAverage /n * (x2-x1) * (y2-y1);
            workingAverage = 0;
        }
    }

    return trueAverage;
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
* log: Whether or not to log the result
*
* returns: The time elapsed in milliseconds from before and after the integralFunc is called.
*/
double _timeFunc(NumericalIntegral integralFunc, Function2D func, struct SimulatorMathRect* bounds, int n, bool log)
{
    clock_t start = clock();
    double d = integralFunc(func, bounds, n);
    clock_t end = clock();

    if (log)
        _logDouble("Result", d);

    return ((double) (end - start)) / CLOCKS_PER_SEC * pow(10,3);
}

int main()
{
    srand((unsigned)time(NULL));

    struct SimulatorMathRect rect = {.x1 = 0, .x2 = 4, .y1 = 0, .y2 = 6};

    /*
    logDouble("Midpoint", MidpointSumIntegral(g, &rect, 10));

    logDouble("Trapezoidal", TrapezoidalSumIntegral(g, &rect, 10));

    logDouble("Simpsons", SimpsonsIntegral(g, &rect, 10));

    logDouble("Monte Carlo", MonteCarloIntegral(g, &rect, 10));

    fillHaltons();

    logDouble("Quasi Monte Carlo", QuasiMontiCarloIntegral(g, &rect, haltonSize));

    */

   _logDouble("Midpoint Time",_timeFunc(MidpointSumIntegral, g, &rect, 1000, true));
   _logDouble("Trapezoidal Time",_timeFunc(TrapezoidalSumIntegral, g, &rect, 1000, true));
   _logDouble("Simpsons Time",_timeFunc(SimpsonsIntegral, g, &rect, 1000, true));
   _logDouble("Monte Carlo Time", _timeFunc(MonteCarloIntegral, g, &rect, 1000 * 1000, true));
   _logDouble("Quasi Monte Carlo Time", _timeFunc(QuasiMonteCarloIntegral, g, &rect,1000*1000, true));

    return 0;
}