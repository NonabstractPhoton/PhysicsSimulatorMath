#include "SimulatorMath.h"

double z = 10;
char buffer[250];


// Monte Carlo vars

int precision = 100;
double maxStdDev = 3;

// Quasi Monte Carlo Halton Sequence vars
double* halton2 = NULL;
double* halton3 = NULL;
int haltonSize = 0;


/*
* Internal Function. Logs a double value with a given label. Uses a predeclared char buffer called buffer
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

//  Internal Primary Function 1 
double f1(double x, double y)
{
    return z/pow((x*x+y*y+z*z), 1.5);
}

// Internal Primary Function 2
double f2(double x, double y)
{
    return exp(-1 * (x*x + y*y));
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
    area += func(x1, y1);
    area += func(x2, y1);
    area += func(x1, y2);
    area += func(x2, y2);

    for (int i = 1; i < n; i++)
    {
        // Sides
        area += 2*func(x1 + i*xStep, y1);
        area += 2*func(x1 + i*xStep, y2);
        area += 2*func(x1, y1 + i*yStep);
        area += 2*func(x2, y1 + i*yStep);

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
    area += func(x1, y1);
    area += func(x2, y1);
    area += func(x1, y2);
    area += func(x2, y2);

    // Odd x Terms
    for (int i = 1; i < n; i+=2)
    {
        // Sides
        area += 4*func(x1+i*xStep, y1);
        area += 4*func(x1+i*xStep, y2);
        area += 4*func(x1, y1+i*yStep);
        area += 4*func(x2, y1+i*yStep);

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
        area += 2*func(x1+i*xStep, y1);
        area += 2*func(x1+i*xStep, y2);
        area += 2*func(x1, y1+i*yStep);
        area += 2*func(x2, y1+i*yStep);

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
* Performs 2 dimensional numerical integration on the given function using the Monte Carlo Method with n samples.
* Uses the Central Limit Theorem with samples of size ~50, thus also handling large values of n.
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
double MonteCarloIntegral2D(Function2D func, struct SimulatorMathRect* bounds, int n)
{
    double x1 = bounds->x1, x2 = bounds->x2, y1 = bounds->y1, y2 = bounds->y2;
    double workingAverage = 0;
    double trueAverage = 0;

    int cycles = n / 50 + 1;


    for (int i = 1; i <= cycles; i++)
    {
        for (int j = (i-1)*n/cycles; j < i*n/cycles; j++)
        {
            workingAverage += func(_genGaussRand(x1, x2, precision, maxStdDev), _genGaussRand(y1, y2, precision, maxStdDev));
        }

        trueAverage += workingAverage /n;
        workingAverage = 0;
    }

    return trueAverage * (x2-x1) * (y2-y1);
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
void fillHaltons(int n)
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

void freeHaltons()
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
* Time Complexity: O(n) after first call with a consistent sample size, O(n^2 * log(n)) on first call with that sample size
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
        fillHaltons(n);


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
        freeHaltons();

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

int main()
{
    srand((unsigned)time(NULL));

    struct SimulatorMathRect rect = {.x1 = -100, .x2 = 100, .y1 = -100, .y2 = 100};

    double output, time;

    
    int samples[] = {256, 512, 1024, 2048};
    int a[] = {10, 100, 1000};
    int b[] = {15, 150, 1500};

    char str[256]; 

    FILE* file = fopen("QuasiMonteCarloData.csv", "w");

    if (file == NULL)
        return -1;

    fputs("func,samples,boundsType,a,b,time,result\n", file);

    for (int i = 0; i < 4; i++)
    {
        fillHaltons(samples[i]);
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                struct SimulatorMathRect b1 = {.x1 = 0, .x2 = a[j], .y1 = 0, .y2 = a[j]}; // 0,a, 0,a
                struct SimulatorMathRect b2 = {.x1 = 0, .x2 = a[j], .y1 = 0, .y2 = b[k]}; // 0,a 0,b
                struct SimulatorMathRect b3 = {.x1 = -1 * a[j], .x2 = a[j], .y1 = 0, .y2 = a[j]}; // -a,a 0, a
                struct SimulatorMathRect b4 = {.x1 = -1 * a[j], .x2 = a[j], .y1 = 0, .y2 = b[k]}; // -a,a 0,b
                struct SimulatorMathRect b5 = {.x1 = -1 * a[j], .x2 = a[j], .y1 = -1 * a[j], .y2 = a[j]}; // -a,-a a,a
                struct SimulatorMathRect b6 = {.x1 = -1 * a[j], .x2 = a[j], .y1 = -1 * b[k], .y2 = b[k]}; // -a,a -b,b

                struct SimulatorMathRect boundsArr[] = {b1, b2, b3, b4, b5, b6};

                for (int l = 0; l < 6; l++)
                {
                    time = _timeFunc(QuasiMonteCarloIntegral2D, f1, &boundsArr[l], samples[i], &output);
                    sprintf(str, "%d,%d,%d,%d,%d,%f,%f\n",1,samples[i],l+1,a[j],b[k],time,output);
                    fputs(str,file);

                    time = _timeFunc(QuasiMonteCarloIntegral2D, f2, &boundsArr[l], samples[i], &output);
                    sprintf(str, "%d,%d,%d,%d,%d,%f,%f\n",2,samples[i],l+1,a[j],b[k],time,output);
                    fputs(str,file);
                }

            }
        }
        freeHaltons();
    }

    return fclose(file);
    
    
    /*
   _logDouble("\nMidpoint Time",_timeFunc(MidpointSumIntegral2D, f1, &rect, 1000, &output));
   _logDouble("Result", output);

   _logDouble("\nTrapezoidal Time",_timeFunc(TrapezoidalSumIntegral2D, f1, &rect, 1000, &output));
   _logDouble("Result", output);

   _logDouble("\nSimpsons Time",_timeFunc(SimpsonsIntegral2D, f1, &rect, 1000, &output));
   _logDouble("Result", output);

   _logDouble("\nMonte Carlo Time", _timeFunc(MonteCarloIntegral2D, f1, &rect, 1000 * 1000, &output));
   _logDouble("Result", output);

    fillHaltons(1000*1000);
   _logDouble("\nQuasi Monte Carlo Time", _timeFunc(QuasiMonteCarloIntegral2D, f1, &rect,1000*1000, &output));
    freeHaltons();
   _logDouble("Result", output);
   */

    return 0;
    
}