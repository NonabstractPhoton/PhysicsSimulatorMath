#include "SimulatorMath.h"

int h = 10;
char buffer[100];


// Monte Carlo vars

int precision = 100;
double maxStdDev = 1.5;

// Quasi Monte Carlo Halton Sequence vars
int haltonSize = 100;
double halton2[100];
double halton3[100];



void logDouble(const char* label, double d)
{
    sprintf(buffer, "%s: %f", label, d);
    puts(buffer);
}

// Test Function
double f(double x, double y)
{
    return h/pow((x*x+h*h), 1.5);
}

double g(double x, double y)
{
    return 5*x+4*y;
}

double MidpointSumIntegral(double (*func)(double, double), double x1, double x2, double y1, double y2, int n)
{    
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

double TrapezoidalSumIntegral(double (*func)(double, double), double x1, double x2, double y1, double y2, int n)
{
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

double _genGaussRand(double lower, double upper, int precision, double maxStdDev)
{   
    // Ex precision of 100 gives rand [0, 100] gives [-50, 50] gives [-1.0, 1.0] with a step of .02 given by 1/(precision/2). This valu eis the amount of stdDevs we can go from the center  
    double stdDev = (rand() % (precision+1) - precision/2) / (double) (precision/2) * maxStdDev;
    
    double fittedValue = 1-exp(-1*stdDev*stdDev);

    double center = (upper+lower)/2;
    double range = (upper-lower);

    // center +- fittedValue * distance
    return center + ((stdDev > 0) - (stdDev < 0)) * fittedValue * range/2;


}

double MonteCarloIntegral(double (*func)(double, double), double x1, double x2, double y1, double y2, int n)
{
    double average = 0;

    for (int i = 0; i < n; i++)
    {
        average += func(_genGaussRand(x1, x2, precision, maxStdDev), _genGaussRand(y1, y2, precision, maxStdDev));
    }

    return average/n * (x2-x1) * (y2-y1);
}

double haltonSeq(int index, int base){
    
    double x = 1, y = 0;
    
    while (index > 0)
    {
        x /= base;
        y += x * (index % base);
        index /= base;
    }
    return y;
}

void fillHaltons()
{
    for (int i = 0; i < haltonSize; i++)
    {
        halton2[i] = haltonSeq(i+1, 2);
        halton3[i] = haltonSeq(i+1, 3);
    }
}

double QuasiMontiCarloIntegral(double (*func)(double, double), double x1, double x2, double y1, double y2, int n)
{
    double average = 0;
    double xRange = x2-x1, yRange = y2-y1;

    if (xRange > yRange)
    {
        for (int i = 0; i < n; i++)
        {
            average += func(x1+halton2[i]*xRange, y1+halton3[i]*yRange);
        }
    }
    else
    {
        for (int i = 0; i < n; i++)
        {
            average += func(x1+halton3[i]*xRange, y1+halton2[i]*yRange);
        }
    }

    return average/n * (xRange) * (yRange);
}

double SimpsonsIntegral(double (*func)(double, double), double x1, double x2, double y1, double y2, int n)
{
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

int main()
{
    printf("Running\n");
    srand((unsigned)time(NULL));

    logDouble("Midpoint", MidpointSumIntegral(g,0,4,0,6,10));

    logDouble("Trapezoidal", TrapezoidalSumIntegral(g,0,4,0,6,10));

    logDouble("Simpsons", SimpsonsIntegral(g,0,4,0,6,10));

    logDouble("Monte Carlo", MonteCarloIntegral(g,0,4,0,6,150));

    fillHaltons();

    logDouble("Quasi Monte Carlo", QuasiMontiCarloIntegral(g,0,4,0,6, haltonSize));

    return 0;
}