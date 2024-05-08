#include "SimulatorMath.h"

int h = 10;

// Test Function
double f(double x, double y)
{
    return h/pow((x*x+h*h), 1.5);
}

double g(double x, double y)
{
    return 5*x+4*y;
}

double MidpointSumIntegral(double x1, double x2, double y1, double y2, int n)
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
            area += g(xStart + i*xStep, yStart + i*yStep);
        }
    }

    return area * xStep * yStep;
}

double TrapezoidalSumIntegral(double x1, double x2, double y1, double y2, int n)
{
    double xStep = (x2-x1)/n;
    double yStep = (y2-y1)/n;

    double area = 0;

    area += g(x1, y1);
    area += g(x2, y1);
    area += g(x1, y2);
    area += g(x2, y2);

    for (int i = 1; i < n-1; i++)
    {
        area += 2*g(x1 + i*xStep, y1);
        area += 2*g(x1 + i*xStep, y2);
        area += 2*g(x1, y1 + i*yStep);
        area += 2*g(x2, y1 + i*yStep);

        for (int j = 1; j < n-1; j++)
        {
            area += 4*g(x1 + i*xStep, y1 + i*yStep);
        }
    }

    return area * xStep * yStep / 4;      
}

int main()
{
    char str1[100], str2[100];
    printf("Running\n");
    sprintf(str1, "%f", MidpointSumIntegral(0,4,0,6,10));
    puts(str1);
    sprintf(str2, "%f", TrapezoidalSumIntegral(0,4,0,6,10));
    puts(str2);
    return 0;
}