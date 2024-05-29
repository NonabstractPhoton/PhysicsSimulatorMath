# Numerical Integration for a Physics Simulator

## Methods
This project tests the Midpoint Rule, Trapezoidal Rule, Simpson's Rule, Monte Carlo Methods, and a Quasi Monte Carlo method for two-dimensional numerical integration.
It compares them based on runtime and error margins across a wide variety of integration regions and applies them to two different functions which will be integrated in the simulator.
All relevant data and a vast set of graphs are available within the data and graphs subdirectories of this repository. The C# code to perform regression analysis and graph the data is located within the
visualization subdirectory.

## Results
The results point to Simpson's Rule as the optimal integration method for our purposes due to its negligible error margins and only slighter higher execution time as compared to other deterministic methods.
A C library containing implementations of each tested numerical integration is available in .dll form within the dlls folder, with its header available in the src folder. All source code is available in the
src folder as well if one wishes to compile it themselves, potentially for a linux based operating system or to create a static object library. A fully documented .NET wrapper for this library is also 
available in .dll form within the dlls folder, with its C# source code available in the dotnetwrapper folder. Both versions of the library come with extensive documentation detailing use, implementation
methodologies, and time complexity. 

## Paper
A technical paper explaining this project in greater detail is attached to the repository. It includes a brief introduction to numerical integration, a deconstruction of its methodology, and a thorough
analysis of our results.
