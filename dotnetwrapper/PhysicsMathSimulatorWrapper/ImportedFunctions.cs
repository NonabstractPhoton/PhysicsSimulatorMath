using System.Runtime.InteropServices;

namespace PhysicsSimulatorMathWrapper
{
    public static class ImportedFunctions
    {
        /// <summary>
        /// A struct containing the rectangular bounds of a 2D integration region. 
        /// This rectangle has corners (x1,y1), (x1,y2), (x2,y1), and (x2,y2). 
        /// <list type="bullet">
        /// <item>
        ///     x1: Smaller 1st dimension bound
        /// </item>
        /// <item>
        ///     x2: Larger 1st dimension bound
        /// </item>
        /// <item>
        ///     y1: Smaller 2nd dimension bound
        /// </item>
        /// <item>
        ///     y2: Larger 2nd dimension bound 
        /// </item>
        /// </list>
        /// </summary>
        [StructLayout(LayoutKind.Sequential)]
        public struct SimulatorMathRect
        {
            public double x1, x2, y1, y2;
        };

        /// <summary>
        /// Represents a mathemtatical function of two variables
        /// </summary>
        public delegate double Function2D(double x, double y);


        /// <summary>
        /// Performs 2 dimensional numerical integration on the given function using a Midpoint Riemann Sum with n rectangles of equal area spanning the integration bounds. 
        /// Time Complexity: O(n^2)
        /// </summary>
        /// <param name="func">
        ///     The function of 2 variables to be integrated
        /// </param>
        /// <param name="bounds">
        ///     The bounds of the integration region in the form of a pointer to a SimulatorMathRect
        /// </param>
        /// <param name="n">
        ///     The amount of rectangles to be used
        /// </param>
        /// <returns>
        ///      The approximated value of the integral
        /// </returns>
        [DllImport("PhysicsSimulatorMath.dll", ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double MidpointSumIntegral2D(Function2D func, in SimulatorMathRect bounds, int n);


        /// <summary>
        /// Performs 2 dimensional numerical integration on the given function using an extension of the Trapezoidal Rule with n shapes spanning the integration bounds.
        /// Time Complexity: O(n^2)
        /// </summary>
        /// <param name="func">
        ///     The function of 2 variables to be integrated
        /// </param>
        /// <param name="bounds">
        ///     The bounds of the integration region in the form of a pointer to a SimulatorMathRect
        /// </param>
        /// <param name="n">
        ///     The amount of rectangles to be used
        /// </param>
        /// <returns>
        ///      The approximated value of the integral
        /// </returns>
        [DllImport("PhysicsSimulatorMath.dll", ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double TrapezoidalSumIntegral2D(Function2D func, in SimulatorMathRect bounds, int n);


        /// <summary>
        /// Performs 2 dimensional numerical integration on the given function using an extension of Simpson's 3/8 Rule with n shapes spanning the integration bounds.
        /// Time Complexity: O(n^2)
        /// </summary>
        /// <param name="func">
        ///     The function of 2 variables to be integrated
        /// </param>
        /// <param name="bounds">
        ///     The bounds of the integration region in the form of a pointer to a SimulatorMathRect
        /// </param>
        /// <param name="n">
        ///     The amount of rectangles to be used
        /// </param>
        /// <returns>
        ///      The approximated value of the integral
        /// </returns>
        [DllImport("PhysicsSimulatorMath.dll", ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double SimpsonsIntegral2D(Function2D func, in SimulatorMathRect bounds, int n);


        /// <summary>
        /// Performs 2 dimensional numerical integration on the given function using the Monte Carlo Method with a modified gaussian probability distribution.
        /// Uses the Central Limit Theorem with n samples of size 50 to approximate true mean result.
        /// Time Complexity: O(n) 
        /// </summary>
        /// <param name="func">
        ///     The function of 2 variables to be integrated
        /// </param>
        /// <param name="bounds">
        ///     The bounds of the integration region in the form of a pointer to a SimulatorMathRect
        /// </param>
        /// <param name="n">
        ///     The amount of rectangles to be used
        /// </param>
        /// <returns>
        ///      The approximated value of the integral
        /// </returns>
        [DllImport("PhysicsSimulatorMath.dll", ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double GaussianMonteCarloIntegral2D(Function2D func, in SimulatorMathRect bounds, int n);


        /// <summary>
        /// Performs 2 dimensional numerical integration on the given function using the Monte Carlo Method with a pseudo-random probability distribution.
        /// Uses the Central Limit Theorem with n samples of size 100 to approximate true mean result.
        /// Time Complexity: O(n) 
        /// </summary>
        /// <param name="func">
        ///     The function of 2 variables to be integrated
        /// </param>
        /// <param name="bounds">
        ///     The bounds of the integration region in the form of a pointer to a SimulatorMathRect
        /// </param>
        /// <param name="n">
        ///     The amount of rectangles to be used
        /// </param>
        /// <returns>
        ///      The approximated value of the integral
        /// </returns>
        [DllImport("PhysicsSimulatorMath.dll", ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double PsuedoRandMonteCarloIntegral2D(Function2D func, in SimulatorMathRect bounds, int n);


        /// <summary>
        /// Pre-computes halton sequence values for bases 2 and 3 for use in the Quasi Monte Carlo Integral.
        /// Must be called before use the Quasi Monte Carlo Integral for max effeciency on repeated use.
        /// Must be freed with FreeHaltons after use.
        /// </summary>
        /// <param name="n">Sequence Length</param>
        [DllImport("PhysicsSimulatorMath.dll", ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern void FillHaltons(int n);

        /// <summary>
        /// Frees the halton sequences generated by FillHaltons
        /// </summary>
        [DllImport("PhysicsSimulatorMath.dll", ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern void FreeHaltons();


        /// <summary>
        /// Performs 2 dimensional numerical integration on the given function using a Quasi Monte Carlo method 
        /// selecting n values from the Hatlon sequence. Must ensure that the size of the Halton sequence arrays matches n if pre-filled.
        /// Uses the Central Limit Theorem with samples of size 50, thus also handling large values of n.
        /// It is highly recommended to use fillHaltons before the execution of this algorithm and freeHaltons after for max effeciency.
        /// Time Complexity: O(n) if fillHaltons is properly called, O(n^2 log(n)) otherwise
        /// </summary>
        /// <param name="func">
        ///     The function of 2 variables to be integrated
        /// </param>
        /// <param name="bounds">
        ///     The bounds of the integration region in the form of a pointer to a SimulatorMathRect
        /// </param>
        /// <param name="n">
        ///     The amount of rectangles to be used
        /// </param>
        /// <returns>
        ///      The approximated value of the integral
        /// </returns>
        [DllImport("PhysicsSimulatorMath.dll", ExactSpelling = true, CallingConvention = CallingConvention.Cdecl)]
        public static extern double QuasiMonteCarloIntegral2D(Function2D func, in SimulatorMathRect bounds, int n);

    }
}
