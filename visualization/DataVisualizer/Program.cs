using ScottPlot;
using System.Data;
using static alglib;

namespace DataVisualizer
{
    internal class Program
    {
        static void Main(string[] args)
        {
            string pathToData = Path.GetFullPath("../../../../../data/");
            // GenerateDetailedGraphs(pathToData);
            GenerateNormalGraphs(pathToData);

        }

        static void GenerateDetailedGraphs(string pathToData)
        {
            string[] approxFileNames =
            {
                "MidpointSumDataDetailed.csv",
                "SimpsonsRuleDataDetailed.csv",
                "TrapezoidalSumDataDetailed.csv"
            };

            string valueFileName = "TrueValues.csv";

            List<TrueValueData> trueValues = loadTrueData(pathToData, valueFileName, 72);

            List<List<ApproximationData>> approxDatas = loadApproxData(pathToData, approxFileNames, 396);

            for (int i = 0; i < approxFileNames.Length; i++)
            {
                foreach (ApproximationData dataObj in approxDatas[i])
                {
                    var trueVal = trueValues.Where((data) => data.integralData == dataObj.integralData).First().result;
                    dataObj.error = (dataObj.result - trueVal) / trueVal * 100;
                }
            }


            for (int i = 0;i < approxFileNames.Length;i++) 
            {
                Plot pErrors = new(), pTimes = new();

                var timeTitle = $"{approxFileNames[i][0..^4]} Times";
                var errorTitle = $"{approxFileNames[i][0..^4]} Errors";

                pTimes.Title(timeTitle);
                pErrors.Title(errorTitle);

                pTimes.XLabel("Samples");
                pErrors.XLabel("Samples");

                pTimes.YLabel("Time (milliseconds)");
                pErrors.YLabel("Error (%)");

                for (int boundsType = 1; boundsType <= 6; boundsType++)
                {
                    for (int func = 1; func <=2; func++)
                    {
                        var lineSet = approxDatas[i].Where((data) => data.integralData.boundsType == boundsType && data.integralData.func == func);
                        var xs = lineSet.Select(data => data.samples).ToList();
                        var yTimes = lineSet.Select(data => data.time).ToList();
                        var yErrors = lineSet.Select(data => data.error).ToList();

                        var timePlot = pTimes.Add.Scatter(xs, yTimes);
                        var errorPlot = pErrors.Add.Scatter(xs, yErrors);

                        timePlot.LegendText = $"Func: {func}, Bounds: {boundsType}";
                        errorPlot.LegendText = $"Func: {func}, Bounds: {boundsType}";
                    }
                }

                pTimes.SavePng($"{timeTitle}.png", 1920, 1080);
                pErrors.SavePng($"{errorTitle}.png", 1920, 1080);


                Plot pAvgErrors = new(), pAvgTimes = new();

                var avgTimeTitle = $"{approxFileNames[i][0..^4]} Average Times Regression";
                var avgErrorTitle = $"{approxFileNames[i][0..^4]} Average Errors Regression";

                pAvgTimes.Title(avgTimeTitle);
                pAvgErrors.Title(avgErrorTitle);

                pAvgTimes.XLabel("Samples");
                pAvgErrors.XLabel("Samples");

                pAvgTimes.YLabel("Average Time (milliseconds)");
                pAvgErrors.YLabel("Average Error (%)");

                List<double> errors = new();

                for (int func = 1; func <= 2; func++)
                {
                    List<double> avgErrors = new(), avgTimes = new(), samples = new();

                    var trials = approxDatas[i].Where((data) => data.integralData.func == func).OrderBy((data) => data.samples).ToList();

                    // Average Values Across Bounds

                    while (trials.Count() > 0)
                    {
                        var current = trials.Take(6);
                        avgTimes.Add(current.Average((data) => data.time));
                        avgErrors.Add(current.Average((data) => data.error));
                        samples.Add(current.First().samples);
                        trials = trials.Skip(6).ToList();
                    }

                    // Plot Avgs Scatter

                    var errPlot = pAvgErrors.Add.Scatter(samples, avgErrors);
                    var timePlot = pAvgTimes.Add.Scatter(samples, avgTimes);

                    int info = 0;

                    ratint.barycentricinterpolant timeInterpolant = new(), errorInterpolant = new();
                    lsfit.polynomialfitreport report = new();
                    xparams _params = new(0);

                    lsfit.polynomialfit(samples.ToArray(), avgTimes.ToArray(), avgTimes.Count(), 3, ref info, timeInterpolant, report, _params);

                    double[] timeCoeffs = new double[3];

                    polint.polynomialbar2pow(timeInterpolant, 0, 1, ref timeCoeffs, _params);

                    double[] errorCoeffs = new double[3];

                    lsfit.polynomialfit(samples.ToArray(), avgErrors.ToArray(), avgErrors.Count(), 2, ref info, errorInterpolant, report, _params);

                    polint.polynomialbar2pow(errorInterpolant, 0, 1, ref errorCoeffs, _params);


                    // Plot Regression Lines

                    errPlot.LegendText = $"Func {func} Averages";
                    timePlot.LegendText = $"Func {func} Averages";

                    var curve = pAvgTimes.Add.Function((x) => timeCoeffs[0] + timeCoeffs[1] * x + timeCoeffs[2] * x * x);

                    curve.LegendText = $"Func {func}, {timeCoeffs[0]} + {timeCoeffs[1]}x + {timeCoeffs[2]}x^2";
                    curve.LinePattern = LinePattern.Dashed;


                    var errLimits = errPlot.GetAxisLimits();

                    curve = pAvgErrors.Add.Function((x) => errorCoeffs[0] + errorCoeffs[1] * x);
                    
                    curve.LinePattern = LinePattern.Dashed;
                    curve.LegendText = $"Func {func}, {errorCoeffs[0]} + {errorCoeffs[1]}x";

                    errors.AddRange(avgErrors);
                }

                // Auto Scale Correction
                double minError = errors.Min(), maxError = errors.Max();
                pAvgErrors.Axes.SetLimitsY(minError < 0 ? minError * 1.25 : minError / 1.25, maxError < 0 ? maxError / 1.25 : maxError * 1.25);

                pAvgErrors.SavePng($"{avgErrorTitle}.png", 1920, 1080);
                pAvgTimes.SavePng($"{avgTimeTitle}.png", 1920, 1080);
            }
        }

        static void GenerateNormalGraphs(string pathToData)
        {
            string[] approxFileNames =
            {
                "MidpointSumData.csv",
                "SimpsonsRuleData.csv",
                "TrapezoidalSumData.csv",
                "GaussianMonteCarloData.csv",
                "PsuedoRandMonteCarloData.csv",
                "QuasiMonteCarloData.csv"
            };

            string valueFileName = "TrueValues.csv";

            List<TrueValueData> trueValues = loadTrueData(pathToData, valueFileName, 72);

            List<List<ApproximationData>> approxDatas = loadApproxData(pathToData, approxFileNames, 288);



            for (int i = 0; i < approxDatas.Count; i++)
            {
                for (int j = 0; j < approxDatas[i].Count; j++)
                {
                    var trueVal = trueValues.Where
                    (
                        (data) => data.integralData == approxDatas[i][j].integralData
                    )
                    .First().result;


                    approxDatas[i][j].error = (approxDatas[i][j].result - trueVal) / trueVal * 100;
                }
            }

            List<PerformanceSummary> summary = new(approxFileNames.Length);

            // Types
            for (int i = 0; i < approxFileNames.Length; i++)
            {
                Plot pTimes = new();
                Plot pErrors = new();

                List<List<ApproximationData>> dataList = new();

                var func1Set = approxDatas[i].Where((data) => data.integralData.func == 1).ToList();
                var func2Set = approxDatas[i].Where((data) => data.integralData.func == 2).ToList();

                // Bounds
                for (int j = 1; j <= 6; j++)
                {
                    var func1Subset = func1Set.Where((data) => data.integralData.boundsType == j).ToList();
                    var func2Subset = func2Set.Where((data) => data.integralData.boundsType == j).ToList();

                    for (int alog = 4; alog > 0; alog /= 2)
                    {
                        for (int blog = 64; blog >= 8; blog /= 2)
                        {
                            IntegralData.ABMarker marker = (IntegralData.ABMarker)alog | (IntegralData.ABMarker)blog;

                            var lineSet1 = func1Subset.Where((data) => data.integralData.abMarker == marker);
                            var lineSet2 = func2Subset.Where((data) => data.integralData.abMarker == marker);

                            if (lineSet1.Count() == 0)
                            {
                                continue;
                            }

                            dataList.Add(lineSet1.ToList());
                            dataList.Add(lineSet2.ToList());
                        }
                    }

                }

                summary.Add(new PerformanceSummary() { typeIndex = i, data = new() });

                for (int samples = 256; samples <= 2048; samples *= 2)
                {
                    var dataOfSamples = approxDatas[i].Where((data) => data.samples == samples);

                    summary[i].data.Add
                    ((
                        samples,
                        dataOfSamples.Select((data) => data.time).Average(),
                        dataOfSamples.Select((data) => data.error).Average()
                    ));
                }

                var sortedByTimes = dataList.OrderBy((data) => data.Select((innerData) => innerData.time).Average()).ToList();
                var sortedByErrors = dataList.OrderBy((data) => data.Select((innerData) => innerData.error).Average()).ToList();

                var timeTitle = approxFileNames[i][0..^4] + " Times";
                var errorTitle = approxFileNames[i][0..^4] + " Errors";

                pTimes.Title(timeTitle);
                pErrors.Title(errorTitle);

                pTimes.XLabel("Samples");
                pErrors.XLabel("Samples");

                pTimes.YLabel("Time (milliseconds)");
                pErrors.YLabel("Error (%)");

                for (int index = 0; index < dataList.Count() / 2; index++)
                {
                    var timesSamplesLower = sortedByTimes[index].Select((data) => data.samples).ToList();
                    var timesSamplesUpper = sortedByTimes[^(index + 1)].Select((data) => data.samples).ToList();

                    var timesLower = sortedByTimes[index].Select((data) => data.time).ToList();
                    var timesUpper = sortedByTimes[^(index + 1)].Select((data) => data.time).ToList();

                    var errorsSamplesLower = sortedByErrors[index].Select((data) => data.samples).ToList();
                    var errorsSamplesUpper = sortedByErrors[^(index + 1)].Select((data) => data.samples).ToList();

                    var errorsLower = sortedByErrors[index].Select((data) => data.error).ToList();
                    var errorsUpper = sortedByErrors[^(index + 1)].Select((data) => data.error).ToList();

                    var timePlotLowerData = sortedByTimes[index].First().integralData;
                    var timePlotUpperData = sortedByTimes[^(index + 1)].First().integralData;

                    var errorPlotLowerData = sortedByErrors[index].First().integralData;
                    var errorPlotUpperData = sortedByErrors[^(index + 1)].First().integralData;


                    var timePlotLower = pTimes.Add.Scatter(timesSamplesLower, timesLower);
                    var timePlotUpper = pTimes.Add.Scatter(timesSamplesUpper, timesLower);

                    var errorPlotLower = pErrors.Add.Scatter(errorsSamplesLower, errorsLower);
                    var errorPlotUpper = pErrors.Add.Scatter(errorsSamplesUpper, errorsLower);

                    if (index < 5)
                    {
                        timePlotLower.LegendText = $"Func: {timePlotLowerData.func}, Bounds: Type {timePlotLowerData.boundsType}, {timePlotLowerData.abMarker}";
                        timePlotUpper.LegendText = $"Func: {timePlotUpperData.func}, Bounds: Type {timePlotUpperData.boundsType}, {timePlotUpperData.abMarker}";

                        errorPlotLower.LegendText = $"Func: {errorPlotLowerData.func}, Bounds: Type {errorPlotLowerData.boundsType}, {errorPlotLowerData.abMarker}";
                        errorPlotUpper.LegendText = $"Func: {errorPlotUpperData.func}, Bounds: Type {errorPlotUpperData.boundsType}, {errorPlotUpperData.abMarker}";
                    }
                }

                pTimes.Legend.IsVisible = false;
                pErrors.Legend.IsVisible = false;

                ScottPlot.Panels.LegendPanel timePanel = new(pTimes.Legend)
                {
                    Edge = Edge.Right,
                    Alignment = Alignment.UpperCenter
                };
                ScottPlot.Panels.LegendPanel errorPanel = new(pErrors.Legend)
                {
                    Edge = Edge.Right,
                    Alignment = Alignment.UpperCenter
                };

                pTimes.Axes.AddPanel(timePanel);
                pErrors.Axes.AddPanel(errorPanel);

                pTimes.SavePng(timeTitle + ".png", 1920, 1080);
                pErrors.SavePng(errorTitle + ".png", 1920, 1080);

            }
            // TODO - Add Regression

            Plot avgErrorsComparisons = new();
            Plot avgTimesComparisons = new();

            for (int i = 0; i < approxFileNames.Length; i++)
            {
                var xs = summary[i].data.Select((triple) => (double)triple.Item1).ToList();

                var yTimes = summary[i].data.Select((triple) => triple.Item2).ToList();
                var timeInfo = avgTimesComparisons.Add.Scatter(xs, yTimes);
                timeInfo.LegendText = approxFileNames[i][0..^4];


                if (i != 3)
                {
                    var yErrors = summary[i].data.Select((triple) => triple.Item3).ToList();
                    var errorInfo = avgErrorsComparisons.Add.Scatter(xs, yErrors);
                    errorInfo.LegendText = approxFileNames[i][0..^4];
                }
            }

            avgTimesComparisons.Title("Average Times");
            avgErrorsComparisons.Title("Average Errors");

            avgTimesComparisons.XLabel("Samples");
            avgErrorsComparisons.XLabel("Samples");

            avgTimesComparisons.YLabel("Time (milliseconds)");
            avgErrorsComparisons.YLabel("Error (%)");

            avgTimesComparisons.SavePng("Average Times.png", 1920, 1080);
            avgErrorsComparisons.SavePng("Average Errors.png", 1920, 1080);
        }

        static List<TrueValueData> loadTrueData(string pathToData, string valueFileName, int lines)
        {
            List<TrueValueData> trueValues = new List<TrueValueData>(72);

            using (StreamReader reader = File.OpenText(pathToData + valueFileName))
            {
                string? line;

                reader.ReadLine();

                while ((line = reader.ReadLine()) != null)
                {
                    trueValues.Add(new TrueValueData(line));
                }
            }

            return trueValues;
        }

        static List<List<ApproximationData>> loadApproxData(string pathToData, string[] fileNames, int lines)
        {
            List<List<ApproximationData>> approxDatas = new List<List<ApproximationData>>(fileNames.Length);

            for (int i = 0; i < fileNames.Length; i++)
            {
                approxDatas.Add(new List<ApproximationData>(lines));

                using (StreamReader reader = File.OpenText(pathToData + fileNames[i]))
                {
                    reader.ReadLine();

                    string? line;

                    while ((line = reader.ReadLine()) != null)
                    {
                        // Quasi Monte Carlo
                        if (i == 5)
                            approxDatas[i].Add(new ApproximationData(line, true));
                        else
                            approxDatas[i].Add(new ApproximationData(line));
                    }
                }
            }

            return approxDatas;
        }
    }

    class IntegralData
    {
        public short func, boundsType, a, b;

        [Flags]
        public enum ABMarker
        {
            a10 = 1,
            a100 = 2,
            a1000 = 4,
            b15 = 8,
            b150 = 16,
            b1500 = 32,
            bx = 64
        }


        public ABMarker abMarker;

        public IntegralData(short func, short boundsType, short a, short b)
        {
            this.func = func;
            this.boundsType = boundsType;
            this.a = a;
            this.b = b;

            ABMarker aMarker = (ABMarker)Math.Pow(2, (int)Math.Log10(a)-1);
            ABMarker bMarker;

            if (b == 0)
                bMarker = ABMarker.bx; 
            else
                bMarker = (ABMarker)Math.Pow(2, Math.Round((Math.Log(b) / Math.Log(15)) + 2));

            abMarker = aMarker | bMarker;
        }

        public static bool operator ==(IntegralData d1, IntegralData d2)
        {
            return d1.func == d2.func
                && d1.boundsType == d2.boundsType
                && d1.abMarker == d2.abMarker;
        }

        public static bool operator !=(IntegralData d1, IntegralData d2)
        {
            return !(d1 == d2);
        }

    }


    class ApproximationData
    {
        public int samples;
        public IntegralData integralData;
        public double time, result, error;

        public ApproximationData(string rawData, bool sqrtSamples = false)
        {
            var data = rawData.Split(",");

            short boundsType = short.Parse(data[2]);
            short b = 0;

            // If b matters
            if (boundsType % 2 == 0)
                b = short.Parse(data[4]);


            integralData = new IntegralData(short.Parse(data[0]), boundsType, short.Parse(data[3]), b);

            samples = sqrtSamples ? (int)Math.Sqrt(int.Parse(data[1])) : int.Parse(data[1]);
            time = double.Parse(data[5]);
            result = double.Parse(data[6]);
        }
        public override string ToString()
        {
            return $"func: {integralData.func}, samples: {samples}, boundsType: {integralData.boundsType}, a: {integralData.a}, b: {integralData.b}, time: {time}, result:  {result}, error: {error}";
        }

    }

    class TrueValueData
    {
        public IntegralData integralData;
        public double result;

        public TrueValueData(string rawData)
        {
            var data = rawData.Split(",");
            short b = 0;

            if (data[3].Length > 0)
                b = short.Parse(data[3]);

            integralData = new IntegralData(short.Parse(data[0]), short.Parse(data[1]), short.Parse(data[2]), b);

            result = double.Parse(data[4]);
        }

        public override string ToString()
        {
            return $"func: {integralData.func}, boundsType: {integralData.boundsType}, a: {integralData.a}, b: {integralData.b}, result:  {result}";
        }
    }

    struct PerformanceSummary
    {
        public int typeIndex;
        
        // Samples, time, error%
        public List <(int, double, double)> data;
    }

}
