using ScottPlot;
using ScottPlot.TickGenerators;
using System.Diagnostics.CodeAnalysis;
using static System.Runtime.InteropServices.JavaScript.JSType;

namespace DataVisualizer
{
    internal class Program
    {
        static void Main(string[] args)
        {
            string pathToData = Path.GetFullPath("../../../../../data/");

            string[] approxFileNames = 
            { 
                "GaussianMonteCarloData.csv", 
                "PsuedoRandMonteCarloData.csv", 
                "QuasiMonteCarloData.csv", 
                "MidpointSumData.csv", 
                "SimpsonsRuleData.csv", 
                "TrapezoidalSumData.csv" 
            };

            string valueFileName = "TrueValues.csv";
            
            List<TrueValueData> trueValues = new List<TrueValueData>(72);

            using (StreamReader reader = File.OpenText(pathToData+valueFileName))
            {
                string? line;

                reader.ReadLine();

                while ((line = reader.ReadLine()) != null)
                {
                    trueValues.Add(new TrueValueData(line));
                }
            }

            List<List<ApproximationData>> approxDatas = new List<List<ApproximationData>>(6);

            int index = 0;

            foreach (string fileName in approxFileNames)
            {
                approxDatas.Add(new List<ApproximationData>(432));

                using (StreamReader reader = File.OpenText(pathToData + fileName))
                {
                    reader.ReadLine();

                    string? line;

                    while ((line = reader.ReadLine()) != null)
                    {
                        approxDatas[index].Add(new ApproximationData(line));
                    }
                }

                index++;
            }

            for (int i = 0; i < approxDatas.Count; i++)
            {
                for (int j = 0; j < approxDatas[i].Count; j++)
                {
                    approxDatas[i][j].error = approxDatas[i][j].result - trueValues.Where
                    (
                        (data) => data.integralData == approxDatas[i][j].integralData
                    )
                    .First().result;
                }
            }

            // Types
            for (int i = 0; i < approxFileNames.Length; i++)
            {
                Plot pTimes = new();
                Plot pErrors = new();

                List<List<int>> xsList = new();
                List<List<double>> ysTimesList = new();
                List<List<double>> ysErrorsList = new();

                var func1Set = approxDatas[i].Where((data) => data.integralData.func == 1).ToList();
                var func2Set = approxDatas[i].Where((data) => data.integralData.func == 2).ToList();

                // Bounds
                for (int j = 0; j < 6; j++)
                {
                    var func1Subset = func1Set.Where((data) => data.integralData.boundsType == j + 1).ToList();
                    var func2Subset = func2Set.Where((data) => data.integralData.boundsType == j + 1).ToList();

                    for (int alog = 4; alog > 0; alog/=2)
                    {
                        for (int blog = 64; blog >= 8; blog /= 2)
                        {
                            IntegralData.ABMarker marker = (IntegralData.ABMarker)alog | (IntegralData.ABMarker)blog;

                            var lineSet = func1Subset.Where((data) => data.integralData.abMarker == marker);

                            xsList.Add(lineSet.Select((data) => data.samples).ToList());
                            ysTimesList.Add(lineSet.Select((data) => data.time).ToList());
                            ysErrorsList.Add(lineSet.Select((data) => data.error).ToList());
                        }
                    }

                }

                var timeTitle = approxFileNames[i][0..^4] + " Times";
                var errorTitle = approxFileNames[i][0..^4] + " Errors";

                pTimes.Title(timeTitle);
                pErrors.Title(errorTitle);
                
                pTimes.XLabel("Samples");
                pErrors.XLabel("Samples");

                pTimes.YLabel("Time (Milliseconds)");
                pErrors.YLabel("Error");

                int listIndex = 0;
                foreach (var xList in xsList)
                {
                    pTimes.Add.Scatter(xList, ysTimesList[listIndex]);
                    pErrors.Add.Scatter(xList, ysErrorsList[listIndex]);
                    listIndex++;
                }

                pTimes.SavePng(timeTitle+ ".png", 800, 450);
                pErrors.SavePng(errorTitle+ ".png", 800, 450);

            }
            // TODO - Add Regression
        }
    }

    struct IntegralData
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

        public ApproximationData(string rawData)
        {
            var data = rawData.Split(",");
            short b = 0;

            if (short.Parse(data[2]) % 2 == 0)
                 b = short.Parse(data[4]);

            integralData = new IntegralData(short.Parse(data[0]), short.Parse(data[2]), short.Parse(data[3]), b);

            samples = int.Parse(data[1]);
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
}
