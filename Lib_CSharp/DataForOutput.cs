using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections.ObjectModel;

namespace Lib_CSharp
{
    public class DataForOutput
    {
        public double[] X { get; set; } // значения абсцисс 
        public List<double[]> YL { get; set; } // значения ординат 

        public ObservableCollection<string> infoList { get; set; }
        public DataForOutput()
        {
            infoList = new ObservableCollection<string>();
        }
        public DataForOutput(string Title, string xTitle, string yTitle, double[] X, List<double[]> YL, List<string> legends)
        {
            infoList = new ObservableCollection<string>();
            this.X = X;
            this.YL = YL;
        }
        public void DataToStringList(double[] X, double[] YL, double[] YL_CPP, string legend)
        {
                infoList.Add(legend);
                for (int j = 0; j < X.Length; j++)
                {
                    infoList.Add($"lam = {X[j].ToString("F7")}   y_CSharp = {YL[j].ToString("F7")}   y_CPP = {YL_CPP[j].ToString("F7")}");
                }
        }
    }
}

