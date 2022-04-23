using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.ComponentModel;


namespace Lib_CSharp
{
    [Serializable]
    public enum DpF { W, Wzz, U, Uz, V, Vzz }

    [Serializable]
    public class LambdaParams : IDataErrorInfo, INotifyPropertyChanged
    {
        public int n { get; set; }
        double _xL; 
        double _xR;
        public double lamL 
        { get { return _xL; }
          set 
            { _xL = value; 
                if(PropertyChanged != null)
                {
                    PropertyChanged(this, new PropertyChangedEventArgs("lamL"));
                    PropertyChanged(this, new PropertyChangedEventArgs("lamR"));
                }
            }
        }
        public double lamR
        {
            get { return _xR; }
            set
            {  _xR = value;
                if (PropertyChanged != null)
                {
                    PropertyChanged(this, new PropertyChangedEventArgs("lamL"));
                    PropertyChanged(this, new PropertyChangedEventArgs("lamR"));
                }
            }
        }
        public DpF F { get; set; }
        public bool lamIsUniform { get; set; }
        public double[] lam { get; set; }
        //public double[] mData { get; set; }
        //public bool IsLambdaDataCreated { get; set; }

        public LambdaParams(int n = 101, double lamL = 0, double lamR = 1, bool lamIsUniform = true, DpF F = DpF.W)
        {
            this.F = F;
            this.n = n;
            this.lamL = lamL;
            this.lamR = lamR;
            this.lamIsUniform = lamIsUniform;
            Create();
        }

        public event PropertyChangedEventHandler PropertyChanged;

        public void Create()
        {
            try
            {
                if (lamIsUniform) createXUniform();
                else createXLog();

                //this.mData = new double[n];
  
                    //Func<double, double> f;
                    //if (F == DpF.W) f = W;
                    //else if (F == DpF.U) f = U;
                    //else f = G; 

                    //for (int j = 0; j < n; j++)
                    //{
                    //    mData[j] = f(lam[j]);
                    //}

                //IsLambdaDataCreated = true;
            }
            catch (Exception ex)
            {
                throw new Exception("Error in LambdaData.Create()", ex);
            }
        }
        public override string ToString()
        {
            return $"F = {F} n = {n} lamL = {lamL} lamR = {lamR}";
        }
        public string Error => "Error";

        public string this[string columnName]
        {
            get
            {
                string str = "";
                switch (columnName)
                {
                    case "n":
                        if (n < 3) str = "n must be greater 2";
                        break;
                    case "lamR":
                    case "lamL":
                        if (lamL >= lamR) str = "lamL must be less then lamR";
                        break;
                }
                return str;
            }
        }
        //double W(double x)
        //{ return x * 0.01 + 2; }
        //double U(double x)
        //{ return x * (x - 1) + 2; }
        //double G(double x)
        //{ return x * (x - 1) * (x - 2) + 2; }
        void createXLog()
        {
            this.lam = new double[n];
            //Random rnd = new Random();
            //for (int j = 0; j < n; j++)
            //{
            //    x[j] = rnd.NextDouble() * (lamR - lamL) + lamL;
            //}
            //Array.Sort(x);
            //x[0] = lamL;
            //x[n - 1] = lamR;
        }
        void createXUniform()
        {
            this.lam = new double[n];
            double h = (lamR - lamL) / (n - 1);
            for (int j = 0; j < n; j++)
            {
                lam[j] = lamL + j * h;
            }
        }
    }
}
