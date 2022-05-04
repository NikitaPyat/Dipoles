using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.ComponentModel;

namespace Lib_CSharp
{
 
    [Serializable]
    // Хранятся ненормированные значения
    public class ModelParams: IDataErrorInfo, INotifyPropertyChanged
    {
        public event PropertyChangedEventHandler PropertyChanged;
        public double eps1 { get; set; } // множитель eps0
        public double eps2 { get; set; } // множитель eps0
        public double h { get; set; } // в длинах волн в слое
        public double ro { get; set; }// в длинах волн в воздухе
         public double z0 { get; set; } // в длинах волн в воздухе
        public double z { get; set; } // в длинах волн в воздухе
        public double wave { get; set; } // в длинах волн в воздухе
        public double _omega; // частота
       
        public double omega
        {
            get { return _omega; }
            set 
            {
                _omega = value;
                 wave = EMC.lam0om1 / _omega;
                if (PropertyChanged != null)
                    PropertyChanged(this, new PropertyChangedEventArgs("wave"));
            }
        }

        public string Error => "Error";
       
        public string this[string columnName]
        {
            get
            {
                string str = "";
                switch (columnName)
                {
                    case "ro":
                        if (ro <= 0) str = "ro must be positive";
                        break;
                    case "h":
                        if (h <= 0) str = "h must be positive";
                        break;
                    case "z0":
                    case "z":
                        if (z0 <= 0 || z < 0 || z > z0) str = "z0 must be greater 0; z in [0, z0]";
                        break;
                    case "eps1":
                        if (eps1 <= 1 ) str = "eps1 must be greater 1";
                        break;
                    case "eps2":
                        if (eps2 <= 1) str = "eps2 must be greater 1";
                        break;
                    case "wave":
                        if (wave <= 0) str = "w must be positive";
                        break;
                    case "omega":
                        if (_omega <= 0) str = "omega must be positive";
                        break;
                }
                return str;
            }
        }

        public ModelParams(double eps1 = 2, double eps2 = 4, double h = 2, 
                           double ro = 1, double z0 = 0.01, double z = 0, double omega = 100000)
        {
            this.eps1 = eps1;
            this.eps2 = eps2;
            this.h = h;
            this.ro = ro;
            this.z0 = z0;
            this.z = z;
            this._omega = omega;
            wave = EMC.lam0om1 / omega;
        }
        public double k0 
        {
            get { return 2 * Math.PI / wave;  }
        }
        public double k1
        {
            get { return 2 * Math.PI / wave * Math.Sqrt(eps1); }
        }
        public double k2
        {
            get { return 2 * Math.PI / wave * Math.Sqrt(eps2); }
        }
        public override string ToString()
        {
            return $"eps1 = {eps1} eps2 = {eps2} h = { h}" +
                $"\nro = {ro}" +
                $"\nz0 = {z0} z = {z}" + 
                $"\nomega = {omega} wave = {wave}";
        }
    }
}
