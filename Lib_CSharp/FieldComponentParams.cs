using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.ComponentModel;

namespace Lib_CSharp
{
    [Serializable]
    public enum FldComp { Ez_vd, Ez_hd }
    [Serializable]
    public class FieldComponentParams: IDataErrorInfo
    {
        public FldComp FC { get; set; }
        public double roMin { get; set; }// в длинах волн в воздухе
        public double roMax { get; set; }// в длинах волн в воздухе
        public int nro { get; set; }
        public string Error => "Error";

        public string this[string columnName]
        {
            get
            {
                string str = "";
                switch (columnName)
                {
                    case "nro":
                        if (nro < 3) str = "nro must be greater 2";
                        break;
                    case "roMin":
                    case "roMax":
                        if (roMin <= 0 || roMin > roMax) str = "roMin must be positive; roMin must be less roMax";
                        break;
                }
                return str;
            }
        }

        public FieldComponentParams(int nro = 10, double roMin = 0.01, double roMax = 5, FldComp FC = FldComp.Ez_vd)
        {
            this.nro = nro;
            this.roMin = roMin;
            this.roMax = roMax;
            this.FC = FC;
        }

        public override string ToString()
        {
            return $"FC = {FC}  nro = {nro}" +
                   $"\nroMin = {roMin} roMax = {roMax}";
        }
    }
}
