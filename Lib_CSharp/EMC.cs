using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Lib_CSharp
{
    public static class EMC
    {
        public static double mu = 4 * Math.PI * 1.0E-7;
        public static double eps0 = 8.8541878e-12;
        public static double lam0om1 = 2 * Math.PI / Math.Sqrt(eps0 * mu);
    }
}
