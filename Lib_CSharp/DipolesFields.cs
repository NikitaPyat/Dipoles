using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;
using System.IO;
using System.Runtime.Serialization.Formatters.Binary;
using System.Numerics;
using Anderson;

namespace Lib_CSharp
{
   
    public class DipolesFields
    {
       public ModelParams modelParams { get; set; }
       public FieldComponentParams FieldComponent { get; set; }
        public List<string> TestString { get; set; }
       
        double eps1;
        double eps2;
        double h;
        double ro; 
        double z0;
        double z;
        double wave ;
        double k0;
        double k1;
        double k2;
        double k0k0;
        double k1k1;
        double k2k2;
        double k10;
        double k21;
      
        public DipolesFields( ModelParams modelParams, FieldComponentParams FieldComponent)
        {
            this.modelParams = new ModelParams(modelParams.eps1, modelParams.eps2, modelParams.h,
                                               modelParams.ro, modelParams.z0, modelParams.z, modelParams.wave);
            eps1 = modelParams.eps1;
            eps2 = modelParams.eps2;
            wave = EMC.lam0om1 / modelParams.omega;
            h = modelParams.h * wave;
            ro = modelParams.ro * wave;
            z0 = modelParams.z0 * wave;
            z = modelParams.z * wave;
          
            k0 = 2 * Math.PI / wave;
            k1 = k0 * Math.Sqrt(eps1);
            k2 = k0 * Math.Sqrt(eps2);
            k0k0 = k0 * k0;
           // k0k0 = 4 * Math.PI * Math.PI / wave / wave;
            k1k1 = k0k0 * eps1;
            k2k2 = k0k0 * eps2;
            k10 = eps1;  // eps1 / eps0
            k21 = eps2 / eps1;
            this.FieldComponent = FieldComponent;
            TestString = new List<string>();
        }
     
        // для отладки
        void CalculateFlam_CPP(double[] lam, List<double[]> res)
        {
            int nlam = lam.Length;
            try
            {
                double[] resCPP = new double[2 * nlam];
                Test(modelParams.eps1, modelParams.eps2, modelParams.h,
                     modelParams.ro, modelParams.z0, modelParams.z, modelParams.wave,
                     nlam, lam, resCPP);

                for (int j = 0; j < nlam; j++)
                {
                    res[0][j] = resCPP[2 * j];
                    res[1][j] = resCPP[2 * j + 1];
                }
            }
            catch (Exception ex)
            {
                throw new Exception("Error in CalculateFlam_CPP", ex);
            }
        }
        // для графики
        public void CalculateFlam_CSharp(double[] lam, List<double[]> res)
        {
            //TestString =
            //     $"\nMath.PI = {Math.PI}" +
            //     $"\nmu = {EMC.mu}" +
            //     $"\neps0 = {EMC.eps0}" +
            //     $"\nlam0om1  = {EMC.lam0om1}" +
            //     $"eps1 = {eps1} eps2 = {eps2} h = { h}" +
            //     $"\nro = {ro}" +
            //     $"\nz0 = {z0} z = {z}" +
            //     $"\n  wave = {wave}" +
            //     $"\n k0k0 = {k0k0}" + $"\n k1k1 = {k1k1}" + $"\n k2k2 = {k2k2}" +
            //     $"\n k10 = {k10}" + $"\n k21 = {k21}";

            int nlam = lam.Length;
            double[] fValues = new double[4];
            try
            {
                for (int j = 0; j < nlam; j++)
                {
                    Ez_hd_int_fun(lam[j], fValues);
                    res[0][j] = fValues[0];
                    res[1][j] = fValues[1];
                    res[2][j] = fValues[2];
                    res[3][j] = fValues[3];
                }
            }
            catch (Exception ex)
            {
                throw new Exception("Error in CalculateFlam_CSharp()", ex);
            }
        }

        public void Ez_vd_int_fun(double lambda, double[] fValues)
        {
            Complex tmp = new Complex(lambda * lambda - k0k0, 0);
            Complex nu0 = Complex.Sqrt(tmp);

            tmp = new Complex(lambda * lambda - k1k1, 0);
            Complex nu1 = Complex.Sqrt(tmp);
            Complex nu1k21 = k21 * nu1;

            tmp = new Complex(lambda * lambda - k2k2, 0);
            Complex nu2 = Complex.Sqrt(tmp);

            Complex e1 = Complex.Exp(-2 * nu1 * h);
            Complex F1 = (nu1k21 - nu2) / (nu1k21 + nu2) * e1;

            Complex G0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex G0N = k10 * (F1 + 1);
                Complex G0D = nu0 * G0N - nu1 * (F1 - 1);
                G0 = 2 * G0N / G0D;

            }
            else
            {
                Complex G0N = k10 * (k21 + h * nu2);
                Complex G0D = 2 * G0N * nu0 + nu2 * (1 - F1);
                G0 = 4 * G0N / G0D;
            }
            Complex e0 = Complex.Exp(-nu0 * (z + z0));
            double nu02 = lambda * lambda - k0k0;
            G0 = wave * G0 * e0 * (1 + nu02 / k0k0);   // k0k0??

            fValues[0] = G0.Real;        // Re
            fValues[1] = G0.Imaginary;   // Im
            fValues[2] = fValues[0];     // Dublicate
            fValues[3] = fValues[1];     // Dublicate
        }
        public void Ez_hd_int_fun(double lambda, double[] fValues)
        {
            Complex tmp = new Complex(lambda * lambda - k0k0, 0);
            Complex nu0 = Complex.Sqrt(tmp);

            tmp = new Complex(lambda * lambda - k1k1, 0);
            Complex nu1 = Complex.Sqrt(tmp);

            tmp = new Complex(lambda * lambda - k2k2, 0);
            Complex nu2 = Complex.Sqrt(tmp);

            Complex e1 = Complex.Exp(-2 * nu1 * h);
            Complex R1 = (nu1 - nu2) / (nu1 + nu2) * e1;
            Complex F1 = (k2k2 * nu1 - k1k1 * nu2) / (k2k2 * nu1 + k1k1 * nu2) * e1;

            Complex F0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex F0N = F1 + (k1k1 * nu0 - k0k0 * nu1) / (k1k1 * nu0 + k0k0 * nu1);
                Complex F0D = (k1k1 * nu0 - k0k0 * nu1) / (k1k1 * nu0 + k0k0 * nu1) * F1 + 1;
                F0 = F0N / F0D;

            }
            else
            {
                Complex F0N = k1k1 * nu0 * (k2k2 / (k1k1 * nu2) + h);
                Complex F0D = 2 * F0N + k0k0 * (1 - F1);
                F0 = -1 + 4 * F0N / F0D;
            }

            Complex R0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex R0N = nu0 - nu1 - R1 * (nu0 + nu1);
                Complex R0D = nu0 + nu1 + R1 * (nu0 - nu1);
                R0 = R0N / R0D;

            }
            else
            {
                Complex R0N = nu0 * (1 + h * nu2);
                Complex R0D = 2 * R0N + nu2 * (1 - R1);
                R0 = -1 + 4 * R0N / R0D;
            }
            Complex e_plus = Complex.Exp(-nu0 * (z0 + z));
            Complex e_minus = Complex.Exp(-nu0 * (z0 - z));

            Complex result = wave * R0 * e_plus / nu0 + e_minus / nu0 + (nu0 * (R0 + F0) * e_plus * (nu0 - 1)) / k0k0;  // k0k0?

            fValues[0] = result.Real;        // Re W
            fValues[1] = result.Imaginary;   // Im W
            fValues[2] = fValues[0];     // Dublicate
            fValues[3] = fValues[1];     // Dublicate
        }

        public Complex Ez_eps2_plus_delta_hd(double delta, double lambda)
        {
            Complex tmp = new Complex(lambda * lambda - k0k0, 0);
            Complex nu0 = Complex.Sqrt(tmp);

            this.k2k2 += delta * this.k0k0;
            this.k21 += delta / eps1;

            tmp = new Complex(lambda * lambda - k1k1, 0);
            Complex nu1 = Complex.Sqrt(tmp);

            tmp = new Complex(lambda * lambda - k2k2, 0);
            Complex nu2 = Complex.Sqrt(tmp);

            Complex e1 = Complex.Exp(-2 * nu1 * h);
            Complex R1 = (nu1 - nu2) / (nu1 + nu2) * e1;
            Complex F1 = (k2k2 * nu1 - k1k1 * nu2) / (k2k2 * nu1 + k1k1 * nu2) * e1;

            Complex F0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex F0N = F1 + (k1k1 * nu0 - k0k0 * nu1) / (k1k1 * nu0 + k0k0 * nu1);
                Complex F0D = (k1k1 * nu0 - k0k0 * nu1) / (k1k1 * nu0 + k0k0 * nu1) * F1 + 1;
                F0 = F0N / F0D;

            }
            else
            {
                Complex F0N = k1k1 * nu0 * (k2k2 / (k1k1 * nu2) + h);
                Complex F0D = 2 * F0N + k0k0 * (1 - F1);
                F0 = -1 + 4 * F0N / F0D;
            }

            Complex R0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex R0N = nu0 - nu1 - R1 * (nu0 + nu1);
                Complex R0D = nu0 + nu1 + R1 * (nu0 - nu1);
                R0 = R0N / R0D;

            }
            else
            {
                Complex R0N = nu0 * (1 + h * nu2);
                Complex R0D = 2 * R0N + nu2 * (1 - R1);
                R0 = -1 + 4 * R0N / R0D;
            }
            Complex e_plus = Complex.Exp(-nu0 * (z0 + z));
            Complex e_minus = Complex.Exp(-nu0 * (z0 - z));

            this.k2k2 = this.k0k0 * this.eps2;
            this.k21 = this.eps2 / this.eps1;

            Complex result = wave * R0 * e_plus / nu0 + e_minus / nu0 + (nu0 * (R0 + F0) * e_plus * (nu0 - 1)) / k0k0;

            return result;
        }

        public Complex Ez_eps1_plus_delta_hd(double delta, double lambda)
        {
            Complex tmp = new Complex(lambda * lambda - k0k0, 0);
            Complex nu0 = Complex.Sqrt(tmp);

            this.k1k1 += delta * this.k0k0;
            this.k10 += delta;
            this.k21 = this.eps2 / (this.eps1 + delta);

            tmp = new Complex(lambda * lambda - k1k1, 0);
            Complex nu1 = Complex.Sqrt(tmp);

            tmp = new Complex(lambda * lambda - k2k2, 0);
            Complex nu2 = Complex.Sqrt(tmp);

            Complex e1 = Complex.Exp(-2 * nu1 * h);
            Complex R1 = (nu1 - nu2) / (nu1 + nu2) * e1;
            Complex F1 = (k2k2 * nu1 - k1k1 * nu2) / (k2k2 * nu1 + k1k1 * nu2) * e1;

            Complex F0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex F0N = F1 + (k1k1 * nu0 - k0k0 * nu1) / (k1k1 * nu0 + k0k0 * nu1);
                Complex F0D = (k1k1 * nu0 - k0k0 * nu1) / (k1k1 * nu0 + k0k0 * nu1) * F1 + 1;
                F0 = F0N / F0D;

            }
            else
            {
                Complex F0N = k1k1 * nu0 * (k2k2 / (k1k1 * nu2) + h);
                Complex F0D = 2 * F0N + k0k0 * (1 - F1);
                F0 = -1 + 4 * F0N / F0D;
            }

            Complex R0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex R0N = nu0 - nu1 - R1 * (nu0 + nu1);
                Complex R0D = nu0 + nu1 + R1 * (nu0 - nu1);
                R0 = R0N / R0D;

            }
            else
            {
                Complex R0N = nu0 * (1 + h * nu2);
                Complex R0D = 2 * R0N + nu2 * (1 - R1);
                R0 = -1 + 4 * R0N / R0D;
            }
            Complex e_plus = Complex.Exp(-nu0 * (z0 + z));
            Complex e_minus = Complex.Exp(-nu0 * (z0 - z));

            this.k1k1 = this.k0k0 * this.eps1;
            this.k10 = this.eps1;  // eps1 / eps0
            this.k21 = this.eps2 / this.eps1;

            Complex result = wave * R0 * e_plus / nu0 + e_minus / nu0 + (nu0 * (R0 + F0) * e_plus * (nu0 - 1)) / k0k0;

            return result;
        }

        public Complex Ez_h_plus_delta_hd(double delta, double lambda)
        {
            Complex tmp = new Complex(lambda * lambda - k0k0, 0);
            Complex nu0 = Complex.Sqrt(tmp);

            this.h += delta;

            tmp = new Complex(lambda * lambda - k1k1, 0);
            Complex nu1 = Complex.Sqrt(tmp);

            tmp = new Complex(lambda * lambda - k2k2, 0);
            Complex nu2 = Complex.Sqrt(tmp);

            Complex e1 = Complex.Exp(-2 * nu1 * h);
            Complex R1 = (nu1 - nu2) / (nu1 + nu2) * e1;
            Complex F1 = (k2k2 * nu1 - k1k1 * nu2) / (k2k2 * nu1 + k1k1 * nu2) * e1;

            Complex F0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex F0N = F1 + (k1k1 * nu0 - k0k0 * nu1) / (k1k1 * nu0 + k0k0 * nu1);
                Complex F0D = (k1k1 * nu0 - k0k0 * nu1) / (k1k1 * nu0 + k0k0 * nu1) * F1 + 1;
                F0 = F0N / F0D;

            }
            else
            {
                Complex F0N = k1k1 * nu0 * (k2k2 / (k1k1 * nu2) + h);
                Complex F0D = 2 * F0N + k0k0 * (1 - F1);
                F0 = -1 + 4 * F0N / F0D;
            }

            Complex R0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex R0N = nu0 - nu1 - R1 * (nu0 + nu1);
                Complex R0D = nu0 + nu1 + R1 * (nu0 - nu1);
                R0 = R0N / R0D;

            }
            else
            {
                Complex R0N = nu0 * (1 + h * nu2);
                Complex R0D = 2 * R0N + nu2 * (1 - R1);
                R0 = -1 + 4 * R0N / R0D;
            }
            Complex e_plus = Complex.Exp(-nu0 * (z0 + z));
            Complex e_minus = Complex.Exp(-nu0 * (z0 - z));

            this.h -= delta;

            Complex result = wave * R0 * e_plus / nu0 + e_minus / nu0 + (nu0 * (R0 + F0) * e_plus * (nu0 - 1)) / k0k0;

            return result;
        }

        public Complex dEzdeps2_hd(double delta, double lambda)
        {
            Complex tmp = new Complex(lambda * lambda - k0k0, 0);
            Complex nu0 = Complex.Sqrt(tmp);

            tmp = new Complex(lambda * lambda - k1k1, 0);
            Complex nu1 = Complex.Sqrt(tmp);

            tmp = new Complex(lambda * lambda - k2k2, 0);
            Complex nu2 = Complex.Sqrt(tmp);

            Complex e1 = Complex.Exp(-2 * nu1 * h);
            Complex R1 = (nu1 - nu2) / (nu1 + nu2) * e1;
            Complex F1 = (k2k2 * nu1 - k1k1 * nu2) / (k2k2 * nu1 + k1k1 * nu2) * e1;

            Complex F0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex F0N = F1 + (k1k1 * nu0 - k0k0 * nu1) / (k1k1 * nu0 + k0k0 * nu1);
                Complex F0D = (k1k1 * nu0 - k0k0 * nu1) / (k1k1 * nu0 + k0k0 * nu1) * F1 + 1;
                F0 = F0N / F0D;

            }
            else
            {
                Complex F0N = k1k1 * nu0 * (k2k2 / (k1k1 * nu2) + h);
                Complex F0D = 2 * F0N + k0k0 * (1 - F1);
                F0 = -1 + 4 * F0N / F0D;
            }

            Complex R0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex R0N = nu0 - nu1 - R1 * (nu0 + nu1);
                Complex R0D = nu0 + nu1 + R1 * (nu0 - nu1);
                R0 = R0N / R0D;

            }
            else
            {
                Complex R0N = nu0 * (1 + h * nu2);
                Complex R0D = 2 * R0N + nu2 * (1 - R1);
                R0 = -1 + 4 * R0N / R0D;
            }
            Complex e_plus = Complex.Exp(-nu0 * (z0 + z));
            Complex e_minus = Complex.Exp(-nu0 * (z0 - z));

            Complex Ez_eps2 = wave * R0 * e_plus / nu0 + e_minus / nu0 + (nu0 * (R0 + F0) * e_plus * (nu0 - 1)) / k0k0;

            Complex Ez_eps2_plus_delta_value = Ez_eps2_plus_delta_hd(delta, lambda);

            if (Math.Abs(Ez_eps2.Real - Ez_eps2_plus_delta_value.Real) < 1e-20)
            {
                return 0;
            }
            else
            {
                return (Ez_eps2_plus_delta_value - Ez_eps2) / delta / Math.Sqrt(Math.Pow(Ez_eps2.Real, 2) + Math.Pow(Ez_eps2.Imaginary, 2));
            }
        }

        public Complex dEzdh_hd(double delta, double lambda)
        {
            Complex tmp = new Complex(lambda * lambda - k0k0, 0);
            Complex nu0 = Complex.Sqrt(tmp);

            tmp = new Complex(lambda * lambda - k1k1, 0);
            Complex nu1 = Complex.Sqrt(tmp);

            tmp = new Complex(lambda * lambda - k2k2, 0);
            Complex nu2 = Complex.Sqrt(tmp);

            Complex e1 = Complex.Exp(-2 * nu1 * h);
            Complex R1 = (nu1 - nu2) / (nu1 + nu2) * e1;
            Complex F1 = (k2k2 * nu1 - k1k1 * nu2) / (k2k2 * nu1 + k1k1 * nu2) * e1;

            Complex F0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex F0N = F1 + (k1k1 * nu0 - k0k0 * nu1) / (k1k1 * nu0 + k0k0 * nu1);
                Complex F0D = (k1k1 * nu0 - k0k0 * nu1) / (k1k1 * nu0 + k0k0 * nu1) * F1 + 1;
                F0 = F0N / F0D;

            }
            else
            {
                Complex F0N = k1k1 * nu0 * (k2k2 / (k1k1 * nu2) + h);
                Complex F0D = 2 * F0N + k0k0 * (1 - F1);
                F0 = -1 + 4 * F0N / F0D;
            }

            Complex R0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex R0N = nu0 - nu1 - R1 * (nu0 + nu1);
                Complex R0D = nu0 + nu1 + R1 * (nu0 - nu1);
                R0 = R0N / R0D;

            }
            else
            {
                Complex R0N = nu0 * (1 + h * nu2);
                Complex R0D = 2 * R0N + nu2 * (1 - R1);
                R0 = -1 + 4 * R0N / R0D;
            }
            Complex e_plus = Complex.Exp(-nu0 * (z0 + z));
            Complex e_minus = Complex.Exp(-nu0 * (z0 - z));

            Complex Ez_h = wave * R0 * e_plus / nu0 + e_minus / nu0 + (nu0 * (R0 + F0) * e_plus * (nu0 - 1)) / k0k0;

            Complex Ez_h_plus_delta_value = Ez_h_plus_delta_hd(delta, lambda);

            if (Math.Abs(Ez_h.Real - Ez_h_plus_delta_value.Real) < 1e-20)
            {
                return 0;
            }
            else
            {
                return (Ez_h_plus_delta_value - Ez_h) / delta / Math.Sqrt(Math.Pow(Ez_h.Real, 2) + Math.Pow(Ez_h.Imaginary, 2));
            }
        }

        public Complex dEzdeps1_hd(double delta, double lambda)
        {
            Complex tmp = new Complex(lambda * lambda - k0k0, 0);
            Complex nu0 = Complex.Sqrt(tmp);

            tmp = new Complex(lambda * lambda - k1k1, 0);
            Complex nu1 = Complex.Sqrt(tmp);

            tmp = new Complex(lambda * lambda - k2k2, 0);
            Complex nu2 = Complex.Sqrt(tmp);

            Complex e1 = Complex.Exp(-2 * nu1 * h);
            Complex R1 = (nu1 - nu2) / (nu1 + nu2) * e1;
            Complex F1 = (k2k2 * nu1 - k1k1 * nu2) / (k2k2 * nu1 + k1k1 * nu2) * e1;

            Complex F0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex F0N = F1 + (k1k1 * nu0 - k0k0 * nu1) / (k1k1 * nu0 + k0k0 * nu1);
                Complex F0D = (k1k1 * nu0 - k0k0 * nu1) / (k1k1 * nu0 + k0k0 * nu1) * F1 + 1;
                F0 = F0N / F0D;

            }
            else
            {
                Complex F0N = k1k1 * nu0 * (k2k2 / (k1k1 * nu2) + h);
                Complex F0D = 2 * F0N + k0k0 * (1 - F1);
                F0 = -1 + 4 * F0N / F0D;
            }

            Complex R0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex R0N = nu0 - nu1 - R1 * (nu0 + nu1);
                Complex R0D = nu0 + nu1 + R1 * (nu0 - nu1);
                R0 = R0N / R0D;

            }
            else
            {
                Complex R0N = nu0 * (1 + h * nu2);
                Complex R0D = 2 * R0N + nu2 * (1 - R1);
                R0 = -1 + 4 * R0N / R0D;
            }
            Complex e_plus = Complex.Exp(-nu0 * (z0 + z));
            Complex e_minus = Complex.Exp(-nu0 * (z0 - z));

            Complex Ez_eps1 = wave * R0 * e_plus / nu0 + e_minus / nu0 + (nu0 * (R0 + F0) * e_plus * (nu0 - 1)) / k0k0;

            Complex Ez_eps1_plus_delta_value = Ez_eps1_plus_delta_hd(delta, lambda);

            if (Math.Abs(Ez_eps1.Real - Ez_eps1_plus_delta_value.Real) < 1e-20)
            {
                return 0;
            }
            else
            {
                return (Ez_eps1_plus_delta_value - Ez_eps1) / delta / Math.Sqrt(Math.Pow(Ez_eps1.Real, 2) + Math.Pow(Ez_eps1.Imaginary, 2));
            }
        }

        public void Ez_deps1_int_fun_hd(double lambda, double[] fValues)
        {
            Complex pr = dEzdeps1_hd(0.0001, lambda);
            double result = Math.Sqrt(Math.Pow(pr.Real, 2) + Math.Pow(pr.Imaginary, 2));

            fValues[0] = result;   // Re W
            fValues[1] = 0;   // Im W
            fValues[2] = fValues[0];     // Dublicate
            fValues[3] = fValues[1];     // Dublicate
        }

        public void Ez_deps2_int_fun_hd(double lambda, double[] fValues)
        {
            Complex pr = dEzdeps2_hd(0.0001, lambda);
            double result = Math.Sqrt(Math.Pow(pr.Real, 2) + Math.Pow(pr.Imaginary, 2));

            fValues[0] = result;   // Re W
            fValues[1] = 0;   // Im W
            fValues[2] = fValues[0];     // Dublicate
            fValues[3] = fValues[1];     // Dublicate
        }

        public void Ez_dh_int_fun_hd(double lambda, double[] fValues)
        {
            Complex pr = dEzdh_hd(0.0001, lambda);
            double result = Math.Sqrt(Math.Pow(pr.Real, 2) + Math.Pow(pr.Imaginary, 2));

            fValues[0] = result;   // Re W
            fValues[1] = 0;   // Im W
            fValues[2] = fValues[0];     // Dublicate
            fValues[3] = fValues[1];     // Dublicate
        }






















        public Complex Ez_eps2_plus_delta_vd(double delta, double lambda)
        {
            this.k2k2 += delta * this.k0k0;
            this.k21 += delta / eps1;

            Complex tmp = new Complex(lambda * lambda - k0k0, 0);
            Complex nu0 = Complex.Sqrt(tmp);

            tmp = new Complex(lambda * lambda - k1k1, 0);
            Complex nu1 = Complex.Sqrt(tmp);
            Complex nu1k21 = k21 * nu1;

            tmp = new Complex(lambda * lambda - k2k2, 0);
            Complex nu2 = Complex.Sqrt(tmp);

            Complex e1 = Complex.Exp(-2 * nu1 * h);
            Complex F1 = (nu1k21 - nu2) / (nu1k21 + nu2) * e1;

            Complex G0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex G0N = k10 * (F1 + 1);
                Complex G0D = nu0 * G0N - nu1 * (F1 - 1);
                G0 = 2 * G0N / G0D;

            }
            else
            {
                Complex G0N = k10 * (k21 + h * nu2);
                Complex G0D = 2 * G0N * nu0 + nu2 * (1 - F1);
                G0 = 4 * G0N / G0D;
            }
            Complex e0 = Complex.Exp(-nu0 * (z + z0));
            double nu02 = lambda * lambda - k0k0;
            Complex result = wave * G0 * e0 * (1 + nu02 / k0k0);   // k0k0??

            this.k2k2 = this.k0k0 * this.eps2;
            this.k21 = this.eps2 / this.eps1;

            return result;
        }

        public Complex Ez_eps1_plus_delta_vd(double delta, double lambda)
        {

            this.k1k1 += delta * this.k0k0;
            this.k10 += delta;
            this.k21 = this.eps2 / (this.eps1 + delta);

            Complex tmp = new Complex(lambda * lambda - k0k0, 0);
            Complex nu0 = Complex.Sqrt(tmp);

            tmp = new Complex(lambda * lambda - k1k1, 0);
            Complex nu1 = Complex.Sqrt(tmp);
            Complex nu1k21 = k21 * nu1;

            tmp = new Complex(lambda * lambda - k2k2, 0);
            Complex nu2 = Complex.Sqrt(tmp);

            Complex e1 = Complex.Exp(-2 * nu1 * h);
            Complex F1 = (nu1k21 - nu2) / (nu1k21 + nu2) * e1;

            Complex G0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex G0N = k10 * (F1 + 1);
                Complex G0D = nu0 * G0N - nu1 * (F1 - 1);
                G0 = 2 * G0N / G0D;

            }
            else
            {
                Complex G0N = k10 * (k21 + h * nu2);
                Complex G0D = 2 * G0N * nu0 + nu2 * (1 - F1);
                G0 = 4 * G0N / G0D;
            }
            Complex e0 = Complex.Exp(-nu0 * (z + z0));
            double nu02 = lambda * lambda - k0k0;

            this.k1k1 = this.k0k0 * this.eps1;
            this.k10 = this.eps1;  // eps1 / eps0
            this.k21 = this.eps2 / this.eps1;

            Complex result = wave * G0 * e0 * (1 + nu02 / k0k0);   // k0k0??

            return result;
        }

        public Complex Ez_h_plus_delta_vd(double delta, double lambda)
        {
            this.h += delta;

            Complex tmp = new Complex(lambda * lambda - k0k0, 0);
            Complex nu0 = Complex.Sqrt(tmp);

            tmp = new Complex(lambda * lambda - k1k1, 0);
            Complex nu1 = Complex.Sqrt(tmp);
            Complex nu1k21 = k21 * nu1;

            tmp = new Complex(lambda * lambda - k2k2, 0);
            Complex nu2 = Complex.Sqrt(tmp);

            Complex e1 = Complex.Exp(-2 * nu1 * h);
            Complex F1 = (nu1k21 - nu2) / (nu1k21 + nu2) * e1;

            Complex G0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex G0N = k10 * (F1 + 1);
                Complex G0D = nu0 * G0N - nu1 * (F1 - 1);
                G0 = 2 * G0N / G0D;

            }
            else
            {
                Complex G0N = k10 * (k21 + h * nu2);
                Complex G0D = 2 * G0N * nu0 + nu2 * (1 - F1);
                G0 = 4 * G0N / G0D;
            }
            Complex e0 = Complex.Exp(-nu0 * (z + z0));

            double nu02 = lambda * lambda - k0k0;
            Complex result = wave * G0 * e0 * (1 + nu02 / k0k0);   // k0k0??

            this.h -= delta;

            return result;
        }

        public Complex dEzdeps2_vd(double delta, double lambda)
        {
            Complex tmp = new Complex(lambda * lambda - k0k0, 0);
            Complex nu0 = Complex.Sqrt(tmp);

            tmp = new Complex(lambda * lambda - k1k1, 0);
            Complex nu1 = Complex.Sqrt(tmp);
            Complex nu1k21 = k21 * nu1;

            tmp = new Complex(lambda * lambda - k2k2, 0);
            Complex nu2 = Complex.Sqrt(tmp);

            Complex e1 = Complex.Exp(-2 * nu1 * h);
            Complex F1 = (nu1k21 - nu2) / (nu1k21 + nu2) * e1;

            Complex G0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex G0N = k10 * (F1 + 1);
                Complex G0D = nu0 * G0N - nu1 * (F1 - 1);
                G0 = 2 * G0N / G0D;

            }
            else
            {
                Complex G0N = k10 * (k21 + h * nu2);
                Complex G0D = 2 * G0N * nu0 + nu2 * (1 - F1);
                G0 = 4 * G0N / G0D;
            }
            Complex e0 = Complex.Exp(-nu0 * (z + z0));
            double nu02 = lambda * lambda - k0k0;
            Complex Ez_eps2 = wave * G0 * e0 * (1 + nu02 / k0k0);   // k0k0??

            Complex Ez_eps2_plus_delta_value = Ez_eps2_plus_delta_vd(delta, lambda);

            if (Math.Abs(Ez_eps2.Real - Ez_eps2_plus_delta_value.Real) < 1e-20)
            {
                return 0;
            }
            else
            {
                return (Ez_eps2_plus_delta_value - Ez_eps2) / delta / Math.Sqrt(Math.Pow(Ez_eps2.Real, 2) + Math.Pow(Ez_eps2.Imaginary, 2));
            }
        }

        public Complex dEzdh_vd(double delta, double lambda)
        {
            Complex tmp = new Complex(lambda * lambda - k0k0, 0);
            Complex nu0 = Complex.Sqrt(tmp);

            tmp = new Complex(lambda * lambda - k1k1, 0);
            Complex nu1 = Complex.Sqrt(tmp);
            Complex nu1k21 = k21 * nu1;

            tmp = new Complex(lambda * lambda - k2k2, 0);
            Complex nu2 = Complex.Sqrt(tmp);

            Complex e1 = Complex.Exp(-2 * nu1 * h);
            Complex F1 = (nu1k21 - nu2) / (nu1k21 + nu2) * e1;

            Complex G0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex G0N = k10 * (F1 + 1);
                Complex G0D = nu0 * G0N - nu1 * (F1 - 1);
                G0 = 2 * G0N / G0D;

            }
            else
            {
                Complex G0N = k10 * (k21 + h * nu2);
                Complex G0D = 2 * G0N * nu0 + nu2 * (1 - F1);
                G0 = 4 * G0N / G0D;
            }
            Complex e0 = Complex.Exp(-nu0 * (z + z0));

            double nu02 = lambda * lambda - k0k0;
            Complex Ez_h = wave * G0 * e0 * (1 + nu02 / k0k0);   // k0k0??

            Complex Ez_h_plus_delta_value = Ez_h_plus_delta_vd(delta, lambda);

            if (Math.Abs(Ez_h.Real - Ez_h_plus_delta_value.Real) < 1e-20)
            {
                return 0;
            }
            else
            {
                return (Ez_h_plus_delta_value - Ez_h) / delta / Math.Sqrt(Math.Pow(Ez_h.Real, 2) + Math.Pow(Ez_h.Imaginary, 2));
            }
        }

        public Complex dEzdeps1_vd(double delta, double lambda)
        {
            Complex tmp = new Complex(lambda * lambda - k0k0, 0);
            Complex nu0 = Complex.Sqrt(tmp);

            tmp = new Complex(lambda * lambda - k1k1, 0);
            Complex nu1 = Complex.Sqrt(tmp);
            Complex nu1k21 = k21 * nu1;

            tmp = new Complex(lambda * lambda - k2k2, 0);
            Complex nu2 = Complex.Sqrt(tmp);

            Complex e1 = Complex.Exp(-2 * nu1 * h);
            Complex F1 = (nu1k21 - nu2) / (nu1k21 + nu2) * e1;

            Complex G0;

            if (Math.Abs(lambda * lambda - k1k1) > 1.0E-16)
            {
                Complex G0N = k10 * (F1 + 1);
                Complex G0D = nu0 * G0N - nu1 * (F1 - 1);
                G0 = 2 * G0N / G0D;

            }
            else
            {
                Complex G0N = k10 * (k21 + h * nu2);
                Complex G0D = 2 * G0N * nu0 + nu2 * (1 - F1);
                G0 = 4 * G0N / G0D;
            }
            Complex e0 = Complex.Exp(-nu0 * (z + z0));

            double nu02 = lambda * lambda - k0k0;
            Complex Ez_eps1 = wave * G0 * e0 * (1 + nu02 / k0k0);   // k0k0??
            Complex Ez_eps1_plus_delta_value = Ez_eps1_plus_delta_vd(delta, lambda);

            if (Math.Abs(Ez_eps1.Real - Ez_eps1_plus_delta_value.Real) < 1e-20)
            {
                return 0;
            }
            else
            {
                return (Ez_eps1_plus_delta_value - Ez_eps1) / delta / Math.Sqrt(Math.Pow(Ez_eps1.Real, 2) + Math.Pow(Ez_eps1.Imaginary, 2));
            }
        }

        public void Ez_deps1_int_fun_vd(double lambda, double[] fValues)
        {
            Complex pr = dEzdeps1_vd(0.0001, lambda);
            double result = Math.Sqrt(Math.Pow(pr.Real, 2) + Math.Pow(pr.Imaginary, 2));

            fValues[0] = result;   // Re W
            fValues[1] = 0;   // Im W
            fValues[2] = fValues[0];     // Dublicate
            fValues[3] = fValues[1];     // Dublicate
        }

        public void Ez_deps2_int_fun_vd(double lambda, double[] fValues)
        {
            Complex pr = dEzdeps2_vd(0.0001, lambda);
            double result = Math.Sqrt(Math.Pow(pr.Real, 2) + Math.Pow(pr.Imaginary, 2));

            fValues[0] = result;   // Re W
            fValues[1] = 0;   // Im W
            fValues[2] = fValues[0];     // Dublicate
            fValues[3] = fValues[1];     // Dublicate
        }

        public void Ez_dh_int_fun_vd(double lambda, double[] fValues)
        {
            Complex pr = dEzdh_vd(0.0001, lambda);
            double result = Math.Sqrt(Math.Pow(pr.Real, 2) + Math.Pow(pr.Imaginary, 2));

            fValues[0] = result;   // Re W
            fValues[1] = 0;   // Im W
            fValues[2] = fValues[0];     // Dublicate
            fValues[3] = fValues[1];     // Dublicate
        }
        public bool CalculateFields_byAnderson(double[] ros, List<double[]> res)
        {
            Integrand func = Ez_dh_int_fun_vd;
            if (FieldComponent.FC == FldComp.Ez_hd) func = Ez_dh_int_fun_hd;
           
            try
            {
                double eps = 1.0e-5; // не используется
                double[] resH = new double[4];
                double[] c = new double[4];
                double[] cmax = new double[4];
                for (int jro = 0; jro < ros.Length; jro++)
                {
                    double rr = ros[jro];
                    HankelTransformsByAnderson.call(func, rr, eps, 4, 0, resH, c, cmax);
                    res[0][jro] = resH[0];
                    res[1][jro] = resH[1];
                    res[2][jro] = Math.Sqrt( resH[0] * resH[0] + resH[1] * resH[1]);
                    //res[2][jro] = resH[2] / k0k0;
                    //res[3][jro] = resH[3] / k0k0; 
                }
                return true;
            }
            catch (Exception ex)
            {
                throw new Exception("Error in CalculateFields()", ex);
            }
        }

        public bool CalculateFields_bySplines(double[] ros, List<double[]> res)
        {
            Integrand func = Ez_vd_int_fun;
            if (FieldComponent.FC == FldComp.Ez_hd) func = Ez_hd_int_fun;

            try
            {
                double eps = 1.0e-5; 
                double[] resH = new double[4];
              
                for (int jro = 0; jro < ros.Length; jro++)
                {
                    double rr = ros[jro];
                    IntegralJ0_BySpline(func, rr, eps, resH);
                    res[0][jro] = Math.Sqrt(resH[0] * resH[0] + resH[1] * resH[1]);
                    res[1][jro] = resH[1];
                    res[2][jro] = Math.Sqrt(resH[0] * resH[0] + resH[1] * resH[1]);
                    //res[2][jro] = resH[2] / k0k0;
                    //res[3][jro] = resH[3] / k0k0; 
                }
                return true;
            }
            catch (Exception ex)
            {
                throw new Exception("Error in CalculateFields_bySplines", ex);
            }
        }
        public override string ToString()
        {
            return $"\nlam0om1 = {EMC.lam0om1} \nh= {h} \nwave = {wave}" +
                $"\nk0k0 = {k0k0} \nk1k1= {k1k1} \nk2k2 = {k2k2}";
        }

        //public void TestIntegral(double ro, ref double res, ref double exact, ref int n)
        //{
        //    try
        //    {
        //        int niter = 10;
        //        double _eps = 1.0E-5;
        //        double a = 0;
        //        double b = 1;
        //        int npt = 5;

        //        double h = (b - a) / (npt - 1);
        //        double[] points = new double[npt];
        //        double[] func = new double[npt];
        //        for (int j = 0; j < npt; j++)
        //        {   points[j] = a + h * j;
        //            func[j] = Math.Exp(points[j]);
        //        }
        //        double _res = 0;
        //        double _res2 = 0;
        //        IntegralSpline(npt, points, func, ref _res);

        //        for (int iter = 0; iter < niter; iter++)
        //        {
        //         int npt2 = 2 * npt;
        //         double h2 = (b - a) / (npt2 - 1);
        //         double[] points2 = new double[npt2];
        //         double[] func2 = new double[npt2];

        //         for (int j = 0; j < npt2; j++)
        //         {
        //            points2[j] = a + h2 * j;
        //            func2[j] = Math.Exp(points2[j]);
        //         }

        //            IntegralSpline(npt2, points2, func2, ref _res2);
        //            if (Math.Abs(_res - _res2) < _eps)
        //            {
        //                n = iter;
        //                break; 
        //            }
        //            else
        //            {
        //                npt = npt2;
        //                _res = _res2;
        //            }
        //        }
        //        res = _res2;
        //        exact = Math.Exp(1) - 1;
        //    }
        //    catch(Exception ex)
        //    {
        //        throw new Exception("DipoleFields.TestIntegral", ex);
        //    }
        //}

        void IntegralJ0_BySpline(Integrand func, double ro, double eps, double[] resH)
        {
            double[] res = new double[2];
            double lam0 = 0;
            double lam1 = k0;
            double lam2 = k1;
            double lam3 = k2;
            if (k1 > k2)
            {   lam2 = k2;
                lam3 = k1;
            }
            double lam4 = 4.0 / eps; 

            int n = 0;
            _Integral_J0(func, ro, lam0, lam1, res, eps, ref n);
            resH[0] = res[0];
            resH[1] = res[1];
            TestString.Add($"lam0 = {lam0} lam1 = {lam1} n = {n}" + $"\nres[0] = {res[0]} res[1] = {res[1]}");
           
            _Integral_J0(func, ro, lam1, lam2, res, eps, ref n);
            resH[0] += res[0];
            resH[1] += res[1];
            TestString.Add($"lam1 = {lam1} lam2 = {lam2} n = {n}" + $"\nres[0] = {res[0]} res[1] = {res[1]}");

            _Integral_J0(func, ro, lam2, lam3, res, eps, ref n);
            resH[0] += res[0];
            resH[1] += res[1];
            TestString.Add($"lam2 = {lam2} lam3 = {lam3} n = {n}" + $"\nres[0] = {res[0]} res[1] = {res[1]}");

            _Integral_J0(func, ro, lam3, lam4, res, eps, ref n);
            resH[0] += res[0];
            resH[1] += res[1];
            TestString.Add($"lam3 = {lam3} lam4 = {lam4} n = {n}" + $"\nres[0] = {res[0]} res[1] = {res[1]}");
        }

        void _Integral_J0 (Integrand func, double ro, double lamA, double lamB, double[] res, double _eps, ref int n)
        {
            int _n = 0;
            try
            {
                const int niter = 5;
                int npt = 20;

                double[] fValues = new double[4];
               
                double h = (lamB - lamA) / (npt - 1);
                double[] points = new double[npt];
                double[] fv1 = new double[npt];
                double[] fv2 = new double[npt];

                for (int j = 0; j < npt; j++)
                {
                    points[j] = lamA + h * j;
                    func(points[j], fValues);
                    fv1[j] = fValues[0];
                    fv2[j] = fValues[1];
                }
                double _res1 = 0;
                double _res2 = 0;
                double _resp1 = 0;
                double _resp2 = 0;

                IntegralSpline_J0(ro, npt, points, fv1, ref _resp1);
                IntegralSpline_J0(ro, npt, points, fv2, ref _resp2);
               
                for (int iter = 0; iter < niter; iter++)
                {
                    int npt2 = 2 * npt;
                    double h2 = (lamB - lamA) / (npt2 - 1);
                    double[] points2 = new double[npt2];
                    fv1 = new double[npt2];
                    fv2 = new double[npt2];

                    for (int j = 0; j < npt2; j++)
                    {
                        points2[j] = lamA + h2 * j;
                        func(points2[j], fValues);
                        fv1[j] = fValues[0];
                        fv2[j] = fValues[1];
                    }

                    IntegralSpline_J0(ro, npt2, points2, fv1, ref _res1);
                    IntegralSpline_J0(ro, npt2, points2, fv2, ref _res2);

                    //    if (Math.Abs(_res1 - _resp1) < _eps && Math.Abs(_res2 - _resp2) < _eps)
                    double del1 = Math.Abs((_res1 - _resp1) / _res1);
                    double del2 = Math.Abs((_res2 - _resp2) / _res2);
                    if (del1 < _eps && del2 < _eps)
                    {
                        _n = iter;
                        break;
                    }
                    else
                    {    npt = npt2;
                        _resp1 = _res1;
                        _resp2 = _res2;
                    }
                }
                res[0] = _res1;
                res[1] = _res2;
                n = _n;
            }
            catch (Exception ex)
            {
                throw new Exception("DipoleFields.TestIntegral", ex);
            }
        }

        [DllImport("..\\..\\..\\x64\\DEBUG\\DLL_CPP.dll", CallingConvention = CallingConvention.Cdecl)]   // NetFramework
        public static extern
        void Test(double _eps1, double _eps2, double _h,
                  double _ro, double _z0, double _z, double _wave,
                  int nlam, double[] lam, double[] res);
        [DllImport("..\\..\\..\\x64\\DEBUG\\DLL_CPP.dll", CallingConvention = CallingConvention.Cdecl)]   // NetFramework
        public static extern
        double IntegralSpline(int length, double[] points, double[] func, ref double integral);
        
        [DllImport("..\\..\\..\\x64\\DEBUG\\DLL_CPP.dll", CallingConvention = CallingConvention.Cdecl)]   // NetFramework
        public static extern
        double IntegralSpline_J0(double ro, int length, double[] points, double[] func, ref double integral);

        [DllImport("..\\..\\..\\x64\\DEBUG\\DLL_CPP.dll", CallingConvention = CallingConvention.Cdecl)]   // NetFramework
        public static extern
       double IntegralSpline_J1(double ro, int length, double[] points, double[] func, ref double integral);
    }
}
