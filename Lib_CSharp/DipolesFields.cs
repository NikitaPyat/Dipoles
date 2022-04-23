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
       public string TestString { get; set; }
       
        double eps1;
        double eps2;
        double h;
        double ro; 
        double z0;
        double z;
        double wave ;
        double k0k0;
        double k1k1;
        double k2k2;
        double k10;
        double k21;
      
        public DipolesFields( ModelParams modelParams)
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
          
            k0k0 = 4 * Math.PI * Math.PI / wave / wave;
            k1k1 = k0k0 * eps1;
            k2k2 = k0k0 * eps2;
            k10 = eps1;  // eps1 / eps0
            k21 = eps2 / eps1;

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

        public void Ez_vd_int_fun (double lambda, double[] fValues)
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

        public Complex Ez_eps2_plus_delta(double delta, double lambda)
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

        public Complex Ez_eps1_plus_delta(double delta, double lambda)
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

        public Complex Ez_h_plus_delta(double delta, double lambda)
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

        public Complex dEzdeps2(double delta, double lambda)
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

            Complex Ez_eps2_plus_delta_value = Ez_eps2_plus_delta(delta, lambda);

            if (Math.Abs(Ez_eps2.Real - Ez_eps2_plus_delta_value.Real) < 1e-20)
            {
                return 0;
            }
            else
            {
                return (Ez_eps2_plus_delta_value - Ez_eps2) / delta / Math.Sqrt(Math.Pow(Ez_eps2.Real, 2) + Math.Pow(Ez_eps2.Imaginary, 2));
            }
        }

        public Complex dEzdh(double delta, double lambda)
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

            Complex Ez_h_plus_delta_value = Ez_h_plus_delta(delta, lambda);

            if (Math.Abs(Ez_h.Real - Ez_h_plus_delta_value.Real) < 1e-20)
            {
                return 0;
            }
            else
            {
                return (Ez_h_plus_delta_value - Ez_h) / delta / Math.Sqrt(Math.Pow(Ez_h.Real, 2) + Math.Pow(Ez_h.Imaginary, 2));
            }
        }

        public Complex dEzdeps1(double delta, double lambda)
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

            Complex Ez_eps1 =  wave * R0 * e_plus / nu0 + e_minus / nu0 + (nu0 * (R0 + F0) * e_plus * (nu0 - 1)) / k0k0;

            Complex Ez_eps1_plus_delta_value = Ez_eps1_plus_delta(delta, lambda);

            if (Math.Abs(Ez_eps1.Real - Ez_eps1_plus_delta_value.Real) < 1e-20)
            {
                return 0;
            } else
            {
                return (Ez_eps1_plus_delta_value - Ez_eps1) / delta / Math.Sqrt(Math.Pow(Ez_eps1.Real, 2) + Math.Pow(Ez_eps1.Imaginary, 2));
            }
        }

        public void Ez_deps1_hd_int_fun(double lambda, double[] fValues)
        {
            Complex pr = dEzdeps1(0.0001, lambda);
            double result = Math.Sqrt(Math.Pow(pr.Real, 2) + Math.Pow(pr.Imaginary, 2)); 

            fValues[0] = result;   // Re W
            fValues[1] = 0;   // Im W
            fValues[2] = fValues[0];     // Dublicate
            fValues[3] = fValues[1];     // Dublicate
        }

        public void Ez_deps2_hd_int_fun(double lambda, double[] fValues)
        {
            Complex pr = dEzdeps2(0.0001, lambda);
            double result = Math.Sqrt(Math.Pow(pr.Real, 2) + Math.Pow(pr.Imaginary, 2));

            fValues[0] = result;   // Re W
            fValues[1] = 0;   // Im W
            fValues[2] = fValues[0];     // Dublicate
            fValues[3] = fValues[1];     // Dublicate
        }

        public void Ez_dh_hd_int_fun(double lambda, double[] fValues)
        {
            Complex pr = dEzdh(0.0001, lambda);
            double result = Math.Sqrt(Math.Pow(pr.Real, 2) + Math.Pow(pr.Imaginary, 2));

            fValues[0] = result;   // Re W
            fValues[1] = 0;   // Im W
            fValues[2] = fValues[0];     // Dublicate
            fValues[3] = fValues[1];     // Dublicate
        }
        public bool CalculateFields(double[] ros, List<double[]> res)
        {
            try
            {
                double eps = 1.0e-5; // не используется
                double[] resH = new double[4];
                double[] c = new double[4];
                double[] cmax = new double[4];
                for (int jro = 0; jro < ros.Length; jro++)
                {
                    double rr = ros[jro];
                    HankelTransformsByAnderson.call(Ez_dh_hd_int_fun, rr, eps, 4, 0, resH, c, cmax);
                    res[0][jro] = resH[0];
                    res[1][jro] = resH[1];
                    res[2][jro] = resH[2] / k0k0;
                    res[3][jro] = resH[3] / k0k0; 
                }
                //TestString = $"IW_re = {res[0]} IW_im = { res[1]}" +
                //             $"IWzz_re = {res[2]} IWzz_im = { res[3]}";
                //Wrap_Splines_Lab2_2022V(var, n, mData.x, mData.mData, bc1_1, bc1_2,
                //                         spParams.nUni, xUni, splinesEnds_1, ref return_value);
                //Wrap_Splines_Lab2_2022V(var, n, mData.x, mData.mData, bc2_1, bc2_2,
                //                    spParams.nUni, xUni, splinesEnds_2, ref return_value);
                //Wrap_Splines_Lab2_2022V(0, n, mData.x, mData.mData, bc1_1, bc1_2,                  // free
                //                         spParams.nUni, xUni, splinesFree, ref return_value);
                return true;
            }
            catch (Exception ex)
            {
                throw new Exception("Error in CalculateFields()", ex);
            }
        }

        public override string ToString()
        {
            return $"\nlam0om1 = {EMC.lam0om1} \nh= {h} \nwave = {wave}" +
                $"\nk0k0 = {k0k0} \nk1k1= {k1k1} \nk2k2 = {k2k2}";
        }

        [DllImport("..\\..\\..\\x64\\DEBUG\\DLL_CPP.dll", CallingConvention = CallingConvention.Cdecl)]   // NetFramework
       // [DllImport("..\\..\\..\\..\\x64\\DEBUG\\DLL_CPP.dll", CallingConvention = CallingConvention.Cdecl)] // Core
        public static extern
        void Wrap_Splines_Lab2_2022(int n, double[] x, double[] mData, double bc1_1, double bc1_2,
                                    int nN, double[] xUni, double[] res, ref int return_value);

        [DllImport("..\\..\\..\\x64\\DEBUG\\DLL_CPP.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern
        void Wrap_Splines_Lab2_2022V(int variant, int n, double[] x, double[] mData, double bc1_1, double bc1_2,
                                     int nN, double[] xUni, double[] res, ref int return_value);

        [DllImport("..\\..\\..\\x64\\DEBUG\\DLL_CPP.dll", CallingConvention = CallingConvention.Cdecl)]   // NetFramework
        public static extern
        void Test(double _eps1, double _eps2, double _h,
                  double _ro, double _z0, double _z, double _wave,
                  int nlam, double[] lam, double[] res);
    }
}
