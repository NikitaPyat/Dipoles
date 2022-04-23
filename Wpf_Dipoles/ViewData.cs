using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.ComponentModel;
using System.Windows;
using Lib_CSharp;
using System.Collections.ObjectModel;
using System.Windows.Forms.DataVisualization.Charting;


namespace Wpf_Dipoles
{

    class ViewData 
    {
        public int variant { get; set; }
        public bool IsUpdated { get; set; }
        public ModelParams modelParams { get; set; }
        public LambdaParams lamParams { get; set; }
        public FieldComponentParams fieldCompParams { get; set; }
        public DipolesFields dipolesFields { get; set; }
        public DataForOutput dataForOutputL { get; set; }
        public DataForOutput dataForOutputR { get; set; }
       
        public ViewData()
        {
            lamParams = new LambdaParams();
            modelParams = new ModelParams();
            fieldCompParams = new FieldComponentParams();
            dipolesFields = new DipolesFields( modelParams);
            IsUpdated = false;
        }
    
        public void CalculateDraw_Flam(Chart chart1, Chart chart2)
        { 
            try
            {   int nlam = lamParams.n;
                double[] lam = new double[nlam];
                double hlam = (lamParams.lamR - lamParams.lamL) / (nlam - 1);
                for (int j = 0; j < nlam; j++) lam[j] = lamParams.lamL + j * hlam;
                List<double[]> rList = new List<double[]>(4);
                rList.Add(new double[nlam]);
                rList.Add(new double[nlam]);
                rList.Add(new double[nlam]);
                rList.Add(new double[nlam]);

                dipolesFields.CalculateFlam_CSharp(lam, rList);
                List<double[]> listL = new List<double[]>();
                List<double[]> listR = new List<double[]>();
                List<string> legendsL = new List<string>();
                List<string> legendsR = new List<string>();
                string Title = "";
                if (lamParams.F == DpF.W)
                {
                    Title = "W";
                    listL.Add(rList[0]);
                    legendsL.Add("Re W");
                    listR.Add(rList[1]);
                    legendsR.Add("Im W");
                }
                else
                {
                    Title = "Wzz";
                    listL.Add(rList[2]);
                    legendsL.Add("Re Wzz");
                    listR.Add(rList[3]);
                    legendsR.Add("Im Wzz");
                }

                ChartView.DrawChart(chart1, Title, "lam", Title, lam, listL, legendsL);
                ChartView.DrawChart(chart2, Title, "lam", Title, lam, listR, legendsR);
                //dataForOutputL = new DataForOutput("Re W", "lam", "Re W", dipolesFields.lamData.lam, ListFL, legendsL);
                //dataForOutputR = new DataForOutput("Im W", "lam", "Im W", dipolesFields.lamData.lam, ListFR, legendsR);

                // dataForOutputL.infoList.Add(dipolesFields.TestString);
                // dataForOutputL.DataToStringList(dipolesFields.lamData.lam, dipolesFields.Re, dipolesFields.Re_CPP, "Re");

                // dataForOutputR.DataToStringList(dipolesFields.lamData.lam, dipolesFields.Im, dipolesFields.Im_CPP, "Im");
                IsUpdated = true;
            }
            catch(Exception ex)
            {
                MessageBox.Show("Error in Calculate_Flam()" + ex.ToString());
            }
        }
       public void CalculateFields(Chart chart1, Chart chart2)
        {
            try
            {
                int nro = fieldCompParams.nro;
                double[] ro = new double[nro];  // в длинах волн в воздухе
                double hro = (fieldCompParams.roMax - fieldCompParams.roMin) / (nro - 1) ;
                double[] ros = new double[nro];
              
                for (int j = 0; j < nro; j++)
                { ro[j] = fieldCompParams.roMin + j * hro;
                  ros[j] = ro[j] * modelParams.wave; // в метрах для вычислений
                }
                List<double[]> rList = new List<double[]>(4);
                rList.Add(new double[nro]);
                rList.Add(new double[nro]);
                rList.Add(new double[nro]);
                rList.Add(new double[nro]);

                dipolesFields.CalculateFields(ros, rList);
                List<double[]> listL = new List<double[]>();
                List<double[]> listR = new List<double[]>();
                List<string> legendsL = new List<string>();
                List<string> legendsR = new List<string>();
                string Title = "";
                if (lamParams.F == DpF.W)
                {
                    Title = "J0_W";
                    listL.Add(rList[0]);
                    legendsL.Add("Re J0_W");
                    listR.Add(rList[1]);
                    legendsR.Add("Im J0_W");
                }
                else
                {
                    Title = "J0_Wzz";
                    listL.Add(rList[2]);
                    legendsL.Add("Re J0_Wzz");
                    listR.Add(rList[3]);
                    legendsR.Add("Im j0_Wzz");
                }

                ChartView.DrawChart(chart1, Title, "ro в длинах волн в воздухе", Title, ro, listL, legendsL);
                ChartView.DrawChart(chart2, Title, "ro в длинах волн в воздухе", Title, ro, listR, legendsR);
                //dataForOutputL = new DataForOutput("Re W", "lam", "Re W", dipolesFields.lamData.lam, ListFL, legendsL);
                //dataForOutputR = new DataForOutput("Im W", "lam", "Im W", dipolesFields.lamData.lam, ListFR, legendsR);

                // dataForOutputL.infoList.Add(dipolesFields.TestString);
                // dataForOutputL.DataToStringList(dipolesFields.lamData.lam, dipolesFields.Re, dipolesFields.Re_CPP, "Re");

                // dataForOutputR.DataToStringList(dipolesFields.lamData.lam, dipolesFields.Im, dipolesFields.Im_CPP, "Im");
                IsUpdated = true;
            }
            catch (Exception ex)
            {
                MessageBox.Show("Error in Calculate_Flam()" + ex.ToString());
            }
        }

        public override string ToString()
        {
            return lamParams.ToString() + "\n" + modelParams.ToString();
        }
    }
}
