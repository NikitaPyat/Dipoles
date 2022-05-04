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
            dipolesFields = new DipolesFields( modelParams, fieldCompParams);
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
       public void _CalculateFields(Chart chart1, Chart chart2, Chart chart3)
        {
            try
            {
                MessageBox.Show("CalculateFields");
                dipolesFields = new DipolesFields(modelParams, fieldCompParams);

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

                dipolesFields.CalculateFields_byAnderson(ros, rList);

                //for (int j = 0; j < nro; j++)
                //{
                //    rList[0][j] = rList[0][j] / ro[j] / ro[j];
                //    rList[1][j] = rList[1][j] / ro[j] / ro[j];
                //    rList[2][j] = rList[2][j] / ro[j] / ro[j];
                //}

                List<double[]> listL = new List<double[]>();
                listL.Add(rList[0]);
                List<double[]> listR = new List<double[]>();
                listR.Add(rList[1]);
                List<double[]> listB = new List<double[]>();
                listB.Add(rList[2]); 
                List<string> legendsL = new List<string>();
                List<string> legendsR = new List<string>();
                List<string> legendsB = new List<string>();

                string Title1 = "Re " + fieldCompParams.FC;
                string Title2 = "Im " + fieldCompParams.FC;
                string Title3 = "Mod " + fieldCompParams.FC;
               
                ChartView.DrawChart(chart1, Title1, "ro в длинах волн в воздухе", "", ro, listL, legendsL);
                ChartView.DrawChart(chart2, Title2, "ro в длинах волн в воздухе", "", ro, listR, legendsR);
                ChartView.DrawChart(chart3, Title3, "ro в длинах волн в воздухе", "", ro, listB, legendsB);
               
                IsUpdated = true;
            }
            catch (Exception ex)
            {
                MessageBox.Show("Error in Calculate_Flam()" + ex.ToString());
            }
        }

        public void CalculateFields_Drv(Chart chart1, Chart chart2, Chart chart3)
        {
          //  MessageBox.Show("CalculateFields_Drv");
            try
            {
                List<double[]> listL = new List<double[]>();
                List<double[]> listR = new List<double[]>();
                List<double[]> listB = new List<double[]>();
               
                List<string> legendsL = new List<string>();
                List<string> legendsR = new List<string>();
                List<string> legendsB = new List<string>();
                
                double delta = 0.001;
               
 
                    double delta_eps1 = 0;
                    double delta_eps2 = 0;
                    double delta_h = 0;

                    dipolesFields = new DipolesFields(modelParams, fieldCompParams);

                    int nro = fieldCompParams.nro;
                    double[] ro = new double[nro];  // в длинах волн в воздухе
                    double hro = (fieldCompParams.roMax - fieldCompParams.roMin) / (nro - 1);
                    double[] ros = new double[nro];

                    for (int j = 0; j < nro; j++)
                    {
                        ro[j] = fieldCompParams.roMin + j * hro;
                        ros[j] = ro[j] * modelParams.wave; // в метрах для вычислений
                    }
                    List<double[]> rList = new List<double[]>(4);
                    rList.Add(new double[nro]);
                    rList.Add(new double[nro]);
                    rList.Add(new double[nro]);
                    rList.Add(new double[nro]);

                    dipolesFields.CalculateFields_byAnderson(ros, rList);
            
                for (int jdelta = 0; jdelta < 3; jdelta++)
                 {
                    if (fieldCompParams.drvParam == DrvParam.eps1) delta_eps1 = delta;
                    else if (fieldCompParams.drvParam == DrvParam.eps2) delta_eps2 = delta;
                    else delta_h = delta;
                 
                    ModelParams modelParams_delta =
                       new ModelParams(modelParams.eps1 * (1 + delta_eps1), modelParams.eps2 * (1 + delta_eps2), modelParams.h * (1 + delta_h),
                                       modelParams.ro, modelParams.z0, modelParams.z, modelParams._omega);

                    DipolesFields dipolesFields_delta = new DipolesFields(modelParams_delta, dipolesFields.FieldComponent);
                  
                    List<double[]> rList_delta = new List<double[]>(4);
                    rList_delta.Add(new double[nro]);
                    rList_delta.Add(new double[nro]);
                    rList_delta.Add(new double[nro]);
                    rList_delta.Add(new double[nro]);
                    dipolesFields_delta.CalculateFields_byAnderson(ros, rList_delta);

                    for (int j = 0; j < nro; j++)
                    {
                        rList_delta[0][j] = Math.Abs((rList_delta[0][j] - rList[0][j]) / delta / rList[0][j]);
                        rList_delta[1][j] = Math.Abs((rList_delta[1][j] - rList[1][j]) / delta / rList[1][j]);
                        //rList_delta[0][j] = (rList_delta[0][j] - rList[0][j]) / delta ;
                        //rList_delta[1][j] = (rList_delta[1][j] - rList[1][j]) / delta ;
                        rList_delta[2][j] = Math.Abs((rList_delta[2][j] - rList[2][j]) / delta / rList[2][j]);
                    }

                    listL.Add(rList_delta[0]);
                    listR.Add(rList_delta[1]);
                    listB.Add(rList_delta[2]);
                    legendsL.Add($"delta = {delta}") ;
                    legendsR.Add($"delta = {delta}");
                    legendsB.Add($"delta = {delta}");
                    delta *= 0.1;
                }
                string Title1 = "Re " + fieldCompParams.FC + " производная по " +  fieldCompParams.drvParam;
                string Title2 = "Im " + fieldCompParams.FC + " производная по " + fieldCompParams.drvParam;
                string Title3 = "Mod " + fieldCompParams.FC + " производная по " + fieldCompParams.drvParam;
                ChartView.DrawChart(chart1, Title1, "ro в длинах волн в воздухе", "", ro, listL, legendsL);
                ChartView.DrawChart(chart2, Title2, "ro в длинах волн в воздухе", "", ro, listR, legendsR);
                ChartView.DrawChart(chart3, Title3, "ro в длинах волн в воздухе", "", ro, listB, legendsB);

                IsUpdated = true;
            }
            catch (Exception ex)
            {
                MessageBox.Show("Error in Calculate_Flam()" + ex.ToString());
            }
        }
        public void CalculateFields(Chart chart1, Chart chart2, Chart chart3)
        {
           // MessageBox.Show("CalculateFields");
            try
            {
                List<double[]> listL = new List<double[]>();
                List<double[]> listR = new List<double[]>();
                List<double[]> listB = new List<double[]>();

                List<string> legendsL = new List<string>();
                List<string> legendsR = new List<string>();
                List<string> legendsB = new List<string>();

                string Title1 = "| d" + fieldCompParams.FC.ToString() + "/dh " + "| (norm)";
                string Title2 = "| d" + fieldCompParams.FC.ToString() + "/dh " + "| (norm)";
                string Title3 = "Mod " + fieldCompParams.FC;


                fieldCompParams.roMax = 0.001;
                fieldCompParams.roMin = 0.1;
                double delta = 1;
                
                int nro = fieldCompParams.nro;
                double[] ro = new double[nro];  // в длинах волн в воздухе
                double hro = (fieldCompParams.roMax - fieldCompParams.roMin) / (nro - 1);
                double[] ros = new double[nro];

                for (int j = 0; j < nro; j++)
                {
                    ro[j] = fieldCompParams.roMin + j * hro;
                    ros[j] = ro[j] * modelParams.wave; // в метрах для вычислений
                }
              
                string legend_string = "";
                double eps1 = 1;
                double eps2 = 1;
                double h = 1;

                for (int jdelta = 0; jdelta < 3; jdelta++)
                {
                    if (fieldCompParams.drvParam == DrvParam.eps1)
                    {
                        eps1 = modelParams.eps1 * delta;
                        legend_string = $"eps1 = {eps1}";
                    }
                    else if (fieldCompParams.drvParam == DrvParam.eps2)
                    {
                        eps2 = modelParams.eps2 * delta;
                        legend_string = $"eps2 = {eps2}";
                    }
                    else
                    {
                        h = modelParams.h * delta;
                        legend_string = $"h = {h}";
                    }
                 //   MessageBox.Show($"h = {h}");
                    ModelParams modelParams_delta =
                       new ModelParams(eps1, eps2, h,
                                       modelParams.ro, modelParams.z0, modelParams.z, modelParams._omega);

                    DipolesFields dipolesFields_delta = new DipolesFields(modelParams_delta, dipolesFields.FieldComponent);

                    List<double[]> rList_delta = new List<double[]>(4);
                    rList_delta.Add(new double[nro]);
                    rList_delta.Add(new double[nro]);
                    rList_delta.Add(new double[nro]);
                    rList_delta.Add(new double[nro]);
                    dipolesFields_delta.CalculateFields_byAnderson(ros, rList_delta);

                    listL.Add(rList_delta[0]);
                    legendsL.Add(legend_string);
                    delta *=2;
                }
                ChartView.DrawChart(chart1, Title1, "ro в длинах волн в воздухе", "", ro, listL, legendsL);













                fieldCompParams.roMax = 0.1;
                fieldCompParams.roMin = 1.5;
                delta = 1;

                nro = fieldCompParams.nro;
                ro = new double[nro];  // в длинах волн в воздухе
                hro = (fieldCompParams.roMax - fieldCompParams.roMin) / (nro - 1);
                ros = new double[nro];

                for (int j = 0; j < nro; j++)
                {
                    ro[j] = fieldCompParams.roMin + j * hro;
                    ros[j] = ro[j] * modelParams.wave; // в метрах для вычислений
                }

                legend_string = "";
                eps1 = 1;
                eps2 = 1;
                h = 1;

                for (int jdelta = 0; jdelta < 3; jdelta++)
                {
                    if (fieldCompParams.drvParam == DrvParam.eps1)
                    {
                        eps1 = modelParams.eps1 * delta;
                        legend_string = $"eps1 = {eps1}";
                    }
                    else if (fieldCompParams.drvParam == DrvParam.eps2)
                    {
                        eps2 = modelParams.eps2 * delta;
                        legend_string = $"eps2 = {eps2}";
                    }
                    else
                    {
                        h = modelParams.h * delta;
                        legend_string = $"h = {h}";
                    }
                    //   MessageBox.Show($"h = {h}");
                    ModelParams modelParams_delta =
                       new ModelParams(eps1, eps2, h,
                                       modelParams.ro, modelParams.z0, modelParams.z, modelParams._omega);

                    DipolesFields dipolesFields_delta = new DipolesFields(modelParams_delta, dipolesFields.FieldComponent);

                    List<double[]> rList_delta = new List<double[]>(4);
                    rList_delta.Add(new double[nro]);
                    rList_delta.Add(new double[nro]);
                    rList_delta.Add(new double[nro]);
                    rList_delta.Add(new double[nro]);
                    dipolesFields_delta.CalculateFields_byAnderson(ros, rList_delta);

                    listR.Add(rList_delta[0]);
                    legendsR.Add(legend_string);
                    delta *= 2;
                }
                ChartView.DrawChart(chart2, Title2, "ro в длинах волн в воздухе", "", ro, listR, legendsR);

                IsUpdated = true;
            }
            catch (Exception ex)
            {
                MessageBox.Show("Error in Calculate_Flam()" + ex.ToString());
            }
        }
        public void CalculateFields_Compare_Anderson_Spline(Chart chart1, Chart chart2, Chart chart3)
        {
            // MessageBox.Show("CalculateFields");
            try
            {
                List<double[]> listL = new List<double[]>();
                List<double[]> listR = new List<double[]>();
                List<double[]> listB = new List<double[]>();

                List<string> legendsL = new List<string>();
                List<string> legendsR = new List<string>();
                List<string> legendsB = new List<string>();

                string Title1 = fieldCompParams.FC.ToString();
                string Title2 = fieldCompParams.FC.ToString();

                fieldCompParams.roMax = 0.001;
                fieldCompParams.roMin = 0.1;

                int nro = fieldCompParams.nro;
                double[] ro = new double[nro];  // в длинах волн в воздухе
                double hro = (fieldCompParams.roMax - fieldCompParams.roMin) / (nro - 1);
                double[] ros = new double[nro];

                for (int j = 0; j < nro; j++)
                {
                    ro[j] = fieldCompParams.roMin + j * hro;
                    ros[j] = ro[j] * modelParams.wave; // в метрах для вычислений
                }
                /// Anderson
                List<double[]> rListA = new List<double[]>(4);
                rListA.Add(new double[nro]);
                rListA.Add(new double[nro]);
                rListA.Add(new double[nro]);
                rListA.Add(new double[nro]);

                dipolesFields.CalculateFields_byAnderson(ros, rListA);
                string legend_string = "Anderson";

                listL.Add(rListA[0]);
                listR.Add(rListA[1]);
                listB.Add(rListA[2]);
                legendsL.Add(legend_string);
                legendsR.Add(legend_string);
                legendsB.Add(legend_string);


                ChartView.DrawChart(chart1, Title1, "ro в длинах волн в воздухе", "", ro, listL, legendsL);





                fieldCompParams.roMax = 0.1;
                fieldCompParams.roMin = 1;

                nro = fieldCompParams.nro;
                ro = new double[nro];  // в длинах волн в воздухе
                hro = (fieldCompParams.roMax - fieldCompParams.roMin) / (nro - 1);
                ros = new double[nro];

                for (int j = 0; j < nro; j++)
                {
                    ro[j] = fieldCompParams.roMin + j * hro;
                    ros[j] = ro[j] * modelParams.wave; // в метрах для вычислений
                }
                /// Anderson
                rListA = new List<double[]>(4);
                rListA.Add(new double[nro]);
                rListA.Add(new double[nro]);
                rListA.Add(new double[nro]);
                rListA.Add(new double[nro]);

                dipolesFields.CalculateFields_byAnderson(ros, rListA);
                legend_string = "Anderson";

                listL.Add(rListA[0]);
                listR.Add(rListA[1]);
                listB.Add(rListA[2]);
                legendsL.Add(legend_string);
                legendsR.Add(legend_string);
                legendsB.Add(legend_string);

                ChartView.DrawChart(chart1, Title2, "ro в длинах волн в воздухе", "", ro, listL, legendsL);

                IsUpdated = true;
            }
            catch (Exception ex)
            {
                MessageBox.Show("Error in Calculate_Flam()" + ex.ToString());
            }
        }
        public void TestIntegral()
        {
            try
            {
                //double ro = 5;
                //double res = 0;
                //double exact = 0;
                //int n = 0;
               
                //dipolesFields.TestIntegral_J0(ro, ref res, ref exact, ref n);
                //MessageBox.Show($"ro = {ro} res = {res} exact = {exact} n = {n}");
            }
            catch (Exception ex)
            {
                MessageBox.Show("Error in TestIntegral()" + ex.ToString());
            }
        }
        public override string ToString()
        {
            return lamParams.ToString() + "\n" + modelParams.ToString() + fieldCompParams.ToString(); ;
        }
    }
}
