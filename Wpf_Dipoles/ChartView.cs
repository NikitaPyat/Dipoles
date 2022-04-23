using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms.DataVisualization.Charting;
using System.Windows;
using Lib_CSharp;

namespace Wpf_Dipoles
{
    static class ChartView
    {
        static System.Drawing.Color[] Colors =
                 new System.Drawing.Color[] { System.Drawing.Color.Red,
                                            System.Drawing.Color.Blue,
                                            System.Drawing.Color.Cyan,
                                            System.Drawing.Color.Magenta,
                                            System.Drawing.Color.Green,
                                            System.Drawing.Color.Brown,
                                            System.Drawing.Color.LightBlue};


        public static void DrawChart(Chart chart, string Title, string xTitle, string yTitle, 
                              double[] X, List<double[]> YList, List<string> legends)
        {
            chart.ChartAreas.Clear();
            chart.Series.Clear();
            chart.Legends.Clear();
            chart.Titles.Clear();
            Legend legend = new Legend("Custom");


            chart.ChartAreas.Add(new ChartArea("chartArea1"));
            chart.ChartAreas[0].AxisX.MajorGrid.LineColor = System.Drawing.Color.LightGray;
            chart.ChartAreas[0].AxisY.MajorGrid.LineColor = System.Drawing.Color.LightGray;
            chart.ChartAreas[0].AxisX.IsMarginVisible = false;
            legend.InsideChartArea = chart.ChartAreas[0].Name;

            chart.Titles.Add(new Title(Title));
            chart.ChartAreas[0].AxisX.Title = xTitle;
            chart.ChartAreas[0].AxisY.Title = yTitle;

            //chart.ChartAreas[0].AxisX.LabelStyle.Format = "F2";
            //chart.ChartAreas[0].AxisY.LabelStyle.Format = "F2";

            for (int js = 0; js < YList.Count; js++)
            {
                System.Drawing.Color jsColor = Colors[js % Colors.Length];
                chart.Series.Add(js.ToString());
                chart.Series[js].Points.DataBindXY(X, YList[js]);

                //chart.Series[js].MarkerStyle = MarkerStyle.Circle;
                //chart.Series[js].MarkerSize = 5;
                //chart.Series[js].MarkerColor = jsColor;

                chart.Series[js].ChartType = SeriesChartType.Spline;
                chart.Series[js].Color = jsColor;

                chart.Series[js].IsVisibleInLegend = false;
                LegendItem legendItem = new LegendItem();
                legendItem.Name = legends[js];
                legendItem.Color = jsColor;
               //legendItem.ImageStyle = LegendImageStyle.Line;
                legend.CustomItems.Add(legendItem);
            }
            chart.Legends.Add(legend);
        }
    }
}
