using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using Microsoft.Win32;
using Lib_CSharp;
using System.Windows.Forms.DataVisualization.Charting;


namespace Wpf_Dipoles
{
    /// <summary>
    /// Логика взаимодействия для MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        Chart chartL = new Chart();
        Chart chartR = new Chart();
        ViewData viewData = new ViewData();
        public MainWindow()
        {
            InitializeComponent();
            winFormsHost_L.Child = chartL;
            winFormsHost_R.Child = chartR;
            DataContext = viewData;
            grid_model.DataContext = viewData.modelParams;
            grid_lam.DataContext = viewData.lamParams;
            grid_fldComp.DataContext = viewData.fieldCompParams;
        }
        private void OpenCommandHandler(object sender, ExecutedRoutedEventArgs e)
        {
            OpenFileDialog ofd = new OpenFileDialog();
            if (ofd.ShowDialog() == true)
            {
            }
        }

        private void CanOpenCommandHandler(object sender, CanExecuteRoutedEventArgs e)
        {

        }

        private void SaveCommandHandler(object sender, ExecutedRoutedEventArgs e)
        {
            SaveFileDialog ofd = new SaveFileDialog();
            if (ofd.ShowDialog() == true)
            {
                try
                {
                    //viewData.splinesData.Save(ofd.FileName);
                    //viewData.IsUpdated = true;
                }
                catch (Exception ex)
                {
                    MessageBox.Show(ex.ToString());
                }
            }
        }

        private void CanSaveCommandHandler(object sender, CanExecuteRoutedEventArgs e)
        {
            if (viewData.IsUpdated) e.CanExecute = true;
            else e.CanExecute = false;
        }
        private void Window_Loaded(object sender, RoutedEventArgs e)
        {
            comboBox_Func.ItemsSource = Enum.GetValues(typeof(DpF));
            comboBox_Field.ItemsSource = Enum.GetValues(typeof(FldComp));
        }

        private void CanFieldsCommandHandler(object sender, CanExecuteRoutedEventArgs e)
        {
            if (Validation.GetHasError(textBox_nro) == true || Validation.GetHasError(textBox_roMin) == true || Validation.GetHasError(textBox_roMax) == true)
            //    || Validation.GetHasError(textBox_ro) == true || Validation.GetHasError(textBox_z0) == true || Validation.GetHasError(textBox_z) == true
            //    || Validation.GetHasError(textBox_eps1) == true || Validation.GetHasError(textBox_eps2) == true || Validation.GetHasError(textBox_h) == true)
            {
                e.CanExecute = false;
            }
            else e.CanExecute = true;
        }

        private void FlamCommandHandler(object sender, ExecutedRoutedEventArgs e)
        {
            try
            {
                viewData.CalculateDraw_Flam(chartL, chartR);
            }
            catch (Exception ex)
            {
                MessageBox.Show(ex.ToString());
            }
        }
        private void CanFlamCommandHandler(object sender, CanExecuteRoutedEventArgs e)
        {
            if (Validation.GetHasError(textBox_n) == true || Validation.GetHasError(textBox_lamL) == true || Validation.GetHasError(textBox_lamR) == true 
                || Validation.GetHasError(textBox_ro) == true || Validation.GetHasError(textBox_z0) == true || Validation.GetHasError(textBox_z) == true
                || Validation.GetHasError(textBox_eps1) == true || Validation.GetHasError(textBox_eps2) == true || Validation.GetHasError(textBox_h) == true)
            {
                e.CanExecute = false;
            }
            else e.CanExecute = true;
        }

        private void FieldsCommandHandler(object sender, ExecutedRoutedEventArgs e)
        {
            try
            {
                viewData.CalculateFields(chartL, chartR);

                //ChartView chartViewL = new ChartView(viewData.dataForOutputL);
                //chartViewL.DrawChart(chartL);
                //listBox_InfoDataL.DataContext = viewData.dataForOutputL;

                //ChartView chartViewR = new ChartView(viewData.dataForOutputR);
                //chartViewR.DrawChart(chartR);
                //listBox_InfoDataR.DataContext = viewData.dataForOutputR;
            }
            catch (Exception ex)
            {
                MessageBox.Show(ex.ToString());
            }
        }

        private void CanDrawCommandHandler(object sender, CanExecuteRoutedEventArgs e)
        {

        }

        private void DrawCommandHandler(object sender, ExecutedRoutedEventArgs e)
        {

        }

        private void radiobutton_0k0_Click(object sender, RoutedEventArgs e)
        {
            viewData.lamParams.lamL = 0;
            viewData.lamParams.lamR = viewData.modelParams.k0;
            //double k0 = viewData.modelParams.k0;
            //double k0_om = viewData.modelParams.omega * Math.Sqrt(EMC.eps0 * EMC.mu);
            //MessageBox.Show($"k0 = {k0} \nk0_om = {k0_om}");
        }

        private void radiobutton_k0k1_Click(object sender, RoutedEventArgs e)
        {
            viewData.lamParams.lamL = viewData.modelParams.k0;
            viewData.lamParams.lamR = viewData.modelParams.k1;
            //double k0 = viewData.modelParams.k0;
            //double k1 = viewData.modelParams.k1;
            //double k0_om = viewData.modelParams.omega * Math.Sqrt(EMC.eps0 * EMC.mu);
            //double k1_om = viewData.modelParams.omega * Math.Sqrt(viewData.modelParams.eps1 * EMC.eps0 * EMC.mu);
            //MessageBox.Show($"k0 = {k0} \nk0_om = {k0_om}" + $"\nk1 = {k1} \nk1_om = {k1_om}");
        }

        private void radiobutton_k1k2_Click(object sender, RoutedEventArgs e)
        {
            viewData.lamParams.lamL = viewData.modelParams.k1;
            viewData.lamParams.lamR = viewData.modelParams.k2;
         
            //double k1 = viewData.modelParams.k1;
            //double k2 = viewData.modelParams.k2;
            //double k1_om = viewData.modelParams.omega * Math.Sqrt(viewData.modelParams.eps1 * EMC.eps0 * EMC.mu);
            //double k2_om = viewData.modelParams.omega * Math.Sqrt(viewData.modelParams.eps2 * EMC.eps0 * EMC.mu);
            //MessageBox.Show($"k1 = {k1} \nk1_om = {k1_om}" + $"\nk2 = {k2} \nk2_om = {k2_om}");
        }

        private void radiobutton_k0near_Click(object sender, RoutedEventArgs e)
        {
            viewData.lamParams.lamL = 0.99 * viewData.modelParams.k0;
            viewData.lamParams.lamR = 1.01 * viewData.modelParams.k0;
        }

        private void radiobutton_k1near_Click(object sender, RoutedEventArgs e)
        {
            viewData.lamParams.lamL = 0.99 * viewData.modelParams.k1;
            viewData.lamParams.lamR = 1.01 * viewData.modelParams.k1;
        }
        private void radiobutton_k2near_Click(object sender, RoutedEventArgs e)
        {
            viewData.lamParams.lamL = 0.99 * viewData.modelParams.k2;
            viewData.lamParams.lamR = 1.01 * viewData.modelParams.k2;
        }

    }
}


