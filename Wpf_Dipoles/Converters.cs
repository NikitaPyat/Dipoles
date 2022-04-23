using System;
using System.Collections.Generic;
using System.Text;
using System.Numerics;
using System.Windows.Data;
using System.Windows;
using Lib_CSharp;

namespace Wpf_Dipoles
{
    class dblConverter : IValueConverter
    {
        public object Convert(object value, Type targetType, object parameter, System.Globalization.CultureInfo culture)
        {
            double item = (double)value;
            return item.ToString("F3");
        }

        public object ConvertBack(object value, Type targetType, object parameter, System.Globalization.CultureInfo culture)
        {
            return value;
        }
    }

    //class measuredDataConverter : IValueConverter
    //{
    //    public object Convert(object value, Type targetType, object parameter, System.Globalization.CultureInfo culture)
    //    {
    //        LambdaParams item = value as LambdaParams;
    //        if (item != null)
    //            return $"n = {item.n} xL = {item.xL} xR = {item.xR} F = {item.F}";
    //        else return "";
    //    }

    //    public object ConvertBack(object value, Type targetType, object parameter, System.Globalization.CultureInfo culture)
    //    {
    //        return value;
    //    }
    //}

    //class splineParametersConverter : IValueConverter
    //{
    //    public object Convert(object value, Type targetType, object parameter, System.Globalization.CultureInfo culture)
    //    {
    //        SplineParameters item = value as SplineParameters;
    //        return $"nUni = {item.nUni} Left = {item.DL1} Right = {item.DR1}";
    //    }

    //    public object ConvertBack(object value, Type targetType, object parameter, System.Globalization.CultureInfo culture)
    //    {
    //        return value;
    //    }
    //}
}
