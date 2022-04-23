using System;
using System.Collections.Generic;
using System.Windows.Input;

namespace Wpf_Dipoles
{
    public static class CustomCommands
    {
        public static RoutedCommand DrawCommand =
            new RoutedCommand("DrawCommand", typeof(Wpf_Dipoles.CustomCommands));
        //public static RoutedCommand AddCustomCommand =
        //    new RoutedCommand("AddCustomCommand", typeof(Wpf_Lab2.CustomCommands));
        public static RoutedCommand ModelCommand =
            new RoutedCommand("ModelCommand", typeof(Wpf_Dipoles.CustomCommands));
      
        public static RoutedCommand FieldsCommand =
        new RoutedCommand("FieldsCommand", typeof(Wpf_Dipoles.CustomCommands));
        public static RoutedCommand FlamCommand =
        new RoutedCommand("FlamCommand", typeof(Wpf_Dipoles.CustomCommands));
    }
}
