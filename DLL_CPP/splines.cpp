#include "pch.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include "mkl.h"
using namespace std;

extern "C"  _declspec(dllexport)
void Wrap_CubicNatural_Uniform(int nX, double xL, double xR, int nF, double* fValues, 
                               int nsite, double* site, int ndorder, double* res, 
                               int nintegral, double* llimit, double* rlimit, double* res_integral,
                               int& return_value)
{
    
    ofstream ofs ("cpp_errors.txt", ios_base::out);
    //ofs.open(); // окрываем файл для записи
    ofs << "===== nX = " << nX << endl;
    ofs << "===== xL = " << xL << "  xR = " << xR << endl;
    ofs << "===== nF = " << nF << endl;
    for (int j = 0; j < nF * nX; j++)  ofs << "fValues[" << j << "] = " << fValues[j] << endl;
    ofs << "===== nsite = " << nsite << endl;
    for (int j = 0; j < nsite; j++)  ofs << "site[" << j << "] = " << site[j] << endl;
    ofs << "===== ndorder = " << ndorder << endl;
    ofs << "===== nintegral = " << nintegral << endl;
    for (int j = 0; j < nintegral; j++)  ofs << "llimit[" << j << "] = " << llimit[j] << endl;
    for (int j = 0; j < nintegral; j++)  ofs << "rlimit[" << j << "] = " << rlimit[j] << endl;
    ofs.close();     // закрываем файл
   
    //double y[nx * ny];                   // значения функции (function values)

    //double left = xL;
    //double right = xR;
    // double freq;
    // double xx[N];
  /*  double left_val[N - 1], right_val[N - 1];
    double left_der1[N - 1], right_der1[N - 1];
    double left_der2[N - 1], right_der2[N - 1];*/

    //int i, j, errcode = 0;
    //int errnums = 0;

    /***** Initializing parameters for Data Fitting task *****/
    MKL_INT sorder = DF_PP_CUBIC;       // порядок сплайна (spline order)- 
    MKL_INT stype = DF_PP_NATURAL;      // тип сплайна (spline type)
  

    /***** Parameters describing interpolation interval *****/
    MKL_INT nx = nX;                        // число узлов сплайна (number of break points)
    MKL_INT xhint = DF_UNIFORM_PARTITION;   // additional info about break points

    /* Limits of interpolation interval are provided in case of uniform grid */
    double xUniform[2];                        // границы интервала интерполирования (limits of the interpolation interval)
    xUniform[0] = xL;
    xUniform[1] = xR;

    /***** Parameters describing function *****/
    MKL_INT ny = nF;                     // размер вектор-функции (number of functions)
    MKL_INT yhint = DF_NO_HINT;          // additional info about function
 
    /***** Parameters describing spline coefficients storage *****/
    MKL_INT nscoeff = ny * (nx - 1) * DF_PP_CUBIC;   // число коэффициентов сплайна (number of spline coefficients)
    MKL_INT scoeffhint = DF_NO_HINT;                 // additional info about spline coefficients
    double * scoeff = new double[nscoeff];             // массив коэффициентов сплайна (array of spline coefficients)

    /***** Parameters describing boundary conditions type *****/
    MKL_INT bc_type = DF_BC_FREE_END;   // тип граничных условий (boundary conditions type)
    MKL_INT nbc;                        // число граничных условий (number of boundary conditions)

    /* No additional parameters are provided for Free-End boundary conditions */
    double* bc;                         // boundary conditions
    bc = 0;
    /***** Parameters describing internal conditions type *****/
    /* No internal conditions are provided for natural cubic spline */
    MKL_INT ic_type = DF_NO_IC;          // тип условий во внутренних точках (internal conditions type)
    double* ic;                         // условия во внутренних точках (internal conditions)
    ic = 0; 
    
        /***** Create Data Fitting task *****/
        DFTaskPtr task;                     // Data Fitting task descriptor
        int errcode = dfdNewTask1D(&task, nx, xUniform, xhint, ny, fValues, yhint);

        printf("\ndfdNewTask1D errcode = %d:\n", errcode); //CheckDfError(errcode);

        /***** Edit task parameters for natural cubic spline construction *****/
       // bc_type = DF_BC_1ST_RIGHT_DER;
       // bc = new double[1];
       //// bc[0] = 0;
       // bc[0] = der_1_at_right_end;

        errcode = dfdEditPPSpline1D(task, sorder, stype, bc_type, bc, ic_type, ic, scoeff, scoeffhint);

       // printf("\ndfdEditPPSpline1D errcode = %d:\n", errcode);  //CheckDfError(errcode);

        /***** Construct natural cubic spline using STD method *****/
        errcode = dfdConstruct1D(task, DF_PP_SPLINE, DF_METHOD_STD);
       // printf("\ndfdConstruct1D errcode = %d:\n", errcode); // CheckDfError(errcode);

        MKL_INT _nsite = nsite;
        // const MKL_INT sitehint = DF_SORTED_DATA;

         //MKL_INT *dorder;
         /*if (ndorder == 1) dorder = new MKL_INT[1]{ 1 };
         else if (ndorder == 2) dorder = new MKL_INT[2]{ 1 , 1};
         else if (ndorder == 3) dorder = new MKL_INT[3]{ 1 , 1 , 1};
         else if (ndorder == 4) dorder = new MKL_INT[4]{ 1 , 1 , 1 , 1};
         else*/
        MKL_INT dorder[5]{ 1 , 1 , 1 , 1 , 1 };
        // вычисляются значения функции и первой и второй производных
       // const double* datahint = NULL;
       // const MKL_INT rhint = DF_MATRIX_STORAGE_ROWS;
       // MKL_INT* cell = NULL;

        // Вычисление значений сплайна и двух производных в точках site
        errcode = dfdInterpolate1D(task, DF_INTERP, DF_METHOD_PP,
            _nsite, site, DF_SORTED_DATA,
            ndorder, dorder, NULL,
            res, DF_MATRIX_STORAGE_ROWS, NULL);

        /***** Parameters decsribing integration limits *****/
        //MKL_INT nlim = 1;                       // число пар пределов интегрирования (number of pairs of integration limits)
        //MKL_INT llimhint = DF_NON_UNIFORM_PARTITION;  // дополнительная информация о структуре левых пределов интегрирования
        //                                              //additional info about the structure
        //                                              // of left integration limits
        //MKL_INT rlimhint = DF_NON_UNIFORM_PARTITION;  // additional info about the structure
                                                      // of right integration limits

        //double* llim = new double[nlim];            // массив левых концов отрезка интегрирования (left  integration limits)
        //double* rlim = new double[nlim];            //  массив правых концов отрезка интегрирования (right integration limits)

        /***** Additional information about the structure of integration limits *****/
        /* No additional info is provided */
        double* ldatahint, * rdatahint;     // additional info about the structure
                                            // of x and integration limits
        ldatahint = 0;
        rdatahint = 0;

        /***** Parameter dascribing integration results storage format *****/
       // MKL_INT rhint = DF_NO_HINT;        // формат хранения результатов интегрирования(integration results storage format)

        //double* r = new double[ny * nlim];                // массив значений интегралов (integration results)
        //double* ref_r = new double[ny * nlim];            // reference integration results
        MKL_INT nlim = nintegral;
     //===== Вычисление интегралов
        errcode = dfdIntegrate1D(task, DF_METHOD_PP,
                                 nlim, llimit, DF_NON_UNIFORM_PARTITION, rlimit, DF_NON_UNIFORM_PARTITION,
                                 ldatahint, rdatahint, res_integral, DF_MATRIX_STORAGE_ROWS);
        ////CheckDfError(errcode);
        //printf("\ndfdIntegrate1D errcode = %d:\n", errcode);

        ///***** Check computed coefficients *****/
        //errcode = dUniformData(xx, left, right, nx);
        //CheckDfError(errcode);

        ///***** Check spline values in break points *****/
        //errcode = dCheckCubBreakPoints(nx, xx, ny, y, scoeff,
        //    left_val, right_val);
        //if (errcode < 0) errnums++;

        ///***** Check that spline 1st derivatives are equal for left
        //       and right piece of the spline for each break point *****/
        //errcode = dCheckCub1stDerConsistency(nx, xx, ny, scoeff,
        //    left_der1, right_der1);
        //if (errcode < 0) errnums++;

        ///***** Check that spline 2nd derivatives are equal for left
        //       and right piece of the spline for each break point *****/
        //errcode = dCheckCub2ndDerConsistency(nx, xx, ny, scoeff,
        //    left_der2, right_der2);
        //if (errcode < 0) errnums++;

        ///***** Check boundary conditions *****/
        //errcode = dCheckCubBC(nx, xx, ny, scoeff, bc_type, bc);
        //if (errcode < 0) errnums++;

        ///***** Check integration results *****/
        //errcode = dCheckCubIntegrRes(nx, xx, ny, scoeff, nlim, llim, rlim, r, ref_r);
        //if (errcode < 0) errnums++;

        ///***** Print results *****/
        //printf("Number of break points : %d\n", (int)nx);

        ///***** Print given function *****/
        //printf("\n  X           Y\n");

     /*   for (int j = 0; j < nx; j++)
        {
            printf(" %+lf   %+lf\n", xx[j], y[j]);
        }*/

       /***** Delete Data Fitting task *****/
        errcode = dfDeleteTask(&task);

        printf("\ndfDeleteTask errcode = %d:\n", errcode); //CheckDfError(errcode);

                                                       /***** Print summary of the test *****/
        if (errcode != 0)
        {
            printf("\n\nError: Computed natural cubic spline coefficients");
            printf(" or integrals are incorrect\n");
            return_value = 1;
        }
        else
        {
            printf("\n\nComputed natural cubic spline coefficients");
            printf(" and integrals are correct\n");
            return_value = 0;
        }
        delete[] scoeff;
        /* delete[] llim;
         delete[] rlim;
         delete[] r;
         delete[] ref_r;*/
    }

    // Кубический сплайн только для одной функции (векторная функция с одним элементом)
    extern "C"  _declspec(dllexport)
        void Wrap_CubicNatural_1Function(int nX, double xL, double xR, double* fValues,
            int nsite, double* site, int ndorder, double* res,
            int& return_value)
    {

        ofstream ofs("cpp_errors.txt", ios_base::out);
        //ofs.open(); // окрываем файл для записи
        ofs << "===== nX = " << nX << endl;
        ofs << "===== xL = " << xL << "  xR = " << xR << endl;
       
        double h = (xR - xL) / (nX - 1);
        for (int j = 0; j < nX; j++)  ofs << "x[" << j << "] = " << (xL + j * h) << "  fValues = " << fValues[j] << endl;
        ofs << "===== nsite = " << nsite << endl;
        for (int j = 0; j < nsite; j++)  ofs << "site[" << j << "] = " << site[j] << endl;
        ofs << "===== ndorder = " << ndorder << endl;
        ofs.close();     // закрываем файл

        //double y[nx * ny];                   // значения функции (function values)

        //double left = xL;
        //double right = xR;
        // double freq;
        // double xx[N];
      /*  double left_val[N - 1], right_val[N - 1];
        double left_der1[N - 1], right_der1[N - 1];
        double left_der2[N - 1], right_der2[N - 1];*/

        //int i, j, errcode = 0;
        //int errnums = 0;

        /***** Initializing parameters for Data Fitting task *****/
        MKL_INT sorder = DF_PP_CUBIC;       // порядок сплайна (spline order)- 
        MKL_INT stype = DF_PP_NATURAL;      // тип сплайна (spline type)


        /***** Parameters describing interpolation interval *****/
        MKL_INT nx = nX;                        // число узлов сплайна (number of break points)
        MKL_INT xhint = DF_UNIFORM_PARTITION;   // additional info about break points
        double xUniform[2];                        // границы интервала интерполирования (limits of the interpolation interval)
        xUniform[0] = xL;
        xUniform[1] = xR;
        /* Limits of interpolation interval are provided in case of uniform grid */
        //double xUniform[2];                        // границы интервала интерполирования (limits of the interpolation interval)
        //xUniform[0] = xL;
        //xUniform[1] = xR;

        /***** Parameters describing function *****/
        MKL_INT ny = 1;                     // размер вектор-функции (number of functions)
        MKL_INT yhint = DF_NO_HINT;          // additional info about function

        /***** Parameters describing spline coefficients storage *****/
        MKL_INT nscoeff = ny * (nx - 1) * DF_PP_CUBIC;   // число коэффициентов сплайна (number of spline coefficients)
        MKL_INT scoeffhint = DF_NO_HINT;                 // additional info about spline coefficients
        double* scoeff = new double[nscoeff];             // массив коэффициентов сплайна (array of spline coefficients)

        /***** Parameters describing boundary conditions type *****/
        MKL_INT bc_type = DF_BC_FREE_END;   // тип граничных условий (boundary conditions type)
        MKL_INT nbc;                        // число граничных условий (number of boundary conditions)

        /* No additional parameters are provided for Free-End boundary conditions */
        double* bc;                         // boundary conditions
        bc = 0;
        /***** Parameters describing internal conditions type *****/
        /* No internal conditions are provided for natural cubic spline */
        MKL_INT ic_type = DF_NO_IC;          // тип условий во внутренних точках (internal conditions type)
        double* ic;                         // условия во внутренних точках (internal conditions)
        ic = 0;

        /***** Create Data Fitting task *****/
        DFTaskPtr task;                     // Data Fitting task descriptor
        int errcode = dfdNewTask1D(&task, nx, xUniform, xhint, ny, fValues, yhint);

        printf("\ndfdNewTask1D errcode = %d:\n", errcode); //CheckDfError(errcode);

        /***** Edit task parameters for natural cubic spline construction *****/
       // bc_type = DF_BC_1ST_RIGHT_DER;
       // bc = new double[1];
       //// bc[0] = 0;
       // bc[0] = der_1_at_right_end;

        errcode = dfdEditPPSpline1D(task, sorder, stype, bc_type, bc, ic_type, ic, scoeff, scoeffhint);

        // printf("\ndfdEditPPSpline1D errcode = %d:\n", errcode);  //CheckDfError(errcode);

         /***** Construct natural cubic spline using STD method *****/
        errcode = dfdConstruct1D(task, DF_PP_SPLINE, DF_METHOD_STD);
        // printf("\ndfdConstruct1D errcode = %d:\n", errcode); // CheckDfError(errcode);

        MKL_INT _nsite = nsite;
        // const MKL_INT sitehint = DF_SORTED_DATA;

         //MKL_INT *dorder;
         /*if (ndorder == 1) dorder = new MKL_INT[1]{ 1 };
         else if (ndorder == 2) dorder = new MKL_INT[2]{ 1 , 1};
         else if (ndorder == 3) dorder = new MKL_INT[3]{ 1 , 1 , 1};
         else if (ndorder == 4) dorder = new MKL_INT[4]{ 1 , 1 , 1 , 1};
         else*/
        MKL_INT dorder[5]{ 1 , 1 , 1 , 1 , 1 };
        // вычисляются значения функции и первой и второй производных
       // const double* datahint = NULL;
       // const MKL_INT rhint = DF_MATRIX_STORAGE_ROWS;
       // MKL_INT* cell = NULL;

        // Вычисление значений сплайна и двух производных в точках site
        errcode = dfdInterpolate1D(task, DF_INTERP, DF_METHOD_PP,
            _nsite, site, DF_SORTED_DATA,
            ndorder, dorder, NULL,
            res, DF_MATRIX_STORAGE_ROWS, NULL);

        /***** Parameters decsribing integration limits *****/
        //MKL_INT nlim = 1;                       // число пар пределов интегрирования (number of pairs of integration limits)
        //MKL_INT llimhint = DF_NON_UNIFORM_PARTITION;  // дополнительная информация о структуре левых пределов интегрирования
        //                                              //additional info about the structure
        //                                              // of left integration limits
        //MKL_INT rlimhint = DF_NON_UNIFORM_PARTITION;  // additional info about the structure
                                                      // of right integration limits

        //double* llim = new double[nlim];            // массив левых концов отрезка интегрирования (left  integration limits)
        //double* rlim = new double[nlim];            //  массив правых концов отрезка интегрирования (right integration limits)

        /***** Additional information about the structure of integration limits *****/
        /* No additional info is provided */
       // double* ldatahint, * rdatahint;     // additional info about the structure
       //                                     // of x and integration limits
       // ldatahint = 0;
       // rdatahint = 0;

       // /***** Parameter dascribing integration results storage format *****/
       //// MKL_INT rhint = DF_NO_HINT;        // формат хранения результатов интегрирования(integration results storage format)

       // //double* r = new double[ny * nlim];                // массив значений интегралов (integration results)
       // //double* ref_r = new double[ny * nlim];            // reference integration results
       // MKL_INT nlim = nintegral;
       // //===== Вычисление интегралов
       // errcode = dfdIntegrate1D(task, DF_METHOD_PP,
       //     nlim, llimit, DF_NON_UNIFORM_PARTITION, rlimit, DF_NON_UNIFORM_PARTITION,
       //     ldatahint, rdatahint, res_integral, DF_MATRIX_STORAGE_ROWS);
        ////CheckDfError(errcode);
        //printf("\ndfdIntegrate1D errcode = %d:\n", errcode);

        ///***** Check computed coefficients *****/
        //errcode = dUniformData(xx, left, right, nx);
        //CheckDfError(errcode);

        ///***** Check spline values in break points *****/
        //errcode = dCheckCubBreakPoints(nx, xx, ny, y, scoeff,
        //    left_val, right_val);
        //if (errcode < 0) errnums++;

        ///***** Check that spline 1st derivatives are equal for left
        //       and right piece of the spline for each break point *****/
        //errcode = dCheckCub1stDerConsistency(nx, xx, ny, scoeff,
        //    left_der1, right_der1);
        //if (errcode < 0) errnums++;

        ///***** Check that spline 2nd derivatives are equal for left
        //       and right piece of the spline for each break point *****/
        //errcode = dCheckCub2ndDerConsistency(nx, xx, ny, scoeff,
        //    left_der2, right_der2);
        //if (errcode < 0) errnums++;

        ///***** Check boundary conditions *****/
        //errcode = dCheckCubBC(nx, xx, ny, scoeff, bc_type, bc);
        //if (errcode < 0) errnums++;

        ///***** Check integration results *****/
        //errcode = dCheckCubIntegrRes(nx, xx, ny, scoeff, nlim, llim, rlim, r, ref_r);
        //if (errcode < 0) errnums++;

        ///***** Print results *****/
        //printf("Number of break points : %d\n", (int)nx);

        ///***** Print given function *****/
        //printf("\n  X           Y\n");

     /*   for (int j = 0; j < nx; j++)
        {
            printf(" %+lf   %+lf\n", xx[j], y[j]);
        }*/

        /***** Delete Data Fitting task *****/
        errcode = dfDeleteTask(&task);

        printf("\ndfDeleteTask errcode = %d:\n", errcode); //CheckDfError(errcode);

                                                       /***** Print summary of the test *****/
        if (errcode != 0)
        {
            printf("\n\nError: Computed natural cubic spline coefficients");
            printf(" or integrals are incorrect\n");
            return_value = 1;
        }
        else
        {
            printf("\n\nComputed natural cubic spline coefficients");
            printf(" and integrals are correct\n");
            return_value = 0;
        }
        delete[] scoeff;
        /* delete[] llim;
         delete[] rlim;
         delete[] r;
         delete[] ref_r;*/
    }
    //////////////////////////////////

// Кубический сплайн только для одной функции (векторная функция с одним элементом)
// Функция и две производные
extern "C"  _declspec(dllexport)

void Wrap_Splines_Lab2_2022(int n, double* x, double* mData, double condLeft, double condRight,
                            int nN, double* xUni, double* res, int& return_value)
{
                 
    ofstream ofs("cpp_errors.txt", ios_base::out);
    //ofs.open(); // окрываем файл для записи
    ofs << "===== n = " << n << "  nN = " << nN <<endl;
    ofs << "===== condLeft = " << condLeft << "  condRight = " << condRight << endl;

    for (int j = 0; j < n; j++)  ofs << "x[" << j << "] = " << x[j] << "  mData = " << mData[j] << endl;
    ofs << "===== nN = " << nN << endl;
    for (int j = 0; j < nN; j++)  ofs << "xUni[" << j << "] = " << xUni[j] << endl;
  /*  ofs << "===== ndorder = " << ndorder << endl;*/
    ofs.close();     // закрываем файл

    //double y[nx * ny];                   // значения функции (function values)

    //double left = xL;
    //double right = xR;
    // double freq;
    // double xx[N];
  /*  double left_val[N - 1], right_val[N - 1];
    double left_der1[N - 1], right_der1[N - 1];
    double left_der2[N - 1], right_der2[N - 1];*/

    //int i, j, errcode = 0;
    //int errnums = 0;

    /***** Initializing parameters for Data Fitting task *****/
    MKL_INT sorder = DF_PP_CUBIC;       // порядок сплайна (spline order)- 
    MKL_INT stype = DF_PP_NATURAL;      // тип сплайна (spline type)


    /***** Parameters describing interpolation interval *****/
    MKL_INT nx = n;                        // число узлов сплайна (number of break points)
    MKL_INT xhint = DF_NON_UNIFORM_PARTITION;   // additional info about break points
    //double xUniform[2];                         // границы интервала интерполирования (limits of the interpolation interval)
    //xUniform[0] = xL;
    //xUniform[1] = xR;
    /* Limits of interpolation interval are provided in case of uniform grid */
    //double xUniform[2];                        // границы интервала интерполирования (limits of the interpolation interval)
    //xUniform[0] = xL;
    //xUniform[1] = xR;

    /***** Parameters describing function *****/
    MKL_INT ny = 1;                     // размер вектор-функции (number of functions)
    MKL_INT yhint = DF_NO_HINT;          // additional info about function

    /***** Parameters describing spline coefficients storage *****/
    MKL_INT nscoeff = ny * (nx - 1) * DF_PP_CUBIC;   // число коэффициентов сплайна (number of spline coefficients)
    MKL_INT scoeffhint = DF_NO_HINT;                 // additional info about spline coefficients
    double* scoeff = new double[nscoeff];             // массив коэффициентов сплайна (array of spline coefficients)

    /***** Parameters describing boundary conditions type *****/
    MKL_INT bc_type = DF_BC_FREE_END;   // тип граничных условий (boundary conditions type)
    MKL_INT nbc;                        // число граничных условий (number of boundary conditions)

    /* No additional parameters are provided for Free-End boundary conditions */
    double* bc = new double[2]{condLeft, condRight};                         // boundary conditions
  /*  bc = 0;*/
  

  /***** Parameters describing internal conditions type *****/
    /* No internal conditions are provided for natural cubic spline */
    MKL_INT ic_type = DF_NO_IC;          // тип условий во внутренних точках (internal conditions type)
    double* ic;                         // условия во внутренних точках (internal conditions)
    ic = 0;

    /***** Create Data Fitting task *****/
    DFTaskPtr task;                     // Data Fitting task descriptor
    int errcode = dfdNewTask1D(&task, nx, x, xhint, ny, mData, yhint);

    printf("\ndfdNewTask1D errcode = %d:\n", errcode); //CheckDfError(errcode);

    /***** Edit task parameters for natural cubic spline construction *****/
   // bc_type = DF_BC_1ST_RIGHT_DER;
   // bc = new double[1];
   //// bc[0] = 0;
   // bc[0] = der_1_at_right_end;

    errcode = dfdEditPPSpline1D(task, sorder, stype, bc_type, bc, ic_type, ic, scoeff, scoeffhint);

    // printf("\ndfdEditPPSpline1D errcode = %d:\n", errcode);  //CheckDfError(errcode);

     /***** Construct natural cubic spline using STD method *****/
    errcode = dfdConstruct1D(task, DF_PP_SPLINE, DF_METHOD_STD);
    // printf("\ndfdConstruct1D errcode = %d:\n", errcode); // CheckDfError(errcode);

    MKL_INT nsite = nN;
    // const MKL_INT sitehint = DF_SORTED_DATA;

     //MKL_INT *dorder;
     /*if (ndorder == 1) dorder = new MKL_INT[1]{ 1 };
     else if (ndorder == 2) dorder = new MKL_INT[2]{ 1 , 1};
     else if (ndorder == 3) dorder = new MKL_INT[3]{ 1 , 1 , 1};
     else if (ndorder == 4) dorder = new MKL_INT[4]{ 1 , 1 , 1 , 1};
     else*/
    MKL_INT dorder[3]{ 1 , 1 , 1 };
    // вычисляются значения функции и первой и второй производных
   // const double* datahint = NULL;
   // const MKL_INT rhint = DF_MATRIX_STORAGE_ROWS;
   // MKL_INT* cell = NULL;

    // Вычисление значений сплайна и двух производных в точках site
    int ndorder = 3;
    errcode = dfdInterpolate1D(task, DF_INTERP, DF_METHOD_PP,
        nsite, xUni, DF_SORTED_DATA,
        ndorder, dorder, NULL,
        res, DF_MATRIX_STORAGE_ROWS, NULL);

    /***** Delete Data Fitting task *****/
    errcode = dfDeleteTask(&task);

    printf("\ndfDeleteTask errcode = %d:\n", errcode); //CheckDfError(errcode);

                                                   /***** Print summary of the test *****/
    if (errcode != 0)
    {
        printf("\n\nError: Computed natural cubic spline coefficients");
        printf(" or integrals are incorrect\n");
        return_value = 1;
    }
    else
    {
        printf("\n\nComputed natural cubic spline coefficients");
        printf(" and integrals are correct\n");
        return_value = 0;
    }
    delete[] scoeff;
    /* delete[] llim;
     delete[] rlim;
     delete[] r;
     delete[] ref_r;*/
}

//////////////////////////////////
// Кубический сплайн только для одной функции (векторная функция с одним элементом)
// Только значения функции; Разные типы сплайнов; Разные условия на границах
// Кубический сплайн только для одной функции (векторная функция с одним элементом)
// Функция и две производные
extern "C"  _declspec(dllexport)

void Wrap_Splines_Lab2_2022V(int variant, int n, double* x, double* mData, double condLeft, double condRight,
    int nN, double* xUni, double* res, int& return_value)
{

    ofstream ofs("cpp_errors.txt", ios_base::out);
    //ofs.open(); // окрываем файл для записи
    ofs << "===== n = " << n << "  nN = " << nN << endl;
    ofs << "===== condLeft = " << condLeft << "  condRight = " << condRight << endl;

    for (int j = 0; j < n; j++)  ofs << "x[" << j << "] = " << x[j] << "  mData = " << mData[j] << endl;
    ofs << "===== nN = " << nN << endl;
    for (int j = 0; j < nN; j++)  ofs << "xUni[" << j << "] = " << xUni[j] << endl;
    /*  ofs << "===== ndorder = " << ndorder << endl;*/
    ofs.close();     // закрываем файл

    //double y[nx * ny];                   // значения функции (function values)

    //double left = xL;
    //double right = xR;
    // double freq;
    // double xx[N];
  /*  double left_val[N - 1], right_val[N - 1];
    double left_der1[N - 1], right_der1[N - 1];
    double left_der2[N - 1], right_der2[N - 1];*/

    //int i, j, errcode = 0;
    //int errnums = 0;

    /***** Initializing parameters for Data Fitting task *****/
    MKL_INT sorder = DF_PP_CUBIC;       // порядок сплайна (spline order)- 
    MKL_INT stype = DF_PP_NATURAL;      // тип сплайна (spline type)


    /***** Parameters describing interpolation interval *****/
    MKL_INT nx = n;                        // число узлов сплайна (number of break points)
    MKL_INT xhint = DF_NON_UNIFORM_PARTITION;   // additional info about break points
    //double xUniform[2];                         // границы интервала интерполирования (limits of the interpolation interval)
    //xUniform[0] = xL;
    //xUniform[1] = xR;
    /* Limits of interpolation interval are provided in case of uniform grid */
    //double xUniform[2];                        // границы интервала интерполирования (limits of the interpolation interval)
    //xUniform[0] = xL;
    //xUniform[1] = xR;

    /***** Parameters describing function *****/
    MKL_INT ny = 1;                     // размер вектор-функции (number of functions)
    MKL_INT yhint = DF_NO_HINT;          // additional info about function

    /***** Parameters describing spline coefficients storage *****/
    MKL_INT nscoeff = ny * (nx - 1) * DF_PP_CUBIC;   // число коэффициентов сплайна (number of spline coefficients)
    MKL_INT scoeffhint = DF_NO_HINT;                 // additional info about spline coefficients
    double* scoeff = new double[nscoeff];             // массив коэффициентов сплайна (array of spline coefficients)

    /***** Parameters describing boundary conditions type *****/
    MKL_INT bc_type = DF_BC_FREE_END;   // тип Free граничных условий (boundary conditions type)
    MKL_INT nbc = 2;// число граничных условий (number of boundary conditions)

    if (variant == 1)
    {
       bc_type = DF_BC_1ST_LEFT_DER | DF_BC_1ST_RIGHT_DER;       // тип граничных условий (boundary conditions type)
    }
    else if (variant == 2)
    {
       bc_type = DF_BC_2ND_LEFT_DER | DF_BC_2ND_RIGHT_DER;
    }
    else if (variant == 3)
    {
        bc_type = DF_BC_1ST_LEFT_DER | DF_BC_2ND_RIGHT_DER;
    }
    else // if (variant == 0)
    {
        bc_type = DF_BC_FREE_END;
    }
    /* No additional parameters are provided for Free-End boundary conditions */
    double* bc = new double[2]{ condLeft, condRight };                         // boundary conditions
  /*  bc = 0;*/


  /***** Parameters describing internal conditions type *****/
    /* No internal conditions are provided for natural cubic spline */
    MKL_INT ic_type = DF_NO_IC;          // тип условий во внутренних точках (internal conditions type)
    double* ic;                         // условия во внутренних точках (internal conditions)
    ic = 0;

    /***** Create Data Fitting task *****/
    DFTaskPtr task;                     // Data Fitting task descriptor
    int errcode = dfdNewTask1D(&task, nx, x, xhint, ny, mData, yhint);

    printf("\ndfdNewTask1D errcode = %d:\n", errcode); //CheckDfError(errcode);

    /***** Edit task parameters for natural cubic spline construction *****/
   // bc_type = DF_BC_1ST_RIGHT_DER;
   // bc = new double[1];
   //// bc[0] = 0;
   // bc[0] = der_1_at_right_end;

    errcode = dfdEditPPSpline1D(task, sorder, stype, bc_type, bc, ic_type, ic, scoeff, scoeffhint);

    // printf("\ndfdEditPPSpline1D errcode = %d:\n", errcode);  //CheckDfError(errcode);

     /***** Construct natural cubic spline using STD method *****/
    errcode = dfdConstruct1D(task, DF_PP_SPLINE, DF_METHOD_STD);
    // printf("\ndfdConstruct1D errcode = %d:\n", errcode); // CheckDfError(errcode);

    MKL_INT nsite = nN;
    // const MKL_INT sitehint = DF_SORTED_DATA;

     //MKL_INT *dorder;
     /*if (ndorder == 1) dorder = new MKL_INT[1]{ 1 };
     else if (ndorder == 2) dorder = new MKL_INT[2]{ 1 , 1};
     else if (ndorder == 3) dorder = new MKL_INT[3]{ 1 , 1 , 1};
     else if (ndorder == 4) dorder = new MKL_INT[4]{ 1 , 1 , 1 , 1};
     else*/
    MKL_INT dorder[3]{ 1 , 1 , 1 };
    // вычисляются значения функции и первой и второй производных
   // const double* datahint = NULL;
   // const MKL_INT rhint = DF_MATRIX_STORAGE_ROWS;
   // MKL_INT* cell = NULL;

    // Вычисление значений сплайна и двух производных в точках site
    int ndorder = 3;
    errcode = dfdInterpolate1D(task, DF_INTERP, DF_METHOD_PP,
        nsite, xUni, DF_SORTED_DATA,
        ndorder, dorder, NULL,
        res, DF_MATRIX_STORAGE_ROWS, NULL);

    /***** Delete Data Fitting task *****/
    errcode = dfDeleteTask(&task);

    printf("\ndfDeleteTask errcode = %d:\n", errcode); //CheckDfError(errcode);

                                                   /***** Print summary of the test *****/
    if (errcode != 0)
    {
        printf("\n\nError: Computed natural cubic spline coefficients");
        printf(" or integrals are incorrect\n");
        return_value = 1;
    }
    else
    {
        printf("\n\nComputed natural cubic spline coefficients");
        printf(" and integrals are correct\n");
        return_value = 0;
    }
    delete[] scoeff;
    /* delete[] llim;
     delete[] rlim;
     delete[] r;
     delete[] ref_r;*/
}

