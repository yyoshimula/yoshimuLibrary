/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 24.2
 * C/C++ source code generated on  : 2024/11/21 17:17:03
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/

/* Include Files */
#include "main.h"
#include "calcSelfShadow_.h"
#include "calcSelfShadow__terminate.h"
#include "calcSelfShadow__types.h"
#include "selfShadow_.h"

/* Function Declarations */
static void argInit_3000x3_real_T(double result[9000]);

static void argInit_3x1_real_T(double result[3]);

static void argInit_3x3_real_T(double result[9]);

static int argInit_d1000x1_real_T(double result_data[]);

static void argInit_d1000x3_real_T(double result_data[], int result_size[2]);

static double argInit_real_T(void);

static void argInit_struct0_T(struct0_T *result);

/* Function Definitions */
/*
 * Arguments    : double result[9000]
 * Return Type  : void
 */
static void argInit_3000x3_real_T(double result[9000])
{
  int i;
  /* Loop over the array to initialize each element. */
  for (i = 0; i < 9000; i++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[i] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[3]
 * Return Type  : void
 */
static void argInit_3x1_real_T(double result[3])
{
  int idx0;
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[9]
 * Return Type  : void
 */
static void argInit_3x3_real_T(double result[9])
{
  int i;
  /* Loop over the array to initialize each element. */
  for (i = 0; i < 9; i++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[i] = argInit_real_T();
  }
}

/*
 * Arguments    : double result_data[]
 * Return Type  : int
 */
static int argInit_d1000x1_real_T(double result_data[])
{
  int idx0;
  int result_size;
  /* Set the size of the array.
Change this size to the value that the application requires. */
  result_size = 2;
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 2; idx0++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result_data[idx0] = argInit_real_T();
  }
  return result_size;
}

/*
 * Arguments    : double result_data[]
 *                int result_size[2]
 * Return Type  : void
 */
static void argInit_d1000x3_real_T(double result_data[], int result_size[2])
{
  int i;
  /* Set the size of the array.
Change this size to the value that the application requires. */
  result_size[0] = 2;
  result_size[1] = 3;
  /* Loop over the array to initialize each element. */
  for (i = 0; i < 6; i++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result_data[i] = argInit_real_T();
  }
}

/*
 * Arguments    : void
 * Return Type  : double
 */
static double argInit_real_T(void)
{
  return 0.0;
}

/*
 * Arguments    : struct0_T *result
 * Return Type  : void
 */
static void argInit_struct0_T(struct0_T *result)
{
  /* Set the value of each structure field.
Change this value to the value that the application requires. */
  argInit_d1000x3_real_T(result->normal.data, result->normal.size);
  argInit_d1000x3_real_T(result->faces.data, result->faces.size);
  argInit_d1000x3_real_T(result->pos.data, result->pos.size);
  argInit_3000x3_real_T(result->vertices);
  result->shadowFlag.size[0] = argInit_d1000x1_real_T(result->shadowFlag.data);
}

/*
 * Arguments    : int argc
 *                char **argv
 * Return Type  : int
 */
int main(int argc, char **argv)
{
  (void)argc;
  (void)argv;
  /* The initialize function is being called automatically from your entry-point
   * function. So, a call to initialize is not included here. */
  /* Invoke the entry-point functions.
You can call entry-point functions multiple times. */
  main_calcSelfShadow_();
  main_selfShadow_();
  /* Terminate the application.
You do not need to do this more than one time. */
  calcSelfShadow__terminate();
  return 0;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void main_calcSelfShadow_(void)
{
  double dv1[9];
  double dv[3];
  double nJ_tmp[3];
  double flag;
  /* Initialize function 'calcSelfShadow_' input arguments. */
  /* Initialize function input argument 'sun'. */
  /* Initialize function input argument 'nJ'. */
  argInit_3x1_real_T(nJ_tmp);
  /* Initialize function input argument 'vertJ'. */
  /* Initialize function input argument 'vertI'. */
  /* Call the entry-point 'calcSelfShadow_'. */
  argInit_3x1_real_T(dv);
  argInit_3x3_real_T(dv1);
  flag = calcSelfShadow_(dv, nJ_tmp, dv1, nJ_tmp);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void main_selfShadow_(void)
{
  struct0_T sat;
  double dv[3];
  /* Initialize function 'selfShadow_' input arguments. */
  /* Initialize function input argument 'sat'. */
  argInit_struct0_T(&sat);
  /* Initialize function input argument 'sun'. */
  /* Call the entry-point 'selfShadow_'. */
  argInit_3x1_real_T(dv);
  selfShadow_(&sat, dv);
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
