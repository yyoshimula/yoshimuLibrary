/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mtimes.c
 *
 * MATLAB Coder version            : 24.2
 * C/C++ source code generated on  : 2024/11/21 17:17:03
 */

/* Include Files */
#include "mtimes.h"

/* Function Definitions */
/*
 * Arguments    : const double A_data[]
 *                const int A_size[2]
 *                const double B[3]
 *                double C_data[]
 * Return Type  : int
 */
int mtimes(const double A_data[], const int A_size[2], const double B[3],
           double C_data[])
{
  int C_size;
  int i;
  C_size = A_size[0];
  for (i = 0; i < C_size; i++) {
    C_data[i] = (A_data[i] * B[0] + A_data[A_size[0] + i] * B[1]) +
                A_data[(A_size[0] << 1) + i] * B[2];
  }
  return C_size;
}

/*
 * File trailer for mtimes.c
 *
 * [EOF]
 */
