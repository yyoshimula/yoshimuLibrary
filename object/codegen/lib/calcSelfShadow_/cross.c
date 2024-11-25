/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: cross.c
 *
 * MATLAB Coder version            : 24.2
 * C/C++ source code generated on  : 2024/11/21 17:17:03
 */

/* Include Files */
#include "cross.h"

/* Function Definitions */
/*
 * Arguments    : const double a[9]
 *                const double b[9]
 *                double c[9]
 * Return Type  : void
 */
void cross(const double a[9], const double b[9], double c[9])
{
  c[0] = a[3] * b[6] - b[3] * a[6];
  c[3] = b[0] * a[6] - a[0] * b[6];
  c[6] = a[0] * b[3] - b[0] * a[3];
  c[1] = a[4] * b[7] - b[4] * a[7];
  c[4] = b[1] * a[7] - a[1] * b[7];
  c[7] = a[1] * b[4] - b[1] * a[4];
  c[2] = a[5] * b[8] - b[5] * a[8];
  c[5] = b[2] * a[8] - a[2] * b[8];
  c[8] = a[2] * b[5] - b[2] * a[5];
}

/*
 * File trailer for cross.c
 *
 * [EOF]
 */
