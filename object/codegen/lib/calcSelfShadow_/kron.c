/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: kron.c
 *
 * MATLAB Coder version            : 24.2
 * C/C++ source code generated on  : 2024/11/21 17:17:03
 */

/* Include Files */
#include "kron.h"

/* Function Definitions */
/*
 * Arguments    : const double B[3]
 *                double K[9]
 * Return Type  : void
 */
void kron(const double B[3], double K[9])
{
  int j2;
  int kidx;
  kidx = -1;
  for (j2 = 0; j2 < 3; j2++) {
    double d;
    d = B[j2];
    K[kidx + 1] = d;
    K[kidx + 2] = d;
    K[kidx + 3] = d;
    kidx += 3;
  }
}

/*
 * File trailer for kron.c
 *
 * [EOF]
 */
