/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: calcSelfShadow_.c
 *
 * MATLAB Coder version            : 24.2
 * C/C++ source code generated on  : 2024/11/21 17:17:03
 */

/* Include Files */
#include "calcSelfShadow_.h"
#include "cross.h"
#include "kron.h"

/* Function Definitions */
/*
 * ----------------------------------------------------------------------
 *
 *     20200815  y.yoshimura
 *     Inputs:
 *    Outputs:
 *    related function files:
 *    note:
 *    cf:
 *    revisions;
 *
 *    (c) 2020 yasuhiro yoshimura
 * ----------------------------------------------------------------------
 *  calculate shadow if j-th face make the shadow on the i-th face
 *  j番目のfacetがi番目のfacetに影を作るかを判定
 *
 * Arguments    : const double sun[3]
 *                const double nJ[3]
 *                const double vertJ[9]
 *                const double vertI[3]
 * Return Type  : double
 */
double calcSelfShadow_(const double sun[3], const double nJ[3],
                       const double vertJ[9], const double vertI[3])
{
  double Q[3];
  double b_Q;
  double d;
  double flag;
  int i;
  /*  3x1 vector */
  /*  Calculate intersection scalar K */
  flag = (((nJ[0] * vertJ[0] + nJ[1] * vertJ[3]) + nJ[2] * vertJ[6]) -
          ((nJ[0] * vertI[0] + nJ[1] * vertI[1]) + nJ[2] * vertI[2])) /
         ((nJ[0] * sun[0] + nJ[1] * sun[1]) + nJ[2] * sun[2]);
  /*  scalar */
  /*  交点計算 */
  /*  3x1 */
  /* 	交点がメッシュ内にあるかどうかの計算 */
  d = vertI[0] + flag * sun[0];
  Q[0] = d;
  b_Q = (d - vertI[0]) * sun[0];
  d = vertI[1] + flag * sun[1];
  Q[1] = d;
  b_Q += (d - vertI[1]) * sun[1];
  d = vertI[2] + flag * sun[2];
  Q[2] = d;
  b_Q += (d - vertI[2]) * sun[2];
  if (b_Q > 0.0) {
    double b_crossV[9];
    double b_vertJ[9];
    double crossV[9];
    /*  3x3 */
    /* v2 = 3x3 vector 3点を順番につなぐベクトル */
    /*  3x3 matrix */
    /*  外積v2×v1（= crossV）を求める // */
    kron(Q, crossV);
    b_vertJ[0] = vertJ[1] - vertJ[0];
    b_vertJ[1] = vertJ[2] - vertJ[1];
    b_vertJ[2] = vertJ[0] - vertJ[2];
    b_vertJ[3] = vertJ[4] - vertJ[3];
    b_vertJ[4] = vertJ[5] - vertJ[4];
    b_vertJ[5] = vertJ[3] - vertJ[5];
    b_vertJ[6] = vertJ[7] - vertJ[6];
    b_vertJ[7] = vertJ[8] - vertJ[7];
    b_vertJ[8] = vertJ[6] - vertJ[8];
    for (i = 0; i < 9; i++) {
      b_crossV[i] = crossV[i] - vertJ[i];
    }
    cross(b_vertJ, b_crossV, crossV);
    /*  3x3 */
    if (((crossV[0] * crossV[1] + crossV[3] * crossV[4]) +
             crossV[6] * crossV[7] >=
         0.0) &&
        ((crossV[0] * crossV[2] + crossV[3] * crossV[5]) +
             crossV[6] * crossV[8] >=
         0.0)) {
      flag = 0.0;
    } else {
      flag = 1.0;
    }
  } else {
    flag = 1.0;
  }
  return flag;
}

/*
 * File trailer for calcSelfShadow_.c
 *
 * [EOF]
 */
