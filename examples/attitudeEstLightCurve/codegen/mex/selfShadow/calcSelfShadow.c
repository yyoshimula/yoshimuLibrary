/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * calcSelfShadow.c
 *
 * Code generation for function 'calcSelfShadow'
 *
 */

/* Include files */
#include "calcSelfShadow.h"
#include "rt_nonfinite.h"

/* Function Definitions */
real_T calcSelfShadow(const real_T sun[3], const real_T nJ[3],
                      const real_T vertJ[9], const real_T vertI[3])
{
  real_T Q[3];
  real_T flag;
  real_T v2_idx_0;
  real_T v2_idx_3;
  int32_T j2;
  /*  scalar */
  /*  Calculate intersection scalar K */
  flag = (((nJ[0] * vertJ[0] + nJ[1] * vertJ[3]) + nJ[2] * vertJ[6]) -
          ((nJ[0] * vertI[0] + nJ[1] * vertI[1]) + nJ[2] * vertI[2])) /
         ((nJ[0] * sun[0] + nJ[1] * sun[1]) + nJ[2] * sun[2]);
  /*  scalar */
  /*  交点計算 */
  /*  3x1 */
  /* 	交点がメッシュ内にあるかどうかの計算 */
  v2_idx_0 = vertI[0] + flag * sun[0];
  Q[0] = v2_idx_0;
  v2_idx_3 = (v2_idx_0 - vertI[0]) * sun[0];
  v2_idx_0 = vertI[1] + flag * sun[1];
  Q[1] = v2_idx_0;
  v2_idx_3 += (v2_idx_0 - vertI[1]) * sun[1];
  v2_idx_0 = vertI[2] + flag * sun[2];
  Q[2] = v2_idx_0;
  v2_idx_3 += (v2_idx_0 - vertI[2]) * sun[2];
  if (v2_idx_3 > 0.0) {
    real_T v1[9];
    real_T crossV_idx_0;
    real_T crossV_idx_3;
    real_T v2_idx_1;
    real_T v2_idx_2;
    real_T v2_idx_4;
    real_T v2_idx_5;
    real_T v2_idx_7;
    real_T v2_idx_8;
    int32_T kidx;
    kidx = -1;
    for (j2 = 0; j2 < 3; j2++) {
      flag = Q[j2];
      v1[kidx + 1] = flag;
      v1[kidx + 2] = flag;
      v1[kidx + 3] = flag;
      kidx += 3;
    }
    for (j2 = 0; j2 < 9; j2++) {
      v1[j2] -= vertJ[j2];
    }
    /*  3x3 */
    /* v2 = 3x3 vector 3点を順番につなぐベクトル */
    v2_idx_0 = vertJ[1] - vertJ[0];
    v2_idx_1 = vertJ[2] - vertJ[1];
    v2_idx_2 = vertJ[0] - vertJ[2];
    v2_idx_3 = vertJ[4] - vertJ[3];
    v2_idx_4 = vertJ[5] - vertJ[4];
    v2_idx_5 = vertJ[3] - vertJ[5];
    flag = vertJ[7] - vertJ[6];
    v2_idx_7 = vertJ[8] - vertJ[7];
    v2_idx_8 = vertJ[6] - vertJ[8];
    /*  3x3 matrix */
    /*  外積v2×v1（= crossV）を求める // */
    crossV_idx_0 = v2_idx_3 * v1[6] - v1[3] * flag;
    crossV_idx_3 = v1[0] * flag - v2_idx_0 * v1[6];
    flag = v2_idx_0 * v1[3] - v1[0] * v2_idx_3;
    /*  3x3 */
    if (((crossV_idx_0 * (v2_idx_4 * v1[7] - v1[4] * v2_idx_7) +
          crossV_idx_3 * (v1[1] * v2_idx_7 - v2_idx_1 * v1[7])) +
             flag * (v2_idx_1 * v1[4] - v1[1] * v2_idx_4) >=
         0.0) &&
        ((crossV_idx_0 * (v2_idx_5 * v1[8] - v1[5] * v2_idx_8) +
          crossV_idx_3 * (v1[2] * v2_idx_8 - v2_idx_2 * v1[8])) +
             flag * (v2_idx_2 * v1[5] - v1[2] * v2_idx_5) >=
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

/* End of code generation (calcSelfShadow.c) */
