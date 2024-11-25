/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: calcSelfShadow__types.h
 *
 * MATLAB Coder version            : 24.2
 * C/C++ source code generated on  : 2024/11/21 17:17:03
 */

#ifndef CALCSELFSHADOW__TYPES_H
#define CALCSELFSHADOW__TYPES_H

/* Include Files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef struct_emxArray_real_T_1000x3
#define struct_emxArray_real_T_1000x3
struct emxArray_real_T_1000x3 {
  double data[3000];
  int size[2];
};
#endif /* struct_emxArray_real_T_1000x3 */
#ifndef typedef_emxArray_real_T_1000x3
#define typedef_emxArray_real_T_1000x3
typedef struct emxArray_real_T_1000x3 emxArray_real_T_1000x3;
#endif /* typedef_emxArray_real_T_1000x3 */

#ifndef struct_emxArray_real_T_1000
#define struct_emxArray_real_T_1000
struct emxArray_real_T_1000 {
  double data[1000];
  int size[1];
};
#endif /* struct_emxArray_real_T_1000 */
#ifndef typedef_emxArray_real_T_1000
#define typedef_emxArray_real_T_1000
typedef struct emxArray_real_T_1000 emxArray_real_T_1000;
#endif /* typedef_emxArray_real_T_1000 */

#ifndef typedef_struct0_T
#define typedef_struct0_T
typedef struct {
  emxArray_real_T_1000x3 normal;
  emxArray_real_T_1000x3 faces;
  emxArray_real_T_1000x3 pos;
  double vertices[9000];
  emxArray_real_T_1000 shadowFlag;
} struct0_T;
#endif /* typedef_struct0_T */

#endif
/*
 * File trailer for calcSelfShadow__types.h
 *
 * [EOF]
 */
