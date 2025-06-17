/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * selfShadow_types.h
 *
 * Code generation for function 'selfShadow'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"

/* Type Definitions */
#ifndef typedef_struct0_T
#define typedef_struct0_T
typedef struct {
  real_T vertices[11430];
  real_T normal[22812];
  real_T faces[30416];
  real_T area[7604];
  real_T pos[22812];
  real_T uu[22812];
  real_T uv[22812];
  real_T qlb[30416];
  real_T F0[7604];
  real_T kappa[7604];
  real_T Ca[7604];
  real_T Cd[7604];
  real_T Cs[7604];
  real_T nu[7604];
  real_T nv[7604];
  real_T fObs[7604];
  real_T MOI[9];
  real_T m;
  real_T sunlitFlag[7604];
  real_T forces[22812];
  real_T torque[22812];
} struct0_T;
#endif /* typedef_struct0_T */

#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T
struct emxArray_real_T {
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};
#endif /* struct_emxArray_real_T */
#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T
typedef struct emxArray_real_T emxArray_real_T;
#endif /* typedef_emxArray_real_T */

#ifndef typedef_b_selfShadow_api
#define typedef_b_selfShadow_api
typedef struct {
  struct0_T sat;
} b_selfShadow_api;
#endif /* typedef_b_selfShadow_api */

#ifndef typedef_selfShadowStackData
#define typedef_selfShadowStackData
typedef struct {
  b_selfShadow_api f0;
} selfShadowStackData;
#endif /* typedef_selfShadowStackData */

/* End of code generation (selfShadow_types.h) */
