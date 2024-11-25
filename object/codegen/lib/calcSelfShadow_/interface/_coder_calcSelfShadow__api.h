/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_calcSelfShadow__api.h
 *
 * MATLAB Coder version            : 24.2
 * C/C++ source code generated on  : 2024/11/21 17:17:03
 */

#ifndef _CODER_CALCSELFSHADOW__API_H
#define _CODER_CALCSELFSHADOW__API_H

/* Include Files */
#include "emlrt.h"
#include "mex.h"
#include "tmwtypes.h"
#include <string.h>

/* Type Definitions */
#ifndef struct_emxArray_real_T_1000x3
#define struct_emxArray_real_T_1000x3
struct emxArray_real_T_1000x3 {
  real_T data[3000];
  int32_T size[2];
};
#endif /* struct_emxArray_real_T_1000x3 */
#ifndef typedef_emxArray_real_T_1000x3
#define typedef_emxArray_real_T_1000x3
typedef struct emxArray_real_T_1000x3 emxArray_real_T_1000x3;
#endif /* typedef_emxArray_real_T_1000x3 */

#ifndef struct_emxArray_real_T_1000
#define struct_emxArray_real_T_1000
struct emxArray_real_T_1000 {
  real_T data[1000];
  int32_T size[1];
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
  real_T vertices[9000];
  emxArray_real_T_1000 shadowFlag;
} struct0_T;
#endif /* typedef_struct0_T */

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
real_T calcSelfShadow_(real_T sun[3], real_T nJ[3], real_T vertJ[9],
                       real_T vertI[3]);

void calcSelfShadow__api(const mxArray *const prhs[4], const mxArray **plhs);

void calcSelfShadow__atexit(void);

void calcSelfShadow__initialize(void);

void calcSelfShadow__terminate(void);

void calcSelfShadow__xil_shutdown(void);

void calcSelfShadow__xil_terminate(void);

void selfShadow_(struct0_T *sat, real_T sun[3]);

void selfShadow__api(const mxArray *const prhs[2], const mxArray **plhs);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for _coder_calcSelfShadow__api.h
 *
 * [EOF]
 */
