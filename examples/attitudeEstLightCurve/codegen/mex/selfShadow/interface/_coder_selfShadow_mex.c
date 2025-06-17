/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_selfShadow_mex.c
 *
 * Code generation for function '_coder_selfShadow_mex'
 *
 */

/* Include files */
#include "_coder_selfShadow_mex.h"
#include "_coder_selfShadow_api.h"
#include "rt_nonfinite.h"
#include "selfShadow.h"
#include "selfShadow_data.h"
#include "selfShadow_initialize.h"
#include "selfShadow_terminate.h"
#include "selfShadow_types.h"
#include "omp.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  static jmp_buf emlrtJBEnviron;
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  selfShadowStackData *selfShadowStackDataGlobal = NULL;
  selfShadowStackDataGlobal = (selfShadowStackData *)emlrtMxCalloc(
      (size_t)1, (size_t)1U * sizeof(selfShadowStackData));
  mexAtExit(&selfShadow_atexit);
  /* Initialize the memory manager. */
  omp_init_lock(&emlrtLockGlobal);
  omp_init_nest_lock(&selfShadow_nestLockGlobal);
  selfShadow_initialize();
  st.tls = emlrtRootTLSGlobal;
  emlrtSetJmpBuf(&st, &emlrtJBEnviron);
  if (setjmp(emlrtJBEnviron) == 0) {
    selfShadow_mexFunction(selfShadowStackDataGlobal, nlhs, plhs, nrhs, prhs);
    selfShadow_terminate();
    omp_destroy_lock(&emlrtLockGlobal);
    omp_destroy_nest_lock(&selfShadow_nestLockGlobal);
  } else {
    omp_destroy_lock(&emlrtLockGlobal);
    omp_destroy_nest_lock(&selfShadow_nestLockGlobal);
    emlrtReportParallelRunTimeError(&st);
  }
  emlrtMxFree(selfShadowStackDataGlobal);
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal,
                           &emlrtLockerFunction, omp_get_num_procs(), NULL,
                           "UTF-8", true);
  return emlrtRootTLSGlobal;
}

void selfShadow_mexFunction(selfShadowStackData *SD, int32_T nlhs,
                            mxArray *plhs[1], int32_T nrhs,
                            const mxArray *prhs[2])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 2, 4,
                        10, "selfShadow");
  }
  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 10,
                        "selfShadow");
  }
  /* Call the function. */
  selfShadow_api(SD, prhs, &outputs);
  /* Copy over outputs to the caller. */
  emlrtReturnArrays(1, &plhs[0], &outputs);
}

/* End of code generation (_coder_selfShadow_mex.c) */
