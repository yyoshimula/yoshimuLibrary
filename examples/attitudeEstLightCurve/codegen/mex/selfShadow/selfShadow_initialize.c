/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * selfShadow_initialize.c
 *
 * Code generation for function 'selfShadow_initialize'
 *
 */

/* Include files */
#include "selfShadow_initialize.h"
#include "_coder_selfShadow_mex.h"
#include "rt_nonfinite.h"
#include "selfShadow_data.h"

/* Function Declarations */
static void selfShadow_once(void);

/* Function Definitions */
static void selfShadow_once(void)
{
  mex_InitInfAndNan();
}

void selfShadow_initialize(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2022b(&st);
  emlrtClearAllocCountR2012b(&st, false, 0U, NULL);
  emlrtEnterRtStackR2012b(&st);
  if (emlrtFirstTimeR2012b(emlrtRootTLSGlobal)) {
    selfShadow_once();
  }
}

/* End of code generation (selfShadow_initialize.c) */
