/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * selfShadow.h
 *
 * Code generation for function 'selfShadow'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "selfShadow_types.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
emlrtCTX emlrtGetRootTLSGlobal(void);

void emlrtLockerFunction(EmlrtLockeeFunction aLockee, emlrtConstCTX aTLS,
                         void *aData);

void selfShadow(const emlrtStack *sp, struct0_T *sat, const real_T sun[3]);

/* End of code generation (selfShadow.h) */
