/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * calcSelfShadow.h
 *
 * Code generation for function 'calcSelfShadow'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
real_T calcSelfShadow(const real_T sun[3], const real_T nJ[3],
                      const real_T vertJ[9], const real_T vertI[3]);

/* End of code generation (calcSelfShadow.h) */
