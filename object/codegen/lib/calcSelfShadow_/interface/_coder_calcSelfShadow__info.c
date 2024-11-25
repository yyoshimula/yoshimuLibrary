/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_calcSelfShadow__info.c
 *
 * MATLAB Coder version            : 24.2
 * C/C++ source code generated on  : 2024/11/21 17:17:03
 */

/* Include Files */
#include "_coder_calcSelfShadow__info.h"
#include "emlrt.h"
#include "tmwtypes.h"

/* Function Declarations */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void);

/* Function Definitions */
/*
 * Arguments    : void
 * Return Type  : const mxArray *
 */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void)
{
  const mxArray *nameCaptureInfo;
  const char_T *data[4] = {
      "789cc553cb4e8340149d1a6cdc54bbf2339c9818f7ad56d3a4f505da85310accada0330c"
      "0cb4c21fb8f31f5cfb01fa45fe86521e0592094d5aebdd1c4e0e33e7"
      "de730135fac30642681b25f5da4cb095f2768a1ba85c55bd21c1ac3691523a97e96f299a"
      "dc09200c12e2e80cf2938433dbd19d408b5c40027c4ea74066cad8a6",
      "a0d90cd422398b193b29483989a5f8f9c802f3599d30242c7fde212d923c8f4fc9bcca82"
      "79a8923cda15fdb67787af7d103e8e22ee5b369b501d1f0bee1a3cc4"
      "c38e36e8745320c29e024e5f1ad886d04584b9f10466807da063d5d2097fb9df63a539dc"
      "25e7d8a99923d34d9d9aeabc8bdcff6149ffa6d43f51089f181456b7",
      "b791d4afacaf666f95d47e7757975b9677dd1cb2ffb185b666f8fda53cc6b82e3ff3fda3"
      "b94ebfacfecb2f94dcb7e877b82bf16b57f4abfd51771879c43bec7b"
      "371e39ef5df60707a7f33e2e6a7ceafa4012fed7f7ff00f12a822b",
      ""};
  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(&data[0], 1696U, &nameCaptureInfo);
  return nameCaptureInfo;
}

/*
 * Arguments    : void
 * Return Type  : mxArray *
 */
mxArray *emlrtMexFcnProperties(void)
{
  mxArray *xEntryPoints;
  mxArray *xInputs;
  mxArray *xResult;
  const char_T *propFieldName[9] = {"Version",
                                    "ResolvedFunctions",
                                    "Checksum",
                                    "EntryPoints",
                                    "CoverageInfo",
                                    "IsPolymorphic",
                                    "PropertyList",
                                    "UUID",
                                    "ClassEntryPointIsHandle"};
  const char_T *epFieldName[8] = {
      "QualifiedName",    "NumberOfInputs", "NumberOfOutputs", "ConstantInputs",
      "ResolvedFilePath", "TimeStamp",      "Constructor",     "Visible"};
  xEntryPoints =
      emlrtCreateStructMatrix(1, 2, 8, (const char_T **)&epFieldName[0]);
  xInputs = emlrtCreateLogicalMatrix(1, 4);
  emlrtSetField(xEntryPoints, 0, "QualifiedName",
                emlrtMxCreateString("calcSelfShadow_"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(4.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(1.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 0, "ResolvedFilePath",
      emlrtMxCreateString("/Users/yyoshimula/Dropbox/MATLAB/MATLABdrive/"
                          "yoshimuLibrary/object/calcSelfShadow_.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp",
                emlrtMxCreateDoubleScalar(739533.57619212964));
  emlrtSetField(xEntryPoints, 0, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 0, "Visible", emlrtMxCreateLogicalScalar(true));
  xInputs = emlrtCreateLogicalMatrix(1, 2);
  emlrtSetField(xEntryPoints, 1, "QualifiedName",
                emlrtMxCreateString("selfShadow_"));
  emlrtSetField(xEntryPoints, 1, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(2.0));
  emlrtSetField(xEntryPoints, 1, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(1.0));
  emlrtSetField(xEntryPoints, 1, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 1, "ResolvedFilePath",
      emlrtMxCreateString("/Users/yyoshimula/Dropbox/MATLAB/MATLABdrive/"
                          "yoshimuLibrary/object/selfShadow_.m"));
  emlrtSetField(xEntryPoints, 1, "TimeStamp",
                emlrtMxCreateDoubleScalar(739577.71975694445));
  emlrtSetField(xEntryPoints, 1, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 1, "Visible", emlrtMxCreateLogicalScalar(true));
  xResult =
      emlrtCreateStructMatrix(1, 1, 9, (const char_T **)&propFieldName[0]);
  emlrtSetField(xResult, 0, "Version",
                emlrtMxCreateString("24.2.0.2712019 (R2024b)"));
  emlrtSetField(xResult, 0, "ResolvedFunctions",
                (mxArray *)c_emlrtMexFcnResolvedFunctionsI());
  emlrtSetField(xResult, 0, "Checksum",
                emlrtMxCreateString("Jl4vrdr989lzswlkNXh7OB"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/*
 * File trailer for _coder_calcSelfShadow__info.c
 *
 * [EOF]
 */
