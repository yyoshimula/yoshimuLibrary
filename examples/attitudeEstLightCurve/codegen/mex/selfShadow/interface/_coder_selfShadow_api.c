/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_selfShadow_api.c
 *
 * Code generation for function '_coder_selfShadow_api'
 *
 */

/* Include files */
#include "_coder_selfShadow_api.h"
#include "rt_nonfinite.h"
#include "selfShadow.h"
#include "selfShadow_data.h"
#include "selfShadow_types.h"

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               struct0_T *y);

static const mxArray *b_emlrt_marshallOut(const real_T u[22812]);

static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               real_T y[11430]);

static const mxArray *c_emlrt_marshallOut(const real_T u[30416]);

static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               real_T y[22812]);

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               real_T y[30416]);

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                             const char_T *identifier, struct0_T *y);

static const mxArray *emlrt_marshallOut(const struct0_T *u);

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               real_T y[7604]);

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[9]);

static real_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[3];

static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[3];

static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               real_T ret[11430]);

static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               real_T ret[22812]);

static void m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               real_T ret[30416]);

static void n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               real_T ret[7604]);

static void o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId, real_T ret[9]);

static real_T p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

static real_T (*q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[3];

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, struct0_T *y)
{
  static const int32_T dims = 0;
  static const char_T *fieldNames[21] = {
      "vertices", "normal", "faces", "area", "pos",        "uu",     "uv",
      "qlb",      "F0",     "kappa", "Ca",   "Cd",         "Cs",     "nu",
      "nv",       "fObs",   "MOI",   "m",    "sunlitFlag", "forces", "torque"};
  emlrtMsgIdentifier thisId;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b((emlrtConstCTX)sp, parentId, u, 21,
                         (const char_T **)&fieldNames[0], 0U,
                         (const void *)&dims);
  thisId.fIdentifier = "vertices";
  c_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 0, "vertices")),
      &thisId, y->vertices);
  thisId.fIdentifier = "normal";
  d_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 1, "normal")),
      &thisId, y->normal);
  thisId.fIdentifier = "faces";
  e_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 2, "faces")),
      &thisId, y->faces);
  thisId.fIdentifier = "area";
  f_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 3, "area")),
      &thisId, y->area);
  thisId.fIdentifier = "pos";
  d_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 4, "pos")),
      &thisId, y->pos);
  thisId.fIdentifier = "uu";
  d_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 5, "uu")),
      &thisId, y->uu);
  thisId.fIdentifier = "uv";
  d_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 6, "uv")),
      &thisId, y->uv);
  thisId.fIdentifier = "qlb";
  e_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 7, "qlb")),
      &thisId, y->qlb);
  thisId.fIdentifier = "F0";
  f_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 8, "F0")),
      &thisId, y->F0);
  thisId.fIdentifier = "kappa";
  f_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 9, "kappa")),
      &thisId, y->kappa);
  thisId.fIdentifier = "Ca";
  f_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 10, "Ca")),
      &thisId, y->Ca);
  thisId.fIdentifier = "Cd";
  f_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 11, "Cd")),
      &thisId, y->Cd);
  thisId.fIdentifier = "Cs";
  f_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 12, "Cs")),
      &thisId, y->Cs);
  thisId.fIdentifier = "nu";
  f_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 13, "nu")),
      &thisId, y->nu);
  thisId.fIdentifier = "nv";
  f_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 14, "nv")),
      &thisId, y->nv);
  thisId.fIdentifier = "fObs";
  f_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 15, "fObs")),
      &thisId, y->fObs);
  thisId.fIdentifier = "MOI";
  g_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 16, "MOI")),
      &thisId, y->MOI);
  thisId.fIdentifier = "m";
  y->m = h_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 17, "m")),
      &thisId);
  thisId.fIdentifier = "sunlitFlag";
  f_emlrt_marshallIn(sp,
                     emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 18,
                                                    "sunlitFlag")),
                     &thisId, y->sunlitFlag);
  thisId.fIdentifier = "forces";
  d_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 19, "forces")),
      &thisId, y->forces);
  thisId.fIdentifier = "torque";
  d_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 20, "torque")),
      &thisId, y->torque);
  emlrtDestroyArray(&u);
}

static const mxArray *b_emlrt_marshallOut(const real_T u[22812])
{
  static const int32_T iv[2] = {7604, 3};
  const mxArray *m;
  const mxArray *y;
  real_T *pData;
  int32_T b_i;
  int32_T c_i;
  int32_T i;
  y = NULL;
  m = emlrtCreateNumericArray(2, (const void *)&iv[0], mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  i = 0;
  for (b_i = 0; b_i < 3; b_i++) {
    for (c_i = 0; c_i < 7604; c_i++) {
      pData[i + c_i] = u[c_i + 7604 * b_i];
    }
    i += 7604;
  }
  emlrtAssign(&y, m);
  return y;
}

static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               real_T y[11430])
{
  k_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *c_emlrt_marshallOut(const real_T u[30416])
{
  static const int32_T iv[2] = {7604, 4};
  const mxArray *m;
  const mxArray *y;
  real_T *pData;
  int32_T b_i;
  int32_T c_i;
  int32_T i;
  y = NULL;
  m = emlrtCreateNumericArray(2, (const void *)&iv[0], mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  i = 0;
  for (b_i = 0; b_i < 4; b_i++) {
    for (c_i = 0; c_i < 7604; c_i++) {
      pData[i + c_i] = u[c_i + 7604 * b_i];
    }
    i += 7604;
  }
  emlrtAssign(&y, m);
  return y;
}

static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               real_T y[22812])
{
  l_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               real_T y[30416])
{
  m_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                             const char_T *identifier, struct0_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId, y);
  emlrtDestroyArray(&nullptr);
}

static const mxArray *emlrt_marshallOut(const struct0_T *u)
{
  static const int32_T iv[2] = {3810, 3};
  static const int32_T iv1[2] = {3, 3};
  static const int32_T i1 = 7604;
  static const int32_T i10 = 7604;
  static const int32_T i2 = 7604;
  static const int32_T i3 = 7604;
  static const int32_T i4 = 7604;
  static const int32_T i5 = 7604;
  static const int32_T i6 = 7604;
  static const int32_T i7 = 7604;
  static const int32_T i8 = 7604;
  static const int32_T i9 = 7604;
  static const char_T *sv[21] = {
      "vertices", "normal", "faces", "area", "pos",        "uu",     "uv",
      "qlb",      "F0",     "kappa", "Ca",   "Cd",         "Cs",     "nu",
      "nv",       "fObs",   "MOI",   "m",    "sunlitFlag", "forces", "torque"};
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *e_y;
  const mxArray *f_y;
  const mxArray *g_y;
  const mxArray *h_y;
  const mxArray *i_y;
  const mxArray *j_y;
  const mxArray *k_y;
  const mxArray *l_y;
  const mxArray *m;
  const mxArray *m_y;
  const mxArray *n_y;
  const mxArray *y;
  real_T *pData;
  int32_T b_i;
  int32_T c_i;
  int32_T i;
  y = NULL;
  emlrtAssign(&y, emlrtCreateStructMatrix(1, 1, 21, (const char_T **)&sv[0]));
  b_y = NULL;
  m = emlrtCreateNumericArray(2, (const void *)&iv[0], mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  i = 0;
  for (b_i = 0; b_i < 3; b_i++) {
    for (c_i = 0; c_i < 3810; c_i++) {
      pData[i + c_i] = u->vertices[c_i + 3810 * b_i];
    }
    i += 3810;
  }
  emlrtAssign(&b_y, m);
  emlrtSetFieldR2017b(y, 0, "vertices", b_y, 0);
  emlrtSetFieldR2017b(y, 0, "normal", b_emlrt_marshallOut(u->normal), 1);
  emlrtSetFieldR2017b(y, 0, "faces", c_emlrt_marshallOut(u->faces), 2);
  c_y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i1, mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  for (b_i = 0; b_i < 7604; b_i++) {
    pData[b_i] = u->area[b_i];
  }
  emlrtAssign(&c_y, m);
  emlrtSetFieldR2017b(y, 0, "area", c_y, 3);
  emlrtSetFieldR2017b(y, 0, "pos", b_emlrt_marshallOut(u->pos), 4);
  emlrtSetFieldR2017b(y, 0, "uu", b_emlrt_marshallOut(u->uu), 5);
  emlrtSetFieldR2017b(y, 0, "uv", b_emlrt_marshallOut(u->uv), 6);
  emlrtSetFieldR2017b(y, 0, "qlb", c_emlrt_marshallOut(u->qlb), 7);
  d_y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i2, mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  for (b_i = 0; b_i < 7604; b_i++) {
    pData[b_i] = u->F0[b_i];
  }
  emlrtAssign(&d_y, m);
  emlrtSetFieldR2017b(y, 0, "F0", d_y, 8);
  e_y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i3, mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  for (b_i = 0; b_i < 7604; b_i++) {
    pData[b_i] = u->kappa[b_i];
  }
  emlrtAssign(&e_y, m);
  emlrtSetFieldR2017b(y, 0, "kappa", e_y, 9);
  f_y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i4, mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  for (b_i = 0; b_i < 7604; b_i++) {
    pData[b_i] = u->Ca[b_i];
  }
  emlrtAssign(&f_y, m);
  emlrtSetFieldR2017b(y, 0, "Ca", f_y, 10);
  g_y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i5, mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  for (b_i = 0; b_i < 7604; b_i++) {
    pData[b_i] = u->Cd[b_i];
  }
  emlrtAssign(&g_y, m);
  emlrtSetFieldR2017b(y, 0, "Cd", g_y, 11);
  h_y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i6, mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  for (b_i = 0; b_i < 7604; b_i++) {
    pData[b_i] = u->Cs[b_i];
  }
  emlrtAssign(&h_y, m);
  emlrtSetFieldR2017b(y, 0, "Cs", h_y, 12);
  i_y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i7, mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  for (b_i = 0; b_i < 7604; b_i++) {
    pData[b_i] = u->nu[b_i];
  }
  emlrtAssign(&i_y, m);
  emlrtSetFieldR2017b(y, 0, "nu", i_y, 13);
  j_y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i8, mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  for (b_i = 0; b_i < 7604; b_i++) {
    pData[b_i] = u->nv[b_i];
  }
  emlrtAssign(&j_y, m);
  emlrtSetFieldR2017b(y, 0, "nv", j_y, 14);
  k_y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i9, mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  for (b_i = 0; b_i < 7604; b_i++) {
    pData[b_i] = u->fObs[b_i];
  }
  emlrtAssign(&k_y, m);
  emlrtSetFieldR2017b(y, 0, "fObs", k_y, 15);
  l_y = NULL;
  m = emlrtCreateNumericArray(2, (const void *)&iv1[0], mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  i = 0;
  for (b_i = 0; b_i < 3; b_i++) {
    pData[i] = u->MOI[3 * b_i];
    pData[i + 1] = u->MOI[3 * b_i + 1];
    pData[i + 2] = u->MOI[3 * b_i + 2];
    i += 3;
  }
  emlrtAssign(&l_y, m);
  emlrtSetFieldR2017b(y, 0, "MOI", l_y, 16);
  m_y = NULL;
  m = emlrtCreateDoubleScalar(u->m);
  emlrtAssign(&m_y, m);
  emlrtSetFieldR2017b(y, 0, "m", m_y, 17);
  n_y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i10, mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  for (b_i = 0; b_i < 7604; b_i++) {
    pData[b_i] = u->sunlitFlag[b_i];
  }
  emlrtAssign(&n_y, m);
  emlrtSetFieldR2017b(y, 0, "sunlitFlag", n_y, 18);
  emlrtSetFieldR2017b(y, 0, "forces", b_emlrt_marshallOut(u->forces), 19);
  emlrtSetFieldR2017b(y, 0, "torque", b_emlrt_marshallOut(u->torque), 20);
  return y;
}

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               real_T y[7604])
{
  n_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[9])
{
  o_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = p_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[3]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[3];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = j_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[3]
{
  real_T(*y)[3];
  y = q_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               real_T ret[11430])
{
  static const int32_T dims[2] = {3810, 3};
  real_T(*r)[11430];
  int32_T i;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                          (const void *)&dims[0]);
  r = (real_T(*)[11430])emlrtMxGetData(src);
  for (i = 0; i < 11430; i++) {
    ret[i] = (*r)[i];
  }
  emlrtDestroyArray(&src);
}

static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               real_T ret[22812])
{
  static const int32_T dims[2] = {7604, 3};
  real_T(*r)[22812];
  int32_T i;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                          (const void *)&dims[0]);
  r = (real_T(*)[22812])emlrtMxGetData(src);
  for (i = 0; i < 22812; i++) {
    ret[i] = (*r)[i];
  }
  emlrtDestroyArray(&src);
}

static void m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               real_T ret[30416])
{
  static const int32_T dims[2] = {7604, 4};
  real_T(*r)[30416];
  int32_T i;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                          (const void *)&dims[0]);
  r = (real_T(*)[30416])emlrtMxGetData(src);
  for (i = 0; i < 30416; i++) {
    ret[i] = (*r)[i];
  }
  emlrtDestroyArray(&src);
}

static void n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               real_T ret[7604])
{
  static const int32_T dims = 7604;
  real_T(*r)[7604];
  int32_T i;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 1U,
                          (const void *)&dims);
  r = (real_T(*)[7604])emlrtMxGetData(src);
  for (i = 0; i < 7604; i++) {
    ret[i] = (*r)[i];
  }
  emlrtDestroyArray(&src);
}

static void o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId, real_T ret[9])
{
  static const int32_T dims[2] = {3, 3};
  real_T(*r)[9];
  int32_T i;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                          (const void *)&dims[0]);
  r = (real_T(*)[9])emlrtMxGetData(src);
  for (i = 0; i < 9; i++) {
    ret[i] = (*r)[i];
  }
  emlrtDestroyArray(&src);
}

static real_T p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  real_T ret;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 0U,
                          (const void *)&dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[3]
{
  static const int32_T dims[2] = {1, 3};
  real_T(*ret)[3];
  int32_T iv[2];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[3])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void selfShadow_api(selfShadowStackData *SD, const mxArray *const prhs[2],
                    const mxArray **plhs)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  real_T(*sun)[3];
  st.tls = emlrtRootTLSGlobal;
  /* Marshall function inputs */
  emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "sat", &SD->f0.sat);
  sun = i_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "sun");
  /* Invoke the target function */
  selfShadow(&st, &SD->f0.sat, *sun);
  /* Marshall function outputs */
  *plhs = emlrt_marshallOut(&SD->f0.sat);
}

/* End of code generation (_coder_selfShadow_api.c) */
