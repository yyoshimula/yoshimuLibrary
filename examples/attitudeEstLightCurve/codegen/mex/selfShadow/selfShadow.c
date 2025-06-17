/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * selfShadow.c
 *
 * Code generation for function 'selfShadow'
 *
 */

/* Include files */
#include "selfShadow.h"
#include "calcSelfShadow.h"
#include "rt_nonfinite.h"
#include "selfShadow_data.h"
#include "selfShadow_emxutil.h"
#include "selfShadow_types.h"
#include "blas.h"
#include "mwmathutil.h"
#include "omp.h"
#include <stddef.h>
#include <string.h>

/* Variable Definitions */
static emlrtBCInfo
    emlrtBCI =
        {
            -1,           /* iFirst */
            -1,           /* iLast */
            18,           /* lineNo */
            16,           /* colNo */
            "pos",        /* aName */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            0                 /* checkKind */
};

static emlrtDCInfo
    emlrtDCI =
        {
            21,           /* lineNo */
            35,           /* colNo */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            1                 /* checkKind */
};

static emlrtBCInfo
    b_emlrtBCI =
        {
            1,              /* iFirst */
            3810,           /* iLast */
            21,             /* lineNo */
            35,             /* colNo */
            "sat.vertices", /* aName */
            "selfShadow",   /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            0                 /* checkKind */
};

static emlrtDCInfo
    b_emlrtDCI =
        {
            22,           /* lineNo */
            30,           /* colNo */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            1                 /* checkKind */
};

static emlrtBCInfo
    c_emlrtBCI =
        {
            1,              /* iFirst */
            3810,           /* iLast */
            22,             /* lineNo */
            30,             /* colNo */
            "sat.vertices", /* aName */
            "selfShadow",   /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            0                 /* checkKind */
};

static emlrtDCInfo
    c_emlrtDCI =
        {
            23,           /* lineNo */
            30,           /* colNo */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            1                 /* checkKind */
};

static emlrtBCInfo
    d_emlrtBCI =
        {
            1,              /* iFirst */
            3810,           /* iLast */
            23,             /* lineNo */
            30,             /* colNo */
            "sat.vertices", /* aName */
            "selfShadow",   /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            0                 /* checkKind */
};

static emlrtBCInfo
    e_emlrtBCI =
        {
            -1,           /* iFirst */
            -1,           /* iLast */
            24,           /* lineNo */
            48,           /* colNo */
            "normal",     /* aName */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            0                 /* checkKind */
};

static emlrtDCInfo
    d_emlrtDCI =
        {
            29,           /* lineNo */
            39,           /* colNo */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            1                 /* checkKind */
};

static emlrtBCInfo
    f_emlrtBCI =
        {
            1,              /* iFirst */
            3810,           /* iLast */
            29,             /* lineNo */
            39,             /* colNo */
            "sat.vertices", /* aName */
            "selfShadow",   /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            0                 /* checkKind */
};

static emlrtDCInfo
    e_emlrtDCI =
        {
            30,           /* lineNo */
            34,           /* colNo */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            1                 /* checkKind */
};

static emlrtBCInfo
    g_emlrtBCI =
        {
            1,              /* iFirst */
            3810,           /* iLast */
            30,             /* lineNo */
            34,             /* colNo */
            "sat.vertices", /* aName */
            "selfShadow",   /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            0                 /* checkKind */
};

static emlrtDCInfo
    f_emlrtDCI =
        {
            31,           /* lineNo */
            34,           /* colNo */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            1                 /* checkKind */
};

static emlrtBCInfo
    h_emlrtBCI =
        {
            1,              /* iFirst */
            3810,           /* iLast */
            31,             /* lineNo */
            34,             /* colNo */
            "sat.vertices", /* aName */
            "selfShadow",   /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            0                 /* checkKind */
};

static emlrtBCInfo
    i_emlrtBCI =
        {
            -1,           /* iFirst */
            -1,           /* iLast */
            32,           /* lineNo */
            52,           /* colNo */
            "normal",     /* aName */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            0                 /* checkKind */
};

static emlrtECInfo
    emlrtECI =
        {
            -1,           /* nDims */
            43,           /* lineNo */
            1,            /* colNo */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx" /* pName */
};

static emlrtBCInfo
    j_emlrtBCI =
        {
            -1,           /* iFirst */
            -1,           /* iLast */
            40,           /* lineNo */
            13,           /* colNo */
            "tmpFlag",    /* aName */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            0                 /* checkKind */
};

static emlrtBCInfo
    k_emlrtBCI =
        {
            -1,           /* iFirst */
            -1,           /* iLast */
            21,           /* lineNo */
            41,           /* colNo */
            "faces",      /* aName */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            0                 /* checkKind */
};

static emlrtBCInfo
    l_emlrtBCI =
        {
            -1,           /* iFirst */
            -1,           /* iLast */
            22,           /* lineNo */
            36,           /* colNo */
            "faces",      /* aName */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            0                 /* checkKind */
};

static emlrtBCInfo
    m_emlrtBCI =
        {
            -1,           /* iFirst */
            -1,           /* iLast */
            23,           /* lineNo */
            36,           /* colNo */
            "faces",      /* aName */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            0                 /* checkKind */
};

static emlrtBCInfo
    n_emlrtBCI =
        {
            -1,           /* iFirst */
            -1,           /* iLast */
            29,           /* lineNo */
            45,           /* colNo */
            "faces",      /* aName */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            0                 /* checkKind */
};

static emlrtBCInfo
    o_emlrtBCI =
        {
            -1,           /* iFirst */
            -1,           /* iLast */
            30,           /* lineNo */
            40,           /* colNo */
            "faces",      /* aName */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            0                 /* checkKind */
};

static emlrtBCInfo
    p_emlrtBCI =
        {
            -1,           /* iFirst */
            -1,           /* iLast */
            31,           /* lineNo */
            40,           /* colNo */
            "faces",      /* aName */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx", /* pName */
            0                 /* checkKind */
};

static emlrtRTEInfo
    emlrtRTEI =
        {
            11,           /* lineNo */
            1,            /* colNo */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx" /* pName */
};

static emlrtRTEInfo
    b_emlrtRTEI =
        {
            12,           /* lineNo */
            1,            /* colNo */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx" /* pName */
};

static emlrtRTEInfo
    c_emlrtRTEI =
        {
            13,           /* lineNo */
            1,            /* colNo */
            "selfShadow", /* fName */
            "/Users/yyoshimula/Dropbox/MATLAB/yoshimuLibrary/object/"
            "selfShadow.mlx" /* pName */
};

/* Function Definitions */
emlrtCTX emlrtGetRootTLSGlobal(void)
{
  return emlrtRootTLSGlobal;
}

void emlrtLockerFunction(EmlrtLockeeFunction aLockee, emlrtConstCTX aTLS,
                         void *aData)
{
  omp_set_lock(&emlrtLockGlobal);
  emlrtCallLockeeFunction(aLockee, aTLS, aData);
  omp_unset_lock(&emlrtLockGlobal);
}

void selfShadow(const emlrtStack *sp, struct0_T *sat, const real_T sun[3])
{
  jmp_buf emlrtJBEnviron;
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  jmp_buf *volatile emlrtJBStack;
  emlrtStack st;
  emxArray_real_T *faces;
  emxArray_real_T *normal;
  emxArray_real_T *pos;
  real_T tmpFlag_data[7604];
  real_T vertJ[9];
  real_T b_normal[3];
  real_T b_sun[3];
  real_T posI[3];
  real_T absxk;
  real_T b_vertJ_tmp;
  real_T c_vertJ_tmp;
  real_T d;
  real_T d1;
  real_T d2;
  real_T d3;
  real_T d_vertJ_tmp;
  real_T e_vertJ_tmp;
  real_T f_vertJ_tmp;
  real_T flag;
  real_T scale;
  real_T t;
  real_T tmpFlagI;
  real_T vertJ_tmp;
  real_T y;
  real_T *faces_data;
  real_T *normal_data;
  real_T *pos_data;
  int32_T b_i;
  int32_T c_i;
  int32_T i;
  int32_T i1;
  int32_T j;
  int32_T n;
  int32_T partialTrueCount;
  int32_T selfShadow_numThreads;
  int32_T trueCount;
  int16_T tmp_data[7604];
  char_T TRANSA1;
  char_T TRANSB1;
  boolean_T b;
  boolean_T b1;
  boolean_T emlrtHadParallelError = false;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  /*  column vector */
  scale = 3.3121686421112381E-170;
  absxk = muDoubleScalarAbs(sun[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }
  absxk = muDoubleScalarAbs(sun[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }
  absxk = muDoubleScalarAbs(sun[2]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }
  y = scale * muDoubleScalarSqrt(y);
  b_sun[0] = sun[0] / y;
  b_sun[1] = sun[1] / y;
  b_sun[2] = sun[2] / y;
  /*  影判定用フラグ: 0のとき，影になっているとする */
  memset(&sat->sunlitFlag[0], 0, 7604U * sizeof(real_T));
  /*  initialization */
  /*  両方のfacetが太陽光を受けていたら計算する */
  TRANSB1 = 'N';
  TRANSA1 = 'N';
  scale = 1.0;
  absxk = 0.0;
  m_t = (ptrdiff_t)7604;
  n_t = (ptrdiff_t)1;
  k_t = (ptrdiff_t)3;
  lda_t = (ptrdiff_t)7604;
  ldb_t = (ptrdiff_t)3;
  ldc_t = (ptrdiff_t)7604;
  dgemm(&TRANSA1, &TRANSB1, &m_t, &n_t, &k_t, &scale, &sat->normal[0], &lda_t,
        &b_sun[0], &ldb_t, &absxk, &tmpFlag_data[0], &ldc_t);
  trueCount = 0;
  partialTrueCount = 0;
  for (i = 0; i < 7604; i++) {
    if (tmpFlag_data[i] > 0.0) {
      trueCount++;
      tmp_data[partialTrueCount] = (int16_T)i;
      partialTrueCount++;
    }
  }
  emxInit_real_T(sp, &faces, &emlrtRTEI);
  partialTrueCount = faces->size[0] * faces->size[1];
  faces->size[0] = trueCount;
  faces->size[1] = 4;
  emxEnsureCapacity_real_T(sp, faces, partialTrueCount, &emlrtRTEI);
  faces_data = faces->data;
  for (i = 0; i < 4; i++) {
    for (b_i = 0; b_i < trueCount; b_i++) {
      faces_data[b_i + faces->size[0] * i] =
          sat->faces[tmp_data[b_i] + 7604 * i];
    }
  }
  emxInit_real_T(sp, &normal, &b_emlrtRTEI);
  partialTrueCount = normal->size[0] * normal->size[1];
  normal->size[0] = trueCount;
  normal->size[1] = 3;
  emxEnsureCapacity_real_T(sp, normal, partialTrueCount, &b_emlrtRTEI);
  normal_data = normal->data;
  emxInit_real_T(sp, &pos, &c_emlrtRTEI);
  partialTrueCount = pos->size[0] * pos->size[1];
  pos->size[0] = trueCount;
  pos->size[1] = 3;
  emxEnsureCapacity_real_T(sp, pos, partialTrueCount, &c_emlrtRTEI);
  pos_data = pos->data;
  for (i = 0; i < 3; i++) {
    for (b_i = 0; b_i < trueCount; b_i++) {
      partialTrueCount = tmp_data[b_i] + 7604 * i;
      normal_data[b_i + normal->size[0] * i] = sat->normal[partialTrueCount];
      pos_data[b_i + pos->size[0] * i] = sat->pos[partialTrueCount];
    }
  }
  if (trueCount == 0) {
    partialTrueCount = 0;
  } else {
    partialTrueCount = muIntScalarMax_sint32(trueCount, 4);
  }
  for (i = 0; i < partialTrueCount; i++) {
    tmpFlag_data[i] = 1.0;
  }
  /*  initialization */
  emlrtEnterParallelRegion((emlrtCTX)sp, omp_in_parallel());
  emlrtPushJmpBuf((emlrtCTX)sp, &emlrtJBStack);
  selfShadow_numThreads = emlrtAllocRegionTLSs(
      sp->tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs());
#pragma omp parallel num_threads(selfShadow_numThreads) private(               \
        flag, vertJ, posI, tmpFlagI, st, emlrtJBEnviron, i1, n, j, d, d1, b,   \
            d2, d3, b1, vertJ_tmp, b_vertJ_tmp, c_vertJ_tmp, d_vertJ_tmp,      \
            e_vertJ_tmp, f_vertJ_tmp, b_normal)                                \
    firstprivate(emlrtHadParallelError)
  {
    if (setjmp(emlrtJBEnviron) == 0) {
      st.prev = sp;
      st.tls = emlrtAllocTLS((emlrtCTX)sp, omp_get_thread_num());
      st.site = NULL;
      emlrtSetJmpBuf(&st, &emlrtJBEnviron);
    } else {
      emlrtHadParallelError = true;
    }
#pragma omp for nowait
    for (c_i = 0; c_i < partialTrueCount; c_i++) {
      if (emlrtHadParallelError) {
        continue;
      }
      if (setjmp(emlrtJBEnviron) == 0) {
        tmpFlagI = 1.0;
        /*  各iごとの一時フラグ */
        i1 = pos->size[0];
        if (c_i + 1 > pos->size[0]) {
          emlrtDynamicBoundsCheckR2012b(c_i + 1, 1, pos->size[0], &emlrtBCI,
                                        &st);
        }
        posI[0] = pos_data[c_i];
        posI[1] = pos_data[c_i + pos->size[0]];
        posI[2] = pos_data[c_i + pos->size[0] * 2];
        if (faces->size[0] == 0) {
          n = 0;
        } else {
          n = muIntScalarMax_sint32(faces->size[0], 4);
        }
        for (j = 0; j < n; j++) {
          if (j != c_i) {
            if (j + 1 > i1) {
              emlrtDynamicBoundsCheckR2012b(j + 1, 1, i1, &k_emlrtBCI, &st);
            }
            d = faces_data[j];
            d1 = (int32_T)muDoubleScalarFloor(d);
            if (d != d1) {
              emlrtIntegerCheckR2012b(d, &emlrtDCI, &st);
            }
            b = (((int32_T)d < 1) || ((int32_T)d > 3810));
            if (b) {
              emlrtDynamicBoundsCheckR2012b((int32_T)d, 1, 3810, &b_emlrtBCI,
                                            &st);
            }
            if (j + 1 > i1) {
              emlrtDynamicBoundsCheckR2012b(j + 1, 1, i1, &l_emlrtBCI, &st);
            }
            flag = faces_data[j + faces->size[0]];
            if (flag != (int32_T)muDoubleScalarFloor(flag)) {
              emlrtIntegerCheckR2012b(flag, &b_emlrtDCI, &st);
            }
            if (((int32_T)flag < 1) || ((int32_T)flag > 3810)) {
              emlrtDynamicBoundsCheckR2012b((int32_T)flag, 1, 3810, &c_emlrtBCI,
                                            &st);
            }
            if (j + 1 > i1) {
              emlrtDynamicBoundsCheckR2012b(j + 1, 1, i1, &m_emlrtBCI, &st);
            }
            d2 = faces_data[j + faces->size[0] * 2];
            d3 = (int32_T)muDoubleScalarFloor(d2);
            if (d2 != d3) {
              emlrtIntegerCheckR2012b(d2, &c_emlrtDCI, &st);
            }
            b1 = (((int32_T)d2 < 1) || ((int32_T)d2 > 3810));
            if (b1) {
              emlrtDynamicBoundsCheckR2012b((int32_T)d2, 1, 3810, &d_emlrtBCI,
                                            &st);
            }
            vertJ_tmp = sat->vertices[(int32_T)d - 1];
            vertJ[0] = vertJ_tmp;
            vertJ[1] = sat->vertices[(int32_T)flag - 1];
            b_vertJ_tmp = sat->vertices[(int32_T)d2 - 1];
            vertJ[2] = b_vertJ_tmp;
            c_vertJ_tmp = sat->vertices[(int32_T)d + 3809];
            vertJ[3] = c_vertJ_tmp;
            vertJ[4] = sat->vertices[(int32_T)flag + 3809];
            d_vertJ_tmp = sat->vertices[(int32_T)d2 + 3809];
            vertJ[5] = d_vertJ_tmp;
            e_vertJ_tmp = sat->vertices[(int32_T)d + 7619];
            vertJ[6] = e_vertJ_tmp;
            vertJ[7] = sat->vertices[(int32_T)flag + 7619];
            f_vertJ_tmp = sat->vertices[(int32_T)d2 + 7619];
            vertJ[8] = f_vertJ_tmp;
            /*  3x3, coordinate of j-th facet */
            if (j + 1 > i1) {
              emlrtDynamicBoundsCheckR2012b(j + 1, 1, i1, &e_emlrtBCI, &st);
            }
            b_normal[0] = normal_data[j];
            b_normal[1] = normal_data[j + normal->size[0]];
            b_normal[2] = normal_data[j + normal->size[0] * 2];
            flag = calcSelfShadow(b_sun, b_normal, vertJ, posI);
            tmpFlagI *= flag;
            /*  OR演算っぽく */
            /*  四角facetでは2回計算 */
            if (!muDoubleScalarIsNaN(faces_data[faces->size[0] * 3])) {
              if (j + 1 > i1) {
                emlrtDynamicBoundsCheckR2012b(j + 1, 1, i1, &n_emlrtBCI, &st);
              }
              if (d2 != d3) {
                emlrtIntegerCheckR2012b(d2, &d_emlrtDCI, &st);
              }
              if (b1) {
                emlrtDynamicBoundsCheckR2012b((int32_T)d2, 1, 3810, &f_emlrtBCI,
                                              &st);
              }
              if (j + 1 > i1) {
                emlrtDynamicBoundsCheckR2012b(j + 1, 1, i1, &o_emlrtBCI, &st);
              }
              flag = faces_data[j + faces->size[0] * 3];
              if (flag != (int32_T)muDoubleScalarFloor(flag)) {
                emlrtIntegerCheckR2012b(flag, &e_emlrtDCI, &st);
              }
              if (((int32_T)flag < 1) || ((int32_T)flag > 3810)) {
                emlrtDynamicBoundsCheckR2012b((int32_T)flag, 1, 3810,
                                              &g_emlrtBCI, &st);
              }
              if (j + 1 > i1) {
                emlrtDynamicBoundsCheckR2012b(j + 1, 1, i1, &p_emlrtBCI, &st);
              }
              if (d != d1) {
                emlrtIntegerCheckR2012b(d, &f_emlrtDCI, &st);
              }
              if (b) {
                emlrtDynamicBoundsCheckR2012b((int32_T)d, 1, 3810, &h_emlrtBCI,
                                              &st);
              }
              vertJ[0] = b_vertJ_tmp;
              vertJ[1] = sat->vertices[(int32_T)flag - 1];
              vertJ[2] = vertJ_tmp;
              vertJ[3] = d_vertJ_tmp;
              vertJ[4] = sat->vertices[(int32_T)flag + 3809];
              vertJ[5] = c_vertJ_tmp;
              vertJ[6] = f_vertJ_tmp;
              vertJ[7] = sat->vertices[(int32_T)flag + 7619];
              vertJ[8] = e_vertJ_tmp;
              /*  3x3, coordinate of j-th facet */
              if (j + 1 > i1) {
                emlrtDynamicBoundsCheckR2012b(j + 1, 1, i1, &i_emlrtBCI, &st);
              }
              b_normal[0] = normal_data[j];
              b_normal[1] = normal_data[j + normal->size[0]];
              b_normal[2] = normal_data[j + normal->size[0] * 2];
              flag = calcSelfShadow(b_sun, b_normal, vertJ, posI);
              tmpFlagI *= flag;
              /*  OR演算っぽく */
            }
          } else {
            /*  do nothing when i == j */
          }
          if (*emlrtBreakCheckR2012bFlagVar != 0) {
            emlrtBreakCheckR2012b(&st);
          }
        }
        if (c_i + 1 > partialTrueCount) {
          emlrtDynamicBoundsCheckR2012b(c_i + 1, 1, partialTrueCount,
                                        &j_emlrtBCI, &st);
        }
        tmpFlag_data[c_i] = tmpFlagI;
        /*  結果を保存 */
        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b(&st);
        }
      } else {
        emlrtHadParallelError = true;
      }
    }
  }
  emlrtPopJmpBuf((emlrtCTX)sp, &emlrtJBStack);
  emlrtExitParallelRegion((emlrtCTX)sp, omp_in_parallel());
  emxFree_real_T(sp, &pos);
  emxFree_real_T(sp, &normal);
  emxFree_real_T(sp, &faces);
  if (trueCount != partialTrueCount) {
    emlrtSubAssignSizeCheck1dR2017a(trueCount, partialTrueCount, &emlrtECI,
                                    (emlrtConstCTX)sp);
  }
  for (i = 0; i < partialTrueCount; i++) {
    sat->sunlitFlag[tmp_data[i]] = tmpFlag_data[i];
  }
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (selfShadow.c) */
