/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: selfShadow_.c
 *
 * MATLAB Coder version            : 24.2
 * C/C++ source code generated on  : 2024/11/21 17:17:03
 */

/* Include Files */
#include "selfShadow_.h"
#include "calcSelfShadow_.h"
#include "calcSelfShadow__types.h"
#include "mtimes.h"
#include "norm.h"
#include <string.h>

/* Function Definitions */
/*
 * Arguments    : struct0_T *sat
 *                double sun[3]
 * Return Type  : void
 */
void selfShadow_(struct0_T *sat, double sun[3])
{
  double tmpFlag_data[1000];
  double b_sun[3];
  double tmpFlagI;
  int i;
  int j;
  int tmpFlag_tmp_tmp;
  int trueCount;
  int varargin_1;
  int y;
  short tmp_data[1000];
  /*  calculate self-shadowing */
  /*  inputs */
  /*  note */
  /*  全componentのメッシュに対して処理を行う */
  /*  影の処理関数(番目のfacetが影を作る方，番目のfacetは影が写る方) */
  /*  references */
  /*  NA */
  /*  revisions */
  /*  20200811  y.yoshimura, y.yoshimula@gmail.com */
  /*  See also q2zyx. */
  sun[2] = 3.0;
  /*  column vector */
  tmpFlagI = b_norm(sun);
  b_sun[0] = sun[0] / tmpFlagI;
  b_sun[1] = sun[1] / tmpFlagI;
  b_sun[2] = 3.0 / tmpFlagI;
  /*  影判定用フラグ: 0のとき，影になっているとする */
  varargin_1 = sat->faces.size[0];
  if (varargin_1 >= 3) {
    y = varargin_1;
  } else {
    y = 3;
  }
  if (varargin_1 == 0) {
    varargin_1 = 0;
  } else {
    varargin_1 = y;
  }
  sat->shadowFlag.size[0] = varargin_1;
  if (varargin_1 - 1 >= 0) {
    memset(&sat->shadowFlag.data[0], 0,
           (unsigned int)varargin_1 * sizeof(double));
  }
  /*  initialization */
  /*  両方のfacetが太陽光を受けていたら計算する */
  varargin_1 = mtimes(sat->normal.data, sat->normal.size, b_sun, tmpFlag_data);
  trueCount = 0;
  y = 0;
  for (i = 0; i < varargin_1; i++) {
    if (tmpFlag_data[i] > 0.0) {
      trueCount++;
      tmp_data[y] = (short)i;
      y++;
    }
  }
  if (trueCount >= 3) {
    y = trueCount;
  } else {
    y = 3;
  }
  if (trueCount == 0) {
    tmpFlag_tmp_tmp = 0;
  } else {
    tmpFlag_tmp_tmp = y;
  }
  for (varargin_1 = 0; varargin_1 < tmpFlag_tmp_tmp; varargin_1++) {
    tmpFlag_data[varargin_1] = 1.0;
  }
  /*  initialization */
  /*  use "for", if parallel toolbox is not installed */
  for (i = 0; i < tmpFlag_tmp_tmp; i++) {
    double posI[3];
    short b_i;
    tmpFlagI = 1.0;
    /*  各iごとの一時フラグ */
    b_i = tmp_data[i];
    posI[0] = sat->pos.data[b_i];
    posI[1] = sat->pos.data[b_i + sat->pos.size[0]];
    posI[2] = sat->pos.data[b_i + sat->pos.size[0] * 2];
    for (j = 0; j < tmpFlag_tmp_tmp; j++) {
      if (j != i) {
        double vertJ[9];
        double b_sat[3];
        double flag;
        b_i = tmp_data[j];
        varargin_1 = (int)sat->faces.data[b_i];
        y = (int)sat->faces.data[b_i + sat->faces.size[0]];
        trueCount = (int)sat->faces.data[b_i + sat->faces.size[0] * 2];
        /*  3x3, coordinate of j-th facet             */
        vertJ[0] = sat->vertices[varargin_1 - 1];
        vertJ[1] = sat->vertices[y - 1];
        vertJ[2] = sat->vertices[trueCount - 1];
        b_sat[0] = sat->normal.data[b_i];
        vertJ[3] = sat->vertices[varargin_1 + 2999];
        vertJ[4] = sat->vertices[y + 2999];
        vertJ[5] = sat->vertices[trueCount + 2999];
        b_sat[1] = sat->normal.data[b_i + sat->normal.size[0]];
        vertJ[6] = sat->vertices[varargin_1 + 5999];
        vertJ[7] = sat->vertices[y + 5999];
        vertJ[8] = sat->vertices[trueCount + 5999];
        b_sat[2] = sat->normal.data[b_i + sat->normal.size[0] * 2];
        flag = calcSelfShadow_(b_sun, b_sat, vertJ, posI);
        tmpFlagI *= flag;
        /*  OR演算っぽく */
      } else {
        /*  do nothing when i == j */
      }
    }
    tmpFlag_data[i] = tmpFlagI;
    /*  結果を保存 */
  }
  for (varargin_1 = 0; varargin_1 < tmpFlag_tmp_tmp; varargin_1++) {
    sat->shadowFlag.data[tmp_data[varargin_1]] = tmpFlag_data[varargin_1];
  }
}

/*
 * File trailer for selfShadow_.c
 *
 * [EOF]
 */
