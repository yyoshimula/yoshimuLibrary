
/** getle.c *******************************************************************
 *                                                                            *
 *  coding:1996/FEB/06/1.00: 1st coding by T.HANADA                           *
 *  update:1998/MAR/11/1.20: Not knowing the two line elements w/ or w/o a    *
 *                           24-characters name, the function get_tle() is    *
 *                           available to know it with the (flag).            *
 *  bugfix:1998/AUG/18/1.21: The bug in interpreting the epoch time was       *
 *                           fixed.  An epoch of 98001.00000000 corressponds  *
 *                           to 0000 UT on 1998 January 01.                   *
 *  update:2000/JAN/22/1.22: Added `orgmma' function that recovers the        *
 *                           original mean motion from the inputs elements    *
 *                           by invertingthe SGP4 formula.                    *
 *  update:2003/MAY/02/1.23: Modified the function get_tle() to read (bstr).  *
 *  bugfix:2003/MAY/16/1.24: The bug in the function meanmm() was fixed.      *
 *  update:2012/FEB/22/1.30: MeanToOsculating() added.  readtle() deleted.    *
 *  update:2012/FEB/22/1.3A: OsculatingToMean() added.  orgmma() and meanmm() *
 *                           deleted.                                         *
 *  update:2012/JUL/31/1.3B: readtle() restored by using MeanToOsculating().  *
 *  update:2012/AUG/20/1.3C: Consitency of astronomical constants improved.   *
 *  bugfix:2013/MAY/17/1.3D: MeanToOsculating() and OsculatingToMean()        *
 *                           confirmed and corrected.                         *
 *  bugfix:2013/JUL/27/1.3E: Calculation of semimajor axis in readtle()       *
 *                           fixed.                                           *
 *  update:2013/DEC/18/1.40: MeanToOsculating() and OsculatingToMean()        *
 *                           removed to be filed as jkwok.c.                  *
 *  bugfix:2014/JUN/30/1.41: get_tle(), get_tlewn(), and get_tlewo() modified *
 *                           to correctly scan old ID such as "68081  E".     *
 *  bugfix:2016/MAR/24/1.42: get_tle(), and get_tlewo() modified to correctly *
 *                           scan new OBJECT_NAME such as "0 VANGUARD 1".     *
 *  update:____/___/__/_.__:                                                  *
 *  bugfix:____/___/__/_.__:                                                  *
 *                                                                            *
 ******************************************************************************/
#include <math.h>          /* maths functions */
#include "getle.h"
#include "mex.h"

#ifndef TRUE
#include <curses.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <AstroLib.h>
#include <ssdl.h>
#include <float.h>

#define GE 2975536351779840.0                     /* GE, in km^3/day^2. ----- */
#define XKE                .74366916E-1           /* (er/min)^(3/2) --------- */
#define CK2               5.413080E-4
#define DE2RA              .174532925E-1
#define TWOPI             6.2831853
#define XMNPDA         1440.
#define XKMPER         6378.135
#define _MAXNUM_       1000
#ifndef FLT_EPSILON
#define FLT_EPSILON       1.0e-5F
#endif
#ifndef DBL_EPSILON
#define DBL_EPSILON       1.0e-9
#endif

#define max(a,b) ((a) >= (b) ? (a) : (b))

double amach ();

/** get_tle *******************************************************************

  DESCRIPTION: Reads two line elements from a given stream associated with a 
               two line elements set file.

     (stream): A stream associated with the two line elements set file.
       (data): Receives the two line elements w/ or w/o a 24-characters name.
       (flag): Indicates whether (data) comes w/ or w/o a 24-characters name.
               flag = 1 (TRUE) name avairable
               flag = 0 (FALSE) name not avairable

      RETURNS: Error code.

 ******************************************************************************/

char *get_tle (FILE *stream, tleset *data, int *flag)
{
  char *buf, c;
  double a;
  long year;

  if (fgets(data->one, sizeof(data->one), stream) == NULL)
    return (char *) NULL;
  *strstr(data->one, "\n") = '\0';
  buf = data->one;
  switch (buf[0]) {
    case '1':
      *flag = FALSE;
      break;
    case '0':
      buf += 2;
    default:
      *flag = TRUE;
      strcpy(data->name, buf); //puts(data->name);
      if (fgets(data->one, sizeof(data->one), stream) == NULL)
        return (char *) NULL;
      *strstr(data->one, "\n") = '\0';
      buf = data->one;
  }
  sscanf(&buf[2], "%5d", &data->num);
  sscanf(&buf[9], "%8c", data->ID); data->ID[8] = '\0';
  c = buf[20]; buf[20] = '\0'; sscanf(&buf[18], "%ld", &year); buf[20] = c;
  if (year > 56L) 
    year += 1900L;
  else
    year += 2000L;
  sscanf(&buf[20], "%lf", &data->epoch);
  data->epoch += DateTimeToJD(year, 1, 0., 0, 0, 0., 0., TRUE);
//sscanf(&buf[33], "%lf", &data->ddm);
//sscanf(&buf[44], "%lf%lf", &data->dddm, &a); data->dddm *= pow(10.0, a - 5.0);
  sscanf(&buf[54], "%lf%lf", &data->bstr, &a); data->bstr *= pow(10.0, a - 5.0);
  if (fgets(data->two, sizeof(data->two), stream) == NULL)
    return (char *) NULL;
  *strstr(data->two, "\n") = '\0';
  buf = data->two;
  sscanf(&buf[8], "%lf", &data->i);
  sscanf(&buf[17], "%lf", &data->n);
  sscanf(&buf[26], "%lf", &data->e); data->e /= 10000000.0;
  sscanf(&buf[34], "%lf", &data->w);
  sscanf(&buf[43], "%lf", &data->m);
  c = buf[63]; buf[63] = '\0'; sscanf(&buf[52], "%11lf", &data->dm); buf[63] = c;

  return (char *) EOF;
}

/** get_tlewn *****************************************************************

  DESCRIPTION: Reads two line elements from a given stream associated with a 
               two line elements set file.

     (stream): A stream associated with the two line elements set file.
       (data): Receives the two line elements with a 24-characters name.

      RETURNS: Error code.

 ******************************************************************************/

char *get_tlewn (FILE *stream, tleset *data)
{
  char *buf, c;
  double a;
  long year;
/*
 Get Line 0 (OBJECT_NAME). -------------------------------------------------- */
  
  if (fgets(data->one, sizeof(data->one), stream) == NULL)
    return (char *) NULL;
  *strstr(data->one, "\n") = '\0';
  buf = data->one;
  if (data->one[0] == '0') buf += 2;
  strcpy(data->name, buf);
/*
  Get Line 1. --------------------------------------------------------------- */

  if (fgets(data->one, sizeof(data->one), stream) == NULL)
    return (char *) NULL;
  *strstr(data->one, "\n") = '\0';
  buf = data->one;
  sscanf(&buf[2], "%5d", &data->num);
  sscanf(&buf[9], "%8c", data->ID); data->ID[8] = '\0';
  c = buf[20]; buf[20] = '\0'; sscanf(&buf[18], "%ld", &year); buf[20] = c;
  if (year > 56L)
    year += 1900L;
  else
    year += 2000L;
  sscanf(&buf[20], "%lf", &data->epoch);
  data->epoch += DateTimeToJD(year, 1, 0., 0, 0, 0., 0., TRUE);
//sscanf(&buf[33], "%lf", &data->ddm);
//sscanf(&buf[44], "%lf%lf", &data->dddm, &a); data->dddm *= pow (10.0, a - 5.0);
  sscanf(&buf[54], "%lf%lf", &data->bstr, &a); data->bstr *= pow (10.0, a - 5.0);
/*
  Get Line 2. --------------------------------------------------------------- */

  if (fgets(data->two, sizeof(data->two), stream) == NULL)
    return (char *) NULL;
  *strstr(data->two, "\n") = '\0';
  buf = data->two;
  sscanf(&buf[8], "%lf", &data->i);
  sscanf(&buf[17], "%lf", &data->n);
  sscanf(&buf[26], "%lf", &data->e); data->e /= 10000000.0;
  sscanf(&buf[34], "%lf", &data->w);
  sscanf(&buf[43], "%lf", &data->m);
  c = buf[63]; buf[63] = '\0'; sscanf(&buf[52], "%11lf", &data->dm); buf[63] = c;

  return (char *) EOF;
}

/** get_tlewo *****************************************************************

  DESCRIPTION: Reads two line elements from a given stream associated with a 
               two line elements set file.

     (stream): A stream associated with the two line elements set file.
       (data): Receives the two line elements without a 24-characters name.

      RETURNS: Error code.

 ******************************************************************************/

char *get_tlewo (FILE *stream, tleset *data)
{
  char *buf, c;
  double a;
  long year;

/*
  Get Line 1. --------------------------------------------------------------- */

  if (fgets(data->one, sizeof(data->one), stream) == NULL)
    return (char *) NULL;
  *strstr(data->one, "\n") = '\0';
  buf = data->one;
  sscanf(&buf[2], "%5d", &data->num);
  sscanf(&buf[9], "%8c", data->ID); data->ID[8] = '\0';
  c = buf[20]; buf[20] = '\0'; sscanf(&buf[18], "%ld", &year); buf[20] = c;
  if (year > 56L)
    year += 1900L;
  else
    year += 2000L;
  sscanf(&buf[20], "%lf", &data->epoch);
  data->epoch += DateTimeToJD(year, 1, 0., 0, 0, 0., 0., TRUE);
//sscanf(&buf[33], "%lf", &data->ddm);
//sscanf(&buf[44], "%lf%lf", &data->dddm, &a); data->dddm *= pow(10.0, a - 5.0);
  sscanf(&buf[54], "%lf%lf", &data->bstr, &a); data->bstr *= pow(10.0, a - 5.0);
/*
  Get Line 2. --------------------------------------------------------------- */

  if (fgets(data->two, sizeof(data->two), stream) == NULL)
    return (char *) NULL;
  *strstr(data->two, "\n") = '\0';
  buf = data->two;
  sscanf(&buf[8], "%lf", &data->i);
  sscanf(&buf[17], "%lf", &data->n);
  sscanf(&buf[26], "%lf", &data->e); data->e /= 10000000.0;
  sscanf(&buf[34], "%lf", &data->w);
  sscanf(&buf[43], "%lf", &data->m);
  c = buf[63]; buf[63] = '\0'; sscanf(&buf[52], "%11lf", &data->dm); buf[63] = c;

  return (char *) EOF;
}

/** readtle *******************************************************************
 
  DESCRIPTION: Reads two line elements from a given stream associated with a 
               two line elements set file.
 
     (stream): Stream associated with two line elements sets file.
        (num): Number of two line elements sets.
      (norad): If (norad) then this function converts "mean" to "osculating" 
               elements, considering J_2 perturbation.  
      RETURNS: A corresponding tleset pointer.
 
        NOTES: Receives mean mortion in radians per day, not in revolutions per
               day.  Receives angles (.i, .n, .w, .m) in radians.  
 
 ******************************************************************************/

tleset *readtle (FILE *stream, int *num, int norad)
{
  char buf[BUFSIZ];
  double f;
  int i, flag;
  tleset *data, tmp;
  
  if (stream == NULL)
    return (tleset *) NULL;
  get_tle(stream, &tmp, &flag);
  rewind(stream);
  for (i = 0; fgets(buf, sizeof (buf), stream) != NULL; i++) ;
  rewind(stream);
  *num = i / (2 + flag);
  if ((data = (tleset *) calloc(*num, sizeof(tleset))) == NULL) {
    perror(""); return (tleset *) NULL;
  }
  for (i = 0; i < *num; i++) {
    get_tle(stream, &data[i], &flag);
    if (norad) {
      MeanToOsculating(&data[i].dm, &data[i].e, &data[i].i,
                       &data[i].n, &data[i].w, &data[i].m, &data[i].a, &f);
      data[i].dm *= TWOPI;
    }
    else {
      data[i].dm *= TWO_PI;
      data[i].a = pow(GE / (data[i].dm * data[i].dm), 1.0 / 3.0);
      data[i].i *= RAD_PER_DEG;
      data[i].n *= RAD_PER_DEG;
      data[i].w *= RAD_PER_DEG;
      data[i].m *= RAD_PER_DEG;
    }
    data[i].q = data[i].a * (1.0 - data[i].e);
  }
  return (tleset *) data;
}


/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
/* variable declarations here */
    FILE *stream;              /* input scalar */
    tleset *data;              /* input scalar */
    int *flag;              /* output vector 3x1 */
    
        /* check for proper number of arguments */
    if(nrhs != 3) {
        mexErrMsgIdAndTxt("get TLE:nrhs","Three inputs required.");
    }
    if(nlhs != 1) {
        mexErrMsgIdAndTxt("get TLE:nlhs","One output required.");
    }

      /* make sure the input arguments are scalar
    if( !mxIsDoublze(prhs[0]) || 
         mxIsComplex(prhs[0]) ||
         mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("IGRF12:nrhs","Input must be a scalar.");
    }
        if( !mxIsDouble(prhs[1]) || 
         mxIsComplex(prhs[1]) ||
         mxGetNumberOfElements(prhs[1])!=1 ) {
        mexErrMsgIdAndTxt("IGRF12:nrhs","Input must be a scalar.");
    }
        if( !mxIsDouble(prhs[2]) || 
         mxIsComplex(prhs[2]) ||
         mxGetNumberOfElements(prhs[2])!=1 ) {
        mexErrMsgIdAndTxt("IGRF12:nrhs","Input must be a scalar.");
    }
        if( !mxIsDouble(prhs[3]) || 
         mxIsComplex(prhs[3]) ||
         mxGetNumberOfElements(prhs[3])!=1 ) {
        mexErrMsgIdAndTxt("IGRF12:nrhs","Input must be a scalar.");
    } */
    
    data = mxGetVector(prhs[0]);
    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1, 3,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    b = mxGetPr(plhs[0]);
    

    /* call the computational routine */
    getle(stream, data, flag);
  
}