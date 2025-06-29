
/** ssdl.h ********************************************************************
 *                                                                            *
 *  coding:FEB/06/1997/1.00: 1st coding by T.HANADA                           *
 *  update:DEC/18/2013/1.0A: Proto type for reduce.c added                    *
 *  update:___/__/____/_.__:                                                  *
 *  bugfix:___/__/____/_.__:                                                  *
 *                                                                            *
 ******************************************************************************/

#ifndef __SSDL_H__
#define __SSDL_H__

/*******************************  header files  *******************************/

#include <math.h>
#include <stdio.h>
#include <time.h>

#ifndef __SSRCSV_H__
#include <AstroLib/ssrcsv.h>
#endif

#ifndef __TLESET_H__
#include <AstroLib/tleset.h>
#endif

/***************************  constants and macros  ***************************/

#define GE 2975536351779840.0                     /* GE, in km^3/day^2. ----- */
#define GP                0.997269566329084       /* GEO period, in days. --- */
#define GA            42164.16962356673645807092  /* GEO altitude, in km. --- */

#define BSTR2AOMCD       12.741620760032537091964164406946
#define AOMCD2BSTR        0.078482951175          /* BC = 1/(12.741621 x B*)  */

/****************************  function prototypes  ***************************/

/* FK5Reduction.c */

void FK5ToMoD (double, double *, double *, double *);
void MoDToToD (double, double *, double *, double *);
void FK5ToToD (double, double *, double *, double *);
void ToDToMoD (double, double *, double *, double *);
void MoDToFK5 (double, double *, double *, double *);
void ToDToFK5 (double, double *, double *, double *);
void Rigorous (double, double, double *, double *);
void EfNutatn (double, double *, double *);

/* TEME.c */

void MoDToTEME (double, double *, double *, double *);
void FK5ToTEME (double, double *, double *, double *);
void TEMEToMoD (double, double *, double *, double *);
void TEMEToFK5 (double, double *, double *, double *);

/* getle.c */

char *get_tle   (FILE *, tleset *, int *);
char *get_tlewn (FILE *, tleset *);
char *get_tlewo (FILE *, tleset *);
tleset *readtle (FILE *, int *, int );

/* getss.c */

char *get_objects_in_orbit (FILE *, ssrcsv *, int *);
char *get_objects_no_longer_in_orbit (FILE *, ssrcsv *, int *);

/* jkwok.c */

void MeanToOsculating (double *, double *, double *, double *, double *, double *, double *, double *);
int OsculatingToMean (double *, double *, double *, double *, double *, double *, double *, double *);

/* spmtn.c */

void OrbitToSphericalMotion (double, double, double, double, double, double, double, double *, double *, double *, double *, double *, double *, double *);

/* store.c */

void storeGabbardDiag (FILE *, double, double, double, double);
void storeInitCond (FILE *, double, double, double, double);
void storeTwoLineElement (FILE *, char *, int, char *, double, double, double, double, double, double, double);
char *writetle (FILE *, tleset *, int);

/* times.c */

clock_t clockm(void);
struct tm *get_UTC(long *, int *, double *, int  *, int *, double *);
struct tm *get_JST(long *, int *, double *, int  *, int *, double *);
double GMST0 (double );

/* odam.c */

double amach ();
double sign (double , double );
void odam (double *, double *, void (*)(), int , double , int *, double *, double *, double *, int *, int *);
void ude (void (*)(), int , double *, double *, double , double *, double *, int *, int *, double *, int *, int , int , int , int , int , int , int *, int *, int *);
void uint1 (double *, double *, double , double *, double *, int , int *, double *, double *);
void uste1 (double *, double *, void (*)(), int , double *, double *, double *, int *, double *, double *, double *, int *, int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *);

#endif /* !__SSDL_H__ */
