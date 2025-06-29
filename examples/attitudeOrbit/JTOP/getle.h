#ifndef	getle_include_
#define	getle_include_

//char *get_tle (FILE *stream, tleset *data, int *flag)
//void igrf12syn(double date, double alt, double colat, double elong, double *x, double *y, double *z);
void get_tle(FILE *, tleset *, int *);
#endif