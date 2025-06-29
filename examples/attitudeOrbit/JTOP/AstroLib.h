
/*****************************  COPYRIGHT NOTICE  *****************************

  Copyright 1994 by Tim DeBenedictis, 306 Kensington Pl., Syracuse NY 13210.

  The user is hereby granted permission to make copies of this source file
  for his/her personal use, and may utilize it for non-commercial projects
  and products only.  GENERAL REDISTRIBUTION OF THIS SOURCE FILE IS PROHIBITED.
  IT MAY NOT BE INCORPORATED INTO ANY COMMERCIAL PRODUCT FOR RESALE WITHOUT
  SPECIFIC PERMISSION OF THE AUTHOR.  To obtain an authorization agreement
  for commerical reuse and sale, contact the author at the above address.

 ******************************************************************************/

/******************************************************************************
 *                                                                            *
 *  update:____/___/__/_.__:                                                  *
 *  bugfix:____/___/__/_.__:                                                  *
 *                                                                            *
 ******************************************************************************/

#ifndef __ASTROLIB_H__
#define __ASTROLIB_H__

#include <math.h>
#include <time.h>

/***************************  constants and macros  ***************************/

#define PI                   M_PI
#define TWO_PI               (2.0 * M_PI)
#define HALF_PI              M_PI_2
#define RAD_PER_DEG          (M_PI / 180.0)
#define DEG_PER_RAD        (180.0 * M_1_PI)
#define HOUR_PER_RAD        (12.0 * M_1_PI)
#define ARCSEC_PER_RAD  (648000.0 * M_1_PI)

#define B1900                2415020.31352
#define B1950                2433282.423
#define J2000                2451545.0
#define DAYS_PER_JULIAN_YEAR 365.25
#define DAYS_PER_BESSEL_YEAR 365.242198781
#define SIDEREAL_PER_SOLAR   1.002737909350795
#define SOLAR_PER_SIDEREAL   0.997269566329084

#define HELIO_GAUSS_CONST     0.01720209895
#define KM_PER_AU             149597870.0
#define LIGHT_DAYS_PER_AU     0.005775518

#define GEO_GAUSS_CONST       0.01743669161
#define EARTHS_PER_SOLAR_MASS 332946.0
#define KM_PER_EARTH_RADII    6378.1363
#define EARTH_FLATTENING      (1.0/298.257)
#define EARTH_J2              0.1082626925638815E-2
#define EARTH_J3             -0.2532307818191774E-5
#define EARTH_J4             -0.1620429990000000E-5
#define EARTH_C21            -0.2414000052222093E-9
#define EARTH_C22             0.1574421758350994E-5
#define EARTH_C31             0.2190922081404716E-5
#define EARTH_C32             0.3089143533816488E-6
#define EARTH_C33             0.1005601040626586E-6
#define EARTH_C41            -0.5088433157745930E-6
#define EARTH_C42             0.7834048953908266E-7
#define EARTH_C43             0.5917924178248455E-7
#define EARTH_C44            -0.3982546443559900E-8
#define EARTH_S21             0.1543099973784379E-8
#define EARTH_S22            -0.9037666669616874E-6
#define EARTH_S31             0.2687418863136855E-6
#define EARTH_S32            -0.2115075122835371E-6
#define EARTH_S33             0.1971780250456937E-6
#define EARTH_S41            -0.4491281704606470E-6
#define EARTH_S42             0.1482219920570510E-6
#define EARTH_S43            -0.1201263975958658E-7
#define EARTH_S44             0.6525548406274755E-8

#define GALACTIC_N_POLE_RA   192.25
#define GALACTIC_N_POLE_DEC   27.4
#define GALACTIC_LON_ASC_NODE 33.0

/***************************  structs and typedefs  ***************************

#define MAX_NUM_ARGS 8

typedef struct VFPTerm
{
  float amp;
  char tcoeff;
  char argfunc;
  char argcoeffs[MAX_NUM_ARGS];
}
VFPTerm;
 */
/****************************  function prototypes  ***************************/

/* Time.c */

double DateTimeToJD(long, int, double, int, int, double, double, int);
void JDToDateTime(double, double, long *, int *, double *, int *, int *, double *, int);
double EpochToJD(double, char);
double JDToEpoch(double, char);
double GMST(double);
double SemiDiurnalArc(double, double, double);
double UTToTDT(double);

/* Angle.c */

double DegMinSecToDecimal(int, int, double, char);
void DecimalToDegMinSec(double, int *, int *, double *, char *);
double Mod2Pi(double);
double Separation(double, double, double, double);
double PositionAngle(double, double, double, double);

/* Vector.c */

void VectorSum(double[], double[], double[]);
void VectorDifference(double[], double[], double[]);
double VectorMagnitude(double[]);
double DotProduct(double[], double[]);
void CrossProduct(double[], double[], double[]);
void SphericalToXYZ(double, double, double, double *, double *, double *);
void XYZToSpherical(double, double, double, double *, double *, double *);
void SphericalToXYZMotion(double, double, double, double, double, double, double *, double *, double *, double *, double *, double *);
void XYZToSphericalMotion(double, double, double, double, double, double, double *, double *, double *, double *, double *, double *);

/* Reduce.c */

void Precession(double, double, double *, double *, double *, double *, double *, double *);
double Obliquity(double);
void Nutation(double, double *, double *);
void EarthVelocity(double, double *, double *, double *);
void Aberration(double, double, double, double, double, double, double *, double *, double *, double *);
void GeodeticToGeocentricXYZ(double, double, double, double, double, double *, double *, double *);
void GeocentricXYZToGeodetic(double, double, double, double, double, double *, double *, double *);
double Refraction(double, double, double);

/* Rotate.c */

void RotateVector(double[], double[]);
void RotateMatrix(double[], double[]);
void UnRotateVector(double[], double[]);
void UnRotateMatrix(double[], double[]);
void TransposeMatrix(double[]);
void SetIdentityMatrix(double[]);
void SetRotationMatrix(double[], int, ...);
void SetEclipticMatrix(double[], double, char);
void SetHorizonMatrix(double[], double, double, char);
void SetGalacticMatrix(double[], double, char);
void SetPrecessionMatrix(double[], double, double, char);
void SetEclipticPrecessionMatrix(double[], double, double);
void SetNutationMatrix(double[], double, double, double, char);
void RotateOrbit(double[], double *, double *, double *);

/* Orbit.c */

double MeanMotion(double, double, double);
void J2MeanMotion(double, double, double, double, double, double, double *, double *, double *);
void SolveKeplersEqn(double, double, double, double *, double *);
void OrbitToSpherical(double, double, double, double, double, double, double *, double *, double *);
void OrbitToXYZ(double, double, double, double, double, double, double *, double *, double *);
void OrbitToXYZMotion(double, double, double, double, double, double, double, double *, double *, double *, double *, double *, double *, double *);
void XYZMotionToOrbit(double, double, double, double, double, double, double, double *, double *, double *, double *, double *, double *, double *);

/* VFPlanet.c

double SumVFPTerms(double[], double, VFPTerm[], short, short);
void Sun(double, double *, double *, double *);
void Moon(double, double *, double *, double *);
void Mercury(double, double *, double *, double *);
void Venus(double, double *, double *, double *);
void Mars(double, double *, double *, double *);
void Jupiter(double, double *, double *, double *);
void Saturn(double, double *, double *, double *);
void Uranus(double, double *, double *, double *);
void Neptune(double, double *, double *, double *);
void Pluto(double, double *, double *, double *);
 */
/* VSOP87.c */

void Sun (double , double *, double *, double *);
void Mercury (double , double *, double *, double *);
void Venus (double , double *, double *, double *);
void Earth (double , double *, double *, double *);
void Mars (double , double *, double *, double *);
void Jupiter (double , double *, double *, double *);
void Saturn (double , double *, double *, double *);
void Uranus (double , double *, double *, double *);
void Neptune (double , double *, double *, double *);
void Pputo (double , double *, double *, double *);

/* ELP2K.c */

void Moon (double , double *, double *, double *);

/* Physical.c */

void PlanetographicCoordinates(double, double, double, double, double, double *, double *, double *);
void SunRotation(double, double *, double *, double *);
void MercuryRotation(double, double *, double *, double *);
void VenusRotation(double, double *, double *, double *);
void EarthRotation(double, double *, double *, double *);
void MarsRotation(double, double *, double *, double *);
void JupiterRotation(double, double *, double *, double *, char);
void SaturnRotation(double, double *, double *, double *, char);
void UranusRotation(double, double *, double *, double *);
void NeptuneRotation(double, double *, double *, double *);
void PlutoRotation(double, double *, double *, double *);
void MoonRotation(double, double *, double *, double *);
double MercuryMagnitude(double, double, double);
double VenusMagnitude(double, double, double);
double MarsMagnitude(double, double, double);
double JupiterMagnitude(double, double, double);
double SaturnMagnitude(double, double, double, double);
double UranusMagnitude(double, double, double);
double NeptuneMagnitude(double, double, double);
double PlutoMagnitude(double, double, double);
double AsteroidMagnitude(double, double, double, double, double);
double CometMagnitude(double, double, double, double);

/* MarsMoon.c */

void SetPhobosMatrix(double[], double);
void SetDeimosMatrix(double[], double);
void PhobosOrbit(double, double *, double *, double *, double *, double *, double *);
void DeimosOrbit(double, double *, double *, double *, double *, double *, double *);
void PhobosRotation(double, double *, double *, double *);
void DeimosRotation(double, double *, double *, double *);

/* JupiMoon.c */

void SetJupiterMoonMatrix(double[]);
void IoXYZ(double, double *, double *, double *);
void EuropaXYZ(double, double *, double *, double *);
void GanymedeXYZ(double, double *, double *, double *);
void CallistoXYZ(double, double *, double *, double *);
void IoRotation(double, double *, double *, double *);
void EuropaRotation(double, double *, double *, double *);
void GanymedeRotation(double, double *, double *, double *);
void CallistoRotation(double, double *, double *, double *);

/* SatuMoon.c */

void SetSaturnMoonMatrix(double[]);
void MimasOrbit(double, double *, double *, double *, double *, double *, double *);
void EnceladusOrbit(double, double *, double *, double *, double *, double *, double *);
void TethysOrbit(double, double *, double *, double *, double *, double *, double *);
void DioneOrbit(double, double *, double *, double *, double *, double *, double *);
void RheaOrbit(double, double *, double *, double *, double *, double *, double *);
void TitanOrbit(double, double *, double *, double *, double *, double *, double *);
void HyperionOrbit(double, double *, double *, double *, double *, double *, double *);
void IapetusOrbit(double, double *, double *, double *, double *, double *, double *);
void PhoebeOrbit(double, double *, double *, double *, double *, double *, double *);
void MimasRotation(double, double *, double *, double *);
void EnceladusRotation(double, double *, double *, double *);
void TethysRotation(double, double *, double *, double *);
void DioneRotation(double, double *, double *, double *);
void RheaRotation(double, double *, double *, double *);
void TitanRotation(double, double *, double *, double *);
void IapetusRotation(double, double *, double *, double *);
void PhoebeRotation(double, double *, double *, double *);

/* UranMoon.c */

void SetUranusMoonMatrix(double[]);
void MirandaOrbit(double, double *, double *, double *, double *, double *, double *);
void ArielOrbit(double, double *, double *, double *, double *, double *, double *);
void UmbrielOrbit(double, double *, double *, double *, double *, double *, double *);
void TitaniaOrbit(double, double *, double *, double *, double *, double *, double *);
void OberonOrbit(double, double *, double *, double *, double *, double *, double *);
void MirandaRotation(double, double *, double *, double *);
void ArielRotation(double, double *, double *, double *);
void UmbrielRotation(double, double *, double *, double *);
void TitaniaRotation(double, double *, double *, double *);
void OberonRotation(double, double *, double *, double *);

/* NeptMoon.c */

void SetTritonMatrix(double[], double);
void SetNereidMatrix(double[]);
void TritonOrbit(double, double *, double *, double *, double *, double *, double *);
void NereidOrbit(double, double *, double *, double *, double *, double *, double *);
void CharonOrbit(double, double *, double *, double *, double *, double *, double *);
void TritonRotation(double, double *, double *, double *);
void CharonRotation(double, double *, double *, double *);

#endif /* !__ASTROLIB_H__ */
