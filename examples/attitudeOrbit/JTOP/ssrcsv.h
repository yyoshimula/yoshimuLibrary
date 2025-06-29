
/** ssrcsv.h ******************************************************************
 *                                                                            *
 *  coding:Feb/08/96/1.00:1st coding by T.HANADA                              *
 *  update:May/03/03/1.01: The length of strings modified.                    *
 *  update:___/__/__/_.__:                                                    *
 *  bugfix:___/__/__/_.__:                                                    *
 *                                                                            *
 ******************************************************************************/

#ifndef __SSRCSV_H__
#define __SSRCSV_H__
//345678901234567890123456789012345678901234567890123456789012345678901234567890

typedef struct
{
  char OBJECT_ID[12];         /* OBJECT_ID ---------------------------------- */
  char OBJECT_NAME[25];       /* OBJECT_NAME -------------------------------- */
  int NORAD_CAT_ID;           /* NORAD_CAT_ID ------------------------------- */
  char COUNTRY[6];            /* COUNTRY ------------------------------------ */
  double PERIOD;              /* PERIOD, in minute -------------------------- */
  double INCLINATION;         /* INCLINATION, in degrees -------------------- */
  double APOGEE;              /* APOGEE, in kilometers ---------------------- */
  double PERIGEE;             /* PERIGEE, om kilometers --------------------- */
  char RCS_SIZE[6];           /* RCS_SIZE ----------------------------------- */
  double RCSVALUE;            /* RCSVALUE, in meters squared ---------------- */
  char LAUNCH[12];            /* LAUNCH, in YYYY-MM-DD ---------------------- */
  char DECAY[12];             /* DECAY, in YYYY-MM-DD ----------------------- */
  char COMMENT[32];           /* COMMENT ------------------------------------ */
}
ssrcsv;

#endif /* !__SSRCSV_H__ */
