#ifndef CALIB_H
#define CALIB_H

#include "sundials_types.h"

/* Functions Setting Global Calibration Paramters */
void fileNameID(const char *);
void setFileName(char *);

void initCalibFile(char *);

// For initialize.c
realtype setsatD_CALIB();
realtype setbr_CALIB();
realtype setporos_CALIB();
realtype seticsat_CALIB();
realtype setrivEle_CALIB();

//For et_is.c
realtype setis_CALIB();
realtype setet0_CALIB();
realtype setmf_CALIB();
realtype settf_CALIB();

//For f.c
realtype setVic_CALIB();
realtype setrivK_CALIB();
realtype setKh_CALIB();
realtype setRec_CALIB();
realtype setet2_CALIB();
realtype setet1_CALIB();
realtype setsat_THRESH();
realtype setmp_MULTFH();
realtype setmp_MULTFV();
realtype setmpArea_CALIB();
realtype setovl_THRESH_H();
realtype setovl_THRESH_V();
realtype setrzd_CALIB();

realtype setroughEle_CALIB();
realtype setroughRiv_CALIB();
realtype setrivCoeff_CALIB();
realtype setrivDepth_CALIB();
realtype setalpha_CALIB();
realtype setset_MP();
realtype setlai_CALIB();
realtype setvegfrac_CALIB();
realtype setalbedo_CALIB();

#endif
