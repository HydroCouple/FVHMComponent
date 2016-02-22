/*******************************************************************************
 * File        : et_is.c                                                       *
 * Function    : for calculation of evapotranspiration and interception        *
 * Programmer  : Yizhong Qu @ Pennsylvania State Univeristy                    *
 * Version     : May, 2004 (1.0)                                               *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * Interception process and evapotranspiration are considered as weakly coupled*
 * processes compared with overland, channel routing and groundwater flow pro- *
 * cesses. Therefore, before each time step, interception is calculated and    *
 * deducted. After each time step, envapotranspiration is calculated and all   *
 * state variables are adjusted accordingly. In this file, ET is calculated    *
 * using Penman-Monteith equation. Most parameters are calculated based on     *
 * Crop Evapotranspiration (FAO No.56).                                        *
 *                                                                             *
 * This code is free for users with research purpose only, if appropriate      *
 * citation is refered. However, there is no warranty in any format for this   *
 * product.                                                                    *
 *                                                                             *
 * For questions or comments, please contact the authors of the reference.     *
 * One who want to use it for other consideration may also contact Dr.Duffy    *
 * at cxd11@psu.edu.                                                           *
 *******************************************************************************/

 /* TODO */
 //eliminate 'ret' variable by using EleTF

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "nvector_serial.h"
#include "sundials_types.h"
#include "ihm.h"
#include "calib.h"
#include "et_is.h"

#define EPSILON 0.05

realtype is_CALIB;
realtype et0_CALIB;
realtype mf_CALIB;
realtype tf_CALIB;

/*
#define is_CALIB	1.0
#define et0_CALIB	1.0
#define mf_CALIB	1.0
#define tf_CALIB	1.0
*/
realtype Interpolation(TSD *Data, realtype t);

void calIS(realtype t, realtype stepsize, void *DS,N_Vector VY)
{
	int i;
  	realtype totEvap;
  	realtype Delta, Gamma;
  	realtype Rn, G, T, Vel, RH, VP,P,LAI,zero_dh,cnpy_h,rl,r_a;
  	realtype isval=0,etval=0;
  	realtype fracSnow,snowRate,MeltRate,MF,Ts=-3.0,Tr=1.0,To=0.0,massMelt=0,ret;

  	is_CALIB = setis_CALIB();
  	et0_CALIB = setet0_CALIB();
  	mf_CALIB=setmf_CALIB();
  	tf_CALIB = settf_CALIB();
	//realtype *et_Y;

  	Model_Data MD;

  	MD = (Model_Data)DS;
	//et_Y = NV_DATA_S(VY);

 	stepsize=stepsize/(24.0*60.0);
  	for(i=0; i<MD->NumEle; i++)
  	{
    	MD->ElePrep[i] = Interpolation(&MD->TSD_Prep[MD->Ele[i].prep-1], t);
    	Rn = Interpolation(&MD->TSD_Rn[MD->Ele[i].Rn-1], t);
		//G = Interpolation(&MD->TSD_G[MD->Ele[i].G-1], t);
    	T = Interpolation(&MD->TSD_Temp[MD->Ele[i].temp-1], t);
    	Vel = Interpolation(&MD->TSD_WindVel[MD->Ele[i].WindVel-1], t);
    	RH = Interpolation(&MD->TSD_Humidity[MD->Ele[i].humidity-1], t);
    	VP = Interpolation(&MD->TSD_Pressure[MD->Ele[i].pressure-1], t);
    	P = 101.325*pow(10,3)*pow((293-0.0065*MD->Ele[i].zmax)/293,5.26);
    	LAI = Interpolation(&MD->TSD_LAI[MD->Ele[i].LC-1], t);
    	MF = Interpolation(&MD->TSD_MeltF[0], t);
    	MF=mf_CALIB*MF;

		/*****************************************Snow Calculation ****************************************************/
    	fracSnow = T<Ts?1.0:T>Tr?0:(Tr-T)/(Tr-Ts);
    	snowRate = fracSnow*MD->ElePrep[i];
    	MD->EleSnow[i]=MD->EleSnow[i]+snowRate*stepsize;
    	MeltRate=(T>To?(T-To)*MF:0);
		//MeltRate=(T>To?Rn*(T-To)*MF:0);
		//printf("\n %f %f %f %f",MF,fracSnow,MD->EleSnow[i],Rn*(T-To)*MF*stepsize);
    	if(MD->EleSnow[i]>MeltRate*stepsize)
    	{
    		MD->EleSnow[i]=MD->EleSnow[i]-MeltRate*stepsize;
    		//MeltRate=(T>To?Rn*(T-To)*MF/60.0:0);
    	}
    	else
    	{
    		MeltRate=MD->EleSnow[i]/stepsize;
    		MD->EleSnow[i]=0;
    	}
		/**************************************************************************************************************/



		/**************************************Evaporation from canopy*************************************************/
    	MD->EleISmax[i] = MD->SIFactor[MD->Ele[i].LC-1]*Interpolation(&MD->TSD_LAI[MD->Ele[i].LC-1], t);
		//MD->EleIS[i] = is_CALIB*MD->EleISmax[i];
		/// Bhatt
		MD->EleISmax[i] = is_CALIB*MD->EleISmax[i];
    	if(LAI>0.0)
    	{
	    	Delta = 2503*pow(10,3)*exp(17.27*T/(T+237.3))/(pow(237.3 + T, 2));
	    	Gamma = P*1.0035*0.92/(0.622*2441);
	    	zero_dh=Interpolation(&MD->TSD_DH[MD->Ele[i].LC-1], t);
			//zero_dh=0;
	    	cnpy_h = zero_dh/(1.1*(0.0000001+log(1+pow(0.007*LAI,0.25))));
			/*if(LAI<2.85)
	    	{
	    		rl= 0.0002 + 0.3*cnpy_h*pow(0.07*LAI,0.5);
	    	}
	    	else
	    	{
	    		rl= 0.3*cnpy_h*(1-(zero_dh/cnpy_h));
	    	}
			*/
			rl=Interpolation(&MD->TSD_DH[MD->Ele[i].LC-1], t);
	    	r_a = log(MD->Ele[i].windH/rl)*log(10*MD->Ele[i].windH/rl)/(Vel*0.16);

			/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ DELETE 0.01 BELOW $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

	    	MD->EleET[i][0] = (LAI/MD->Ele[i].LAImax)*(pow((MD->EleIS[i]<0?0:MD->EleIS[i])/MD->EleISmax[i],2.0/3.0))*(Rn*(1-MD->Ele[i].Albedo)*Delta+(1.2*1003.5*((VP/RH)-VP)/r_a))/(1000*2441000.0*(Delta+Gamma));
	    	MD->EleTF[i]=tf_CALIB*5.65*pow(10,-2)*MD->EleISmax[i]*exp(3.89*(MD->EleIS[i]<0?0:MD->EleIS[i])/MD->EleISmax[i]);
			//printf("\n %f %f %f %f %f %f %f %f %f",MD->EleIS[i],LAI,MD->Ele[i].LAImax,r_a,rl,cnpy_h,Delta,Gamma,MD->EleET[i][0]);

		}
     	else
     	{
     		MD->EleET[i][0]=0.0;
			MD->EleTF[i]=0.0;
     	}

		/**********************************************************************************************************************/
		/*if((t>144000)&&(t<322560))
		{
			MD->EleET[i][0]=0.1*et0_CALIB*MD->EleET[i][0];
		}
   		else
		{
			MD->EleET[i][0]=1.8*et0_CALIB*MD->EleET[i][0];
		}
		*/
		MD->EleET[i][0]=et0_CALIB*MD->EleET[i][0];
		//printf("\n%f %f",MD->EleIS[i],MD->EleISmax[i]);
		//getchar();
    	if(MD->EleIS[i] >= MD->EleISmax[i])
    	{
      		if((1-fracSnow)*MD->ElePrep[i]>=MD->EleET[i][0]+MD->EleTF[i])
      		{
      			MD->Ele2IS[i] = MD->EleET[i][0];
				ret = MD->EleTF[i];
      		}
     		else if(((1-fracSnow)*MD->ElePrep[i]<MD->EleET[i][0]+MD->EleTF[i])&&(MD->EleIS[i]+stepsize*((1-fracSnow)*MD->ElePrep[i]-MD->EleET[i][0]-MD->EleTF[i])<=0))
     		{
     			MD->Ele2IS[i] =(1-fracSnow)*MD->ElePrep[i];
     			MD->EleET[i][0]=(MD->EleET[i][0]/(MD->EleET[i][0]+MD->EleTF[i]))*(MD->EleIS[i]/stepsize+(1-fracSnow)*MD->ElePrep[i]);
				MD->EleIS[i] = 0;
				ret =(MD->EleTF[i]/(MD->EleET[i][0]+MD->EleTF[i]))*(MD->EleIS[i]/stepsize+(1-fracSnow)*MD->ElePrep[i]);
     		}
     		else
     		{
     			MD->Ele2IS[i] =(1-fracSnow)*MD->ElePrep[i];
     			MD->EleIS[i] = MD->EleIS[i]+stepsize*((1-fracSnow)*MD->ElePrep[i]-MD->EleET[i][0]-MD->EleTF[i]);
				ret =  MD->EleTF[i];
  			}



			/*MD->Ele2IS[i] = MD->ElePrep[i]>=MD->EleET[i][0]?MD->EleET[i][0]:MD->ElePrep[i];
      		MD->EleIS[i] = MD->ElePrep[i]>=MD->EleET[i][0]?MD->EleIS[i]:(MD->EleIS[i]+MD->ElePrep[i]-MD->EleET[i][0]<0?0:MD->EleIS[i]+MD->ElePrep[i]-MD->EleET[i][0]);
      		MD->EleET[i][0]=MD->EleIS[i]/stepsize+MD->ElePrep[i];
			*/
			//printf("\n here haha");
    	}
    	else if(((MD->EleIS[i] + ((1-fracSnow)*MD->ElePrep[i]-MD->EleET[i][0]-MD->EleTF[i])*stepsize) >= MD->EleISmax[i]))
    	{
      		MD->Ele2IS[i] =  (MD->EleISmax[i]/stepsize - (MD->EleIS[i]/stepsize-MD->EleET[i][0]-MD->EleTF[i]));
      		MD->EleIS[i] = MD->EleISmax[i];
      		ret =  MD->EleTF[i];
			//printf("\n here haha1");
    	}
    	else if(((MD->EleIS[i] + ((1-fracSnow)*MD->ElePrep[i]-MD->EleET[i][0]-MD->EleTF[i])*stepsize) <=0))
    	{
      		MD->Ele2IS[i] =  (1-fracSnow)*MD->ElePrep[i];
      		if((MD->EleET[i][0]>0)||(MD->EleTF[i]>0))
			{
        		MD->EleET[i][0] = (MD->EleET[i][0]/(MD->EleET[i][0]+MD->EleTF[i]))*(MD->EleIS[i]/stepsize+(1-fracSnow)*MD->ElePrep[i]);
        		ret =(MD->EleTF[i]/(MD->EleET[i][0]+MD->EleTF[i]))*(MD->EleIS[i]/stepsize+(1-fracSnow)*MD->ElePrep[i]);
			}
      		else
			{
				MD->EleET[i][0] = 0;
				ret = 0;
			}
      		MD->EleIS[i] = 0;
			//printf("\n here haha1");
    	}
    	else
    	{
      		MD->Ele2IS[i] = (1-fracSnow)*MD->ElePrep[i];
      		MD->EleIS[i] = MD->EleIS[i] + ((1-fracSnow)*MD->ElePrep[i]-MD->EleET[i][0]-MD->EleTF[i])*stepsize;
       		ret =  MD->EleTF[i];
			//printf("\n here haha2");
    	}
    	massMelt=massMelt+MeltRate*stepsize;
    	MD->EleNetPrep[i] = (1-MD->Ele[i].VegFrac)*(1-fracSnow)*MD->ElePrep[i] + ((1-fracSnow)*MD->ElePrep[i]+ret - MD->Ele2IS[i])*MD->Ele[i].VegFrac+MeltRate;
		MD->EleTF[i] = ret;
		/*if(MD->Ele2IS[i]>0)
		{
			printf("\ntime %lf %lf",t,MD->Ele2IS[i]);
			exit(1);
		}
		*/


		//MD->EleNetPrep[i] =MD->ElePrep[i];


		//printf("\n%d %f %f %f %f %f %f %e %e %e",i,MD->EleIS[i],MD->ElePrep[i],MD->Ele2IS[i],et_Y[i + MD->NumEle],LAI,zero_dh,MD->EleET[i][0],MD->EleET[i][1],(Rn*(1-MD->Ele[i].Albedo)*Delta+(1.2*1003.5*((VP/RH)-VP)/r_a))/(1000*2441000.0*(Delta+Gamma)));
		//printf("\n%d %f %f %f %f %f",i,MD->ElePrep[i],ret,MD->EleTF[i],MD->Ele2IS[i],MD->EleNetPrep[i]);
		//gtchar();
  	}
	//printf("\t%e",massMelt);
	//getchar();
}

void calET(realtype t, realtype stepsize, N_Vector VY, void *DS)
{
  int i;
/*  realtype Delta, Gamma, Es, Ea;
  realtype Rn, G, T, Vel, H, P;
  realtype ET_Value, ET_Remain;
  realtype *et_Y;
  */
 }
