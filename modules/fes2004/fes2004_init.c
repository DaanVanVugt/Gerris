#define LIB
#include "fes2004_lib.h"


/*####################################################*/
/*                                                    */
/*    received : 1 name                               */
/*    return : 2 pointers (type FILE)                 */
/*                                                    */
/*    the 2 files are open in r (in) and w (out)      */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/*                                                    */
/*####################################################*/

void in_out_file_open(char *name,FILE **in,FILE **out)
{
  char *outfilename;
  int len;

  *in=fopen(name,"r");
  len=strlen(name);
  outfilename=(char *)malloc( (len+15)*sizeof(char));
  sprintf(outfilename,"%s.output_file",name);    
  *out=fopen(outfilename,"w");

}

/*####################################################*/
/*                                                    */
/*    received and return: 1 spectrum                 */
/*                                                    */
/*    all the information in the spectrum are         */
/*    readable in the prediction.h                    */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/*                                                    */
/*####################################################*/

int init_spectrum(spectrum_struct *spectrum, int CTO)
{
  int rstatus;
  int i;

  //the rstatus return code is not test YET 

  // spectrum waves initialisation
  /*------------FES2004-------------------*/
  rstatus=Wave_select(8,&spectrum[0].wave,2); /*2N2*/
  rstatus=Wave_select(5,&spectrum[1].wave,2); /*K1*/
  rstatus=Wave_select(4,&spectrum[2].wave,2); /*K2*/
  rstatus=Wave_select(1,&spectrum[3].wave,2); /*M2*/
  rstatus=Wave_select(37,&spectrum[4].wave,2); /*M4*/
  rstatus=Wave_select(27,&spectrum[5].wave,2);/*Mf*/
  rstatus=Wave_select(28,&spectrum[6].wave,2);/*Mm*/
  rstatus=Wave_select(30,&spectrum[7].wave,2);/*Msqm*/
  rstatus=Wave_select(29,&spectrum[8].wave,2);/*Mtm*/
  rstatus=Wave_select(3,&spectrum[9].wave,2); /*N2*/
  rstatus=Wave_select(6,&spectrum[10].wave,2); /*O1*/
  rstatus=Wave_select(16,&spectrum[11].wave,2); /*P1*/
  rstatus=Wave_select(7,&spectrum[12].wave,2); /*Q1*/
  rstatus=Wave_select(2,&spectrum[13].wave,2); /*S2*/
  /*------------ADMITTED------------------*/
  /*semidiurne*/
  rstatus=Wave_select(10,&spectrum[14].wave,2);/*nu2*/
  rstatus=Wave_select(9, &spectrum[15].wave,2);/*mu2*/
  rstatus=Wave_select(11,&spectrum[16].wave,2);/*L2*/
  rstatus=Wave_select(12,&spectrum[17].wave,2);/*T2*/
  rstatus=Wave_select(13,&spectrum[18].wave,2);/*LA1*/
  rstatus=Wave_select(14,&spectrum[19].wave,2);/*KJ2*/
  rstatus=Wave_select(15,&spectrum[20].wave,2);/*R2*/
  /*diurne*/
  rstatus=Wave_select(17,&spectrum[21].wave,2);/*OO1*/
  rstatus=Wave_select(18,&spectrum[22].wave,2);/*J1*/
  rstatus=Wave_select(19,&spectrum[23].wave,2);/*PHI1*/
  rstatus=Wave_select(20,&spectrum[24].wave,2);/*PI1*/
  rstatus=Wave_select(21,&spectrum[25].wave,2);/*PSI1*/
  rstatus=Wave_select(22,&spectrum[26].wave,2);/*RO1*/
  rstatus=Wave_select(23,&spectrum[27].wave,2);/*SIG1*/
  rstatus=Wave_select(24,&spectrum[28].wave,2);/*TTA1*/
  rstatus=Wave_select(25,&spectrum[29].wave,2);/*2Q1*/
  rstatus=Wave_select(26,&spectrum[30].wave,2);/*Ki1*/
  /*long*/
  rstatus=Wave_select(33,&spectrum[31].wave,2);/*Msm*/
  rstatus=Wave_select(34,&spectrum[32].wave,2);/*Msf*/
  rstatus=Wave_select(35,&spectrum[33].wave,2);/*Mqm*/
  rstatus=Wave_select(36,&spectrum[34].wave,2);/*Mstm*/
  /*------------computed------------------*/
  rstatus=Wave_select(31,&spectrum[35].wave,2);/*SSA*/

  //spectrum data initialisation
  for(i=0;i<13;i++) 
    {
      spectrum[i].buffer=NULL;
      spectrum[i].rstatus=-1;
      sprintf(spectrum[i].PATH,"../data/%s.nc",spectrum[i].wave.name);
    } 
  for(i=0;i<36;i++) spectrum[i].CTO=CTO;
  return(0);
}


/*####################################################*/
/*                                                    */
/*                                                    */
/*   init the mega struct P to prepare the            */
/*   multithreaded part of the prediction             */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/*                                                    */
/*####################################################*/

void  init_thread_struct(int CPU, mega_struct *P, int CTO,double *lat,double *lon, double *time, double *Otide, spectrum_struct *spectrum,  double **amplitude,double **phase)
{
  int i,j;
  for(i=0;i<CPU;i++) 
    {
      P[i].CTO=CTO;
      P[i].lat=lat;
      P[i].lon=lon;
      P[i].time=time;
      P[i].Otide=Otide;
      P[i].amplitude=amplitude;
      P[i].phase=phase;
      for(j=0;j<36;j++)P[i].spectrum[j]=spectrum[j];
    }
}




/*####################################################*/
/*                                                    */
/*                                                    */
/*   init select the wave in prediction.h             */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/* this routine is a part of the aktarus lib LEGOS    */
/*####################################################*/
int Wave_select(int i, tidal_wave *carac,int detided)
{

  int verbose=1;

switch (i)
{ 

 case 1 : {*carac=wM2; break;}
 case 2 : {*carac=wS2; break;}
 case 3 : {*carac=wN2; break;}
 case 4 : {*carac=wK2; break;}
 case 5 : {*carac=wK1; break;}
 case 6 : {*carac=wO1; break;}
 case 7 : {*carac=wQ1; break;}
 case 8 : {*carac=w2N2; break;}
/*-----------------------------*/
 case 9 : {*carac=wMu2; break;}
 case 10 : { *carac=wNu2;break;}
 case 11 : {*carac=wL2; break;}
 case 12 : {*carac=wT2; break;}
 case 13 : {*carac=wLa2; break;} 
 case 14 :
   {
     if (detided==2)  {*carac=wKJ2; break;} 
     else {verbose=0;break;}
   }
 case 15 : 
   {
     if (detided==2)  {*carac=wR2; break;}
     else {verbose=0;break;}
   }
/*-----------------------------*/
 case 16 : {*carac=wP1; break;}
 case 17 : {*carac=wOO1; break;}
 case 18 : {*carac=wJ1; break;}
 case 19 : {*carac=wPhi1; break;}
 case 20 : {*carac=wPi1; break;}
 case 21 : 
  {
    if (detided==0)  {verbose=0;break; }
    else {*carac=wPsi1; break;}
  }
 case 22 : {*carac=wRo1; break;}
 case 23 : {*carac=wSig1; break;}
 case 24 : {*carac=wTta1; break;}
 case 25 : {*carac=w2Q1; break;}
 case 26 : {*carac=wKi1; break;}
/*-----------------------------*/
 case 27 : {*carac=wMf; break;}
 case 28 : {*carac=wMm; break;}
 case 29 : {*carac=wMtm; break;}
 case 30 : 
   {
     if (detided==2) {*carac=wMSqm; break;}
     else {verbose=0;break;}
   }
 case 31 : {*carac=wSsa; break;}
 case 32 : {*carac=wSa; break;}
 case 33 : 
   {
     if (detided==2)  {*carac=wMSm; break;}
     else {verbose=0;break;}
   }
 case 34 : {*carac=wMSf; break;}
 case 35 : {*carac=wMqm; break;}
 case 36 : {*carac=wMStm; break;}

/*-----------------------------*/
 case 37 : {*carac=wM4; break;}
 case 38 : {*carac=wMS4; break;}
 case 39 : {*carac=wMN4; break;}
 case 40 : {*carac=wS4; break;}
 case 41 : {*carac=wN4; break;}

/*-----------------------------*/
 case 42 : {*carac=wS1; break;}

 default : verbose=0;}/*switch*/


return(verbose);


}/*end*/ 



/*####################################################*/
/*                                                    */
/*       init the coef for the splien admittance      */
/*                                                    */
/*           Thierry LETELLIER 04Oct2005              */
/*                                                    */
/*####################################################*/

void  init_admittance_coeff(int *sindice,int *aindice,int *nb,int GROUPE)
{
  switch(GROUPE)
    {
    case SEMI_DIURNAL :
      { 
	sindice[0]=2;
	sindice[1]=3;
	sindice[2]=9;
	aindice[0]=14;
	aindice[1]=15;
	aindice[2]=16;
	aindice[3]=17;
	aindice[4]=18;
	aindice[5]=19;
	aindice[6]=20;
	*nb=7;
	break;
      }
    case DIURNAL :
      {
	sindice[0]=1;
	sindice[1]=10;
	sindice[2]=12;
	aindice[0]=21;
	aindice[1]=22;
	aindice[2]=23;
	aindice[3]=24;
	aindice[4]=25;
	aindice[5]=26;
	aindice[6]=27;
	aindice[7]=28;
	aindice[8]=29;
	aindice[9]=30;
	*nb=10;
	break;
      }
    case LONG :
      {
	sindice[0]=5;
	sindice[1]=6;
	sindice[2]=8;
	aindice[0]=31;
	aindice[1]=32;// LR, modif: 2008/12/18
	aindice[2]=33;// LR, modif: 2008/12/18
	aindice[3]=34;// LR, modif: 2008/12/18
	*nb=4;
	break;
      }
    }

}
