


/*####################################################*/
/*                                                    */
/*      Includes for the prediction program           */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/*                                                    */
/*####################################################*/
// standard include //
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>   
#include <strings.h>

// standard multithread include //
//#include <pthread.h>

// Gnu Scientific library include //
//    www.gnu.org/software/gsl/   //
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>


/*####################################################*/
/*                                                    */
/*             constantes definitions                 */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/*                                                    */
/*####################################################*/
#define coef_spline 3
#define NB_WAVE_IN_SPECTRUM 36
#define RESOLUTION 0.125
#define SEMI_DIURNAL 1
#define DIURNAL 2
#define LONG 3

#define MASK -9999.000

/*####################################################*/
/*                                                    */
/*           variables struct definitions             */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/*                                                    */
/*####################################################*/

//wave struct contains all the tials waves parameters
typedef struct {
  float Ap;                   /*  specify tidal potentiel amplitude */
  int nT,ns,nh,np,np1;        /*  specify the main wave fequency V  */
  int nksi,nnu0,nnu1,nnu2,R;    /*  specify nodal argument u          */
  float shift;                /*  specify nodal argument u          */
  int formula,code;           /*  formula index for nodal factor f + code for JMM */
  double omega;               /*  specify pulsation, redondant with nT,ns,nh,np,np1  */
  char   name[10];            /*   */
  char   spec;
  int admit;
} tidal_wave;

//simple complex definition
typedef struct 
   {
     float  reel,imag;
   } fcomplex;

//the struct that contains all the predicted wave parameters
typedef struct 
{
  fcomplex *buffer;
  tidal_wave wave;
  fcomplex mask,Z;
  char PATH[256];
  int rstatus, headN ,founded,CTO;
  double lat,lon,time, prediction;
} spectrum_struct;



typedef struct 
{
  double amp,phi;
} wave_cst_struct;

/* typedef struct  */
/* { */
/*   wave_cst_struct M2,M4,K1,K2,O1,P1,Q1,S2,N2,DN2,MM,MF,MTM,MSQM; */
/*   wave_cst_struct NU2,MU2,L2,T2,LA2,KJ2,R2,OO1,J1,PHI1,PI1,PSI1,RO1,SIG1,TTA1,DQ1,KI1,MSM,MSF,MQM,MSTM,SSA; */
  
/* } fes2004_extract_struct; */




//the multithread struct
typedef struct
{
  int cnt,CTO,ncid;
  double *lat,*lon,*time,*Otide;
  spectrum_struct *spectrum;
  gsl_matrix *Mat;
  gsl_vector *sec_r,*sec_i,*res_r,*res_i;
  gsl_permutation *Perm;
  //fes2004_extract_struct *constante;
  double **amplitude,**phase; 
  //allocatable variable need in the thread
  float *data_amp;
  float *data_phi;
  double *weight;
  int *sindice,*aindice;
} mega_struct;

typedef struct 
{
  double sh_T,sh_h,sh_s,sh_p1,sh_p;
  double sh_xi,sh_nu,sh_R,sh_x1ra,sh_nuprim,sh_nusec;
  double sh_I,sh_N;
}astro_ang_struct;

typedef struct 
{
  int   day,month,year;   /* Gregorian date (define uniquely a given day)         */
  float second;           /* elapsed time in second from 0h00)*/
} date_t;




/*####################################################*/
/*                                                    */
/*               functions prototype                  */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/*                                                    */
/*####################################################*/
//fes2004_alloc.c
extern void alloc_tide_spectrum(spectrum_struct ** );
extern void alloc_prediction_threads(mega_struct *,int );
extern void alloc_extraction_threads(mega_struct *,int );
extern void free_threads(mega_struct *,int );

//fes2004_extraction.c
extern int fes2004_extraction (char *,int ,double *,double *,double **,double **,int );

//fes2004_prediction.c
extern int fes2004_prediction (char *,int ,int ,double *,double *,double *,double *,int );

//fes2004_io.c
extern void  load_netcdf_fes2004_data(char *,mega_struct *,int );

//fes2004_init.c
extern void in_out_file_open(char *,FILE **,FILE **);
extern int init_spectrum(spectrum_struct *, int );
extern void  init_thread_struct(int , mega_struct *, int ,double *,double *, double *, double *, spectrum_struct *,  double ** ,double **);
extern int Wave_select(int i, tidal_wave *,int );
extern void  init_admittance_coeff(int *,int *,int *,int );

//fes2004_error.c
extern void print_error_1() ;
extern void print_error_2() ;
extern void print_error_3(char *) ;
extern void print_error_4(char *) ;
extern void print_error_5(char *) ;

//fes2004_kernel.c
extern void multi_t_prediction(int , int , mega_struct *);
extern void multi_t_extraction(int , int , mega_struct *);
extern void * pred_coeur(void *);
extern void * extract_coeur(void *);
extern double predic_and_admit(float *,float *,int ,mega_struct *,double ,double );
extern int interpolation_w_mask(double ,double ,int ,int ,double *,float *);
extern void compute_admittance(gsl_matrix *,gsl_permutation *,gsl_vector *,gsl_vector *,gsl_vector *,gsl_vector *,spectrum_struct *,int ,int * ,int *,double);
extern double geo_d_km(double ,double ,double ,double );
//##############################################################################################################################
//          the next functions are a part of the Aktarus tide library developped 
//                          By Thierry LETELLIER [LEGOS]
//                          and Laurent ROBLOU [Noveltis]
//
//       These functions are the prediction kernel YOU DON'T NEED TO CHANGE THEM ...
//                                                ----------------------------------

//     You don't have the right to  use these functions in other programs than the FES2004 prediction
//     If you want more information please contact --- thierry.letellier@free.fr ---
//##############################################################################################################################
extern double Tide_prediction(double time,tidal_wave wave,fcomplex Z,int verbose,int CTO);
extern  void init_argument(double first_time, int verbose, int CTO,astro_ang_struct *astro_ang);
extern  void astronomic_angle(double tj, int verbose, astro_ang_struct *astro_ang);
extern  double greenwhich_argument(tidal_wave w,astro_ang_struct *astro_ang);
extern  double nodal_phase(tidal_wave w,astro_ang_struct *astro_ang);
extern int julian_day(int mm,int id,int iyyy);
extern double nodal_factort(int formula, astro_ang_struct *astro_ang);
extern void calendary(int njd,date_t *actual);
extern void getcnesdate(double t,date_t *actual);
extern double pulsation( tidal_wave wave);
extern void tidal_potential(tidal_wave wave,double lat,double lon, double *Amp, double *G,int terrestre);
//##############################################################################################################################





/*####################################################*/
/*                                                    */
/*               waves definitions                    */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/*                                                    */
/*####################################################*/
#ifdef LIB

/*----------------potentiel---T---s---h---p--p1-ksi--nu-nu1-nu2-R--shift-form-code-frq--name-spec-admit*/

 tidal_wave wmean = { 0.0000,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0,  0,   0,  187, 0.0, "mean",'p',999};/*LR: code=187*/

 tidal_wave wSa   = { 0.0000,  0,  0,  1,  0,  0,  0,  0,  0,  0, 0,   0,   0, 47, 0.0, "Sa",'p',999};
 tidal_wave wSsa  = { 1.9416,  0,  0,  2,  0,  0,  0,  0,  0,  0, 0,   0,   0, 46, 0.0, "Ssa",'p',999};
 tidal_wave wMm   = { 2.2056,  0, +1,  0, -1,  0,  0,  0,  0,  0, 0,   0,  73, 38, 0.0, "Mm",'p',999};
 tidal_wave wMSm  = { 0.3094,  0, +1, -2, +1,  0,  0,  0,  0,  0, 0,   0,  73, 41, 0.0, "MSm",'p',38};/*LR: 31/03/03 */
 tidal_wave wMSf  = { 0.2240,  0, +2, -2,  0,  0,  0,  0,  0,  0, 0,   0,  73, 39, 0.0, "MSf",'p',40};/*LR: 31/03/03  */
 tidal_wave wMf   = { 4.1765,  0, +2,  0,  0,  0, -2,  0,  0,  0, 0,   0,  74, 40, 0.0, "Mf",'p',999};
 tidal_wave wMtm  = { 0.8081,  0, +3,  0, -1,  0, -2,  0,  0,  0, 0,   0,  74, 42, 0.0, "Mtm",'p',999};/*LR: 24/01/02  */
 tidal_wave wMqm  = { 0.0000,  0, +4,  0, -2,  0, -2,  0,  0,  0, 0,   0,  74, 43, 0.0, "Mqm",'p',999};
 tidal_wave wMStm = { 0.1147,  0, +3, -2,  1,  0, -2,  0,  0,  0, 0,   0,  74, 44, 0.0, "MStm",'p',42};/*LR: 31/03/03  */
 tidal_wave wMSqm = { 0.0667,  0, +4, -2,  0,  0, -2,  0,  0,  0, 0,   0,  74, 45, 0.0, "MSqm",'p',999};/*LR: 24/01/02  */

/*----------------potentiel---T---s---h---p--p1-ksi--nu-nu1-nu2-R--shift-form-code-frq--name-spec-admit*/
 tidal_wave wO1   = {10.0573, +1, -2, +1,  0,  0, +2, -1,  0,  0, 0, +90,  75,  1, 0.0, "O1",'p',999};
 tidal_wave wQ1   = { 1.9469, +1, -3, +1, +1,  0, +2, -1,  0,  0, 0, +90,  75, 27, 0.0, "Q1",'p',999};
 tidal_wave w2Q1  = { 0.2587, +1, -4, +1, +2,  0, +2, -1,  0,  0, 0, +90,  75, 67, 0.0, "2Q1",'p',999};/*LR: 24/01/02  */
 tidal_wave wRo1  = { 0.3787, +1, -3, +3, -1,  0, +2, -1,  0,  0, 0, +90,  75, 69, 0.0, "Ro1",'p',999};/*LR: 24/01/02  */
 tidal_wave wSig1 = { 0.1627, +1, -4, +3,  0,  0, +2, -1,  0,  0, 0, +90,  75, 68, 0.0, "Sig1",'p',999};/*LR: 24/01/02  */
 tidal_wave wJ1   = { 0.7921, +1, +1, +1, -1,  0,  0, -1,  0,  0, 0, -90,  76, 29, 0.0, "J1",'p',999};/*avant s=-2*/
 tidal_wave wKi1  = { 0.1120, +1, -1, +3, -1,  0,  0, -1,  0,  0, 0, -90,  76, 70, 0.0, "Ki1",'p',999};/*LR: 24/01/02  */
 tidal_wave wTta1 = { 0.1120, +1, +1, -1, +1,  0,  0, -1,  0,  0, 0, -90,  76, 73, 0.0, "Tta1",'p',999};/*LR: 24/01/02  */
 tidal_wave wS1   = { 0.0000, +1,  0,  0,  0,  0,  0,  0,  0,  0, 0,   0,   0, 26, 0.0, "S1",'p',3};/* LR: 28/03/03 */
 tidal_wave wPi1  = { 0.2747, +1,  0, -2,  0, +1,  0,  0,  0,  0, 0, +90,   0, 71, 0.0, "Pi1",'p',999};
 tidal_wave wK1   = {14.1484, +1,  0, +1,  0,  0,  0,  0, -1,  0, 0, -90, 227,  3, 0.0, "K1",'p',999};
 tidal_wave wP1   = { 4.6806, +1,  0, -1,  0,  0,  0,  0,  0,  0, 0, +90,   0,  2, 0.0, "P1",'p',3};
 tidal_wave wPsi1 = { 0.1120, +1,  0, +2,  0, -1,  0,  0,  0,  0, 0, -90,   0,  0, 0.0, "Psi1",'p',3};/*LR: 31/03/03  */
 tidal_wave wPhi1 = { 0.2027, +1,  0, +3,  0,  0,  0,  0,  0,  0, 0, -90,   0, 72, 0.0, "Phi1",'p',999};/*LR: 24/01/02  */
 tidal_wave wOO1  = { 0.4347, +1, +2, +1,  0,  0, -2, -1,  0,  0, 0, -90,  77, 28, 0.0, "OO1",'p',999};
 tidal_wave wMP1  = { 0.0000, +1, -2, +3,  0,  0,  0, -1,  0,  0, 0, -90,  76, 49, 0.0, "MP1",'n',999};/*LR: 2/04/02*/
 tidal_wave wSO1  = { 0.0000, +1, +2, -1,  0,  0,  0, -1,  0,  0, 0, -90,  76, 48, 0.0, "SO1",'n',999};/*LR: 2/04/02*/
 tidal_wave wKQ1  = { 0.0000, +1, +3, +1, -1,  0, -2, -1,  0,  0, 0, -90,  15, 76, 0.0, "KQ1",'n',999};/*LR: 2/04/02*/
 tidal_wave wM1   = { 0.9788, +1, -1, +1,  0,  0, +1, -1,  0,  0, 0,   0, 144, 74, 0.0, "M1",'p',999};/*LR: 2/04/02*/

/*----------------potentiel---T---s---h---p--p1-ksi--nu-nu1-nu2-R--shift-form-code-frq--name-spec-admit*/
 tidal_wave wE2   = { 0.1789, +2, -5, +4, +1,  0, +2, -2,  0,  0, 0,   0,  78,  4, 0.0, "E2",'p',5};/*LR: 2/04/02*/
 tidal_wave w2N2  = { 0.6267, +2, -4, +2, +2,  0, +2, -2,  0,  0, 0,   0,  78,  5, 0.0, "2N2",'p',999};
 tidal_wave wMu2  = { 0.5841, +2, -4, +4,  0,  0, +2, -2,  0,  0, 0,   0,  78,  6, 0.0, "Mu2",'p',5};/*avant p=-1*/
 tidal_wave wN2   = { 4.6313, +2, -3, +2, +1,  0, +2, -2,  0,  0, 0,   0,  78,  7, 0.0, "N2",'p',999};
 tidal_wave wNu2  = { 0.9094, +2, -3, +4, -1,  0, +2, -2,  0,  0, 0,   0,  78,  8, 0.0, "Nu2",'p',7};
 tidal_wave wL2   = { 0.6694, +2, -1, +2, -1,  0, +2, -2,  0,  0, 1,+180, 215, 11, 0.0, "L2",'p',14};
 tidal_wave wLa2  = { 0.1760, +2, -1,  0, +1,  0, +2, -2,  0,  0, 0,+180,  78, 10, 0.0, "La2",'p',999};/*LR: 24/01/02 */
 tidal_wave wT2   = { 0.6614, +2,  0, -1,  0, +1,  0,  0,  0,  0, 0,   0,   0, 12, 0.0, "T2",'p',14};/*LR: 31/03/03  */
 tidal_wave wS2   = {11.2734, +2,  0,  0,  0,  0,  0,  0,  0,  0, 0,   0,   0, 13, 0.0, "S2",'p',999};
 tidal_wave wM2   = {24.2297, +2, -2, +2,  0,  0, +2, -2,  0,  0, 0,   0,  78,  9, 0.0, "M2",'p',999};
 tidal_wave wK2   = { 3.0697, +2,  0, +2,  0,  0,  0,  0,  0, -2, 0,   0, 235, 14, 0.0, "K2",'p',13};/*LR: 06/06/04  */
 tidal_wave wKJ2  = { 0.1707, +2, +1, +2, -1,  0,  0,  0,  0, -2, 0,   0,  79, 77, 0.0, "KJ2",'n',999};/*LR: 24/01/02 */
 tidal_wave wR2   = { 0.0933, +2,  0, +1,  0, -1,  0,  0,  0,  0, 0,   0,   0, 50, 0.0, "R2",'p',14};/*LR: 31/03/03  */
 tidal_wave wOQ2  = { 0.0000, +2, -5, +2, +1,  0,  0,  0,  0,  0, 0,+180,   7, 51, 0.0, "OQ2",'n',999};/*LR: 2/04/02*/
 tidal_wave w2MK2 = { 0.0000, +2, +4, +2,  0,  0,  4, -4,  0, +2, 0,   0,   5, 65, 0.0, "2MK2",'n',999};/*LR: 2/04/02*/
 tidal_wave wMSK2 = { 0.0000, +2, -2,  0,  0,  0, +2, -2,  0, +2, 0,   0,   8, 31, 0.0, "MSK2",'n',999};/*LR: 2/04/02*/
 tidal_wave wMSN2 = { 0.0000, +2, +1,  0, +1,  0, +2, -2,  0, +2, 0,   0,   6, 15, 0.0, "MSN2",'n',999};/*LR: 2/04/02*/
 tidal_wave w2SM2 = { 0.0000, +2, +2, -2,  0,  0, -2, +2,  0,  0, 0,   0,  16, 16, 0.0, "2SM2",'n',999};/*LR: 2/04/02*/
 tidal_wave wM_SK_2={ 0.0000, +2, -2, +1,  0,  0, +2, -2,  1,  0, 0, +90,   9, 37, 0.0, "M(SK)2",'n',999};/*LR: 2/04/02*/
 tidal_wave wM_KS_2={ 0.0000, +2, -2, +3,  0,  0, +2, -2, -1,  0, 0, -90,   9, 36, 0.0, "M(KS)2",'n',999};/*LR: 2/04/02*/
 tidal_wave wMKS2 = { 0.0000, +2, -2, +4,  0,  0, -2, -2,  0,  0, 0,   0,   8, 30, 0.0, "MKS2",'n',999};/*LR: 2/04/02*/
 tidal_wave wOP2  = { 0.0000, +2, -2,  0,  0,  0, +2, -1,  0,  0, 0,+180,  11,100, 0.0, "OP2",'n',999};/*LR: 2/04/02*/
 tidal_wave wMNS2 = { 0.0000, +2, -5, +4, +1,  0, +4, -4,  0,  0, 0,   0,   6,101, 0.0, "MNS2",'n',999};/*LR: 24/01/02 */

/*----------------potentiel---T---s---h---p--p1-ksi--nu-nu1-nu2-R--shift-form-code-frq--name-spec-admit*/
 tidal_wave wM3   = { 0.4747, +3, -3, +3,  0,  0, +3, -3,  0,  0, 0,   0,   1, 34, 0.0, "M3",'n',999};/*LR: 24/01/02 */
 tidal_wave wS3   = { 0.0000, +3,  0,  0,  0,  0,  0,  0,  0,  0, 0,   0,   1, 35, 0.0, "S3",'n',999};/*LR: 2/04/02*/
 tidal_wave w2MK3 = { 0.0000, +3, -4, +3,  0,  0, +4, -4, +1,  0, 0, +90,  10, 25, 0.0, "2MK3",'n',999};/*LR: 2/04/02*/
 tidal_wave wSO3  = { 0.0000, +3, -2, +1,  0,  0, +2, -1,  0,  0, 0, +90,  11, 53, 0.0, "SO3",'n',999};/*LR: 2/04/02*/
 tidal_wave wMK3  = { 0.0000, +3, -2, +3,  0,  0, +2, -2, -1,  0, 0, -90,  10, 24, 0.0, "MK3",'n',999};/*LR: 2/04/02*/
 tidal_wave wSK3  = { 0.0000, +3,  0, +1,  0,  0,  0,  0, -1,  0, 0, -90,  17, 54, 0.0, "SK3",'n',999};/*LR: 2/04/02*/
 tidal_wave wMO3  = { 0.0000, +3, -4, +1,  0,  0, +2, -2,  0,  0, 0, +90,  13,102, 0.0, "MO3",'n',999};/*LR: 2/04/02*/

/*----------------potentiel---T---s---h---p--p1-ksi--nu-nu1-nu2-R--shift-form-code-frq--name-spec-admit*/
 tidal_wave wN4   = { 0.0000, +4, -6, +4, +2,  0, +4, -4,  0,  0, 0,   0,   1, 33, 0.0, "N4",'n',999};
 tidal_wave wM4   = { 0.0000, +4, -4, +4,  0,  0, +4, -4,  0,  0, 0,   0,   1, 18, 0.0, "M4",'n',999};
 tidal_wave wS4   = { 0.0000, +4,  0,  0,  0,  0,  0,  0,  0,  0, 0,   0,   0, 56, 0.0, "S4",'n',999};
 tidal_wave wMN4  = { 0.0000, +4, -5, +4,  1,  0, +4, -4,  0,  0, 0,   0,   1, 17, 0.0, "MN4",'n',999};
 tidal_wave wMS4  = { 0.0000, +4, -2, +2,  0,  0, +2, -2,  0,  0, 0,   0,   2, 19, 0.0, "MS4",'n',999};
 tidal_wave wMK4  = { 0.0000, +4, -2, +4,  0,  0, +2, -2,  0, -2, 0,   0,   4, 20, 0.0, "MK4",'n',999};
 tidal_wave wSN4  = { 0.0000, +4, -3, +2, +1,  0, +2, -2,  0,  0, 0,   0,   2, 55, 0.0, "SN4",'n',999};/*LR: 29/03/02 */
 tidal_wave w3MS4 = { 0.0000, +4, -6, +6,  0,  0, +6, -6,  0,  0, 0,   0,  12, 58, 0.0, "3MS4",'n',999};/*LR: 29/03/02 */
 tidal_wave wSK4  = { 0.0000, +4,  0, +2,  0,  0,  0,  0,  0, -2, 0,   0,  14,103, 0.0, "SK4",'n',999};

/*----------------potentiel---T---s---h---p--p1-ksi--nu-nu1-nu2-R--shift-form-code-frq--name-spec-admit*/
 tidal_wave wM6   = { 0.0000, +6, -6, +6,  0,  0, +6, -6,  0,  0, 0,   0,  18, 22, 0.0, "M6",'n',999};
 tidal_wave w2MN6 = { 0.0000, +6, -7, +6,  1,  0, +6, -6,  0,  0, 0,   0,  18, 21, 0.0, "2MN6",'n',999};
 tidal_wave w2MS6 = { 0.0000, +6, -4, +4,  0,  0, +4, -4,  0,  0, 0,   0,   6, 59, 0.0, "2MS6",'n',999};
 tidal_wave w2MK6 = { 0.0000, +6, -4, +6,  0,  0, +4, -4,  0, -2, 0,   0,   5, 60, 0.0, "2MK6",'n',999};
 tidal_wave wMSN6 = { 0.0000, +6, -5, +4, +1,  0, +4, -4,  0,  0, 0,   0,   6, 23, 0.0, "MSN6",'n',999};/*LR: 29/03/02 */
 tidal_wave w2SM6 = { 0.0000, +6, -2, +2,  0,  0, +2, -2,  0,  0, 0,   0,  16,104, 0.0, "2SM6",'n',999};/*LR: 2/04/02*/
 tidal_wave wMSK6 = { 0.0000, +6, -2, +4,  0,  0, +2, -2,  0, -2, 0,   0,   8,105, 0.0, "MSK6",'n',999};/*LR: 2/04/02*/

/*----------------potentiel---T---s---h---p--p1-ksi--nu-nu1-nu2-R--shift-form-code-frq--name-spec-admit*/
 tidal_wave w3MS8 = { 0.0000, +8, -6, +6,  0,  0, +6, -6,  0,  0, 0,   0,  19, 61, 0.0, "3MS8",'n',999};/*LR: 2/04/02*/

/*----------------potentiel---T---s---h---p--p1-ksi--nu-nu1-nu2-R--shift-form-code-frq--name-spec-admit*/
 tidal_wave wRH5  = { 0.0000,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0,   0,   0, 78, 3.0, "RH5",'a',999};
  
  
/*----------------potentiel---T---s---h---p--p1--ksi-nu-nu1-nu2-shift-form-code-frq--name-linear*/
 tidal_wave wNUL   = {0.0000, 0,  0,  0,  0,  0,  0,  0, 0,  0,  0,0,   0,  0,  0.0, "NUL",'p',999};  /*TL: 08/06/2004*/
 
  double  pi = 3.14159265358979323844;    
  double  pi2 = 3.14159265358979323844*2.;    
  double  deg_to_rad = 3.14159265358979323844/180.;   

#else
  extern tidal_wave wmean; 

  extern tidal_wave wSa   ;
  extern tidal_wave wSsa  ;
  extern tidal_wave wMm   ;
  extern tidal_wave wMSm  ;
  extern tidal_wave wMSf  ;
  extern tidal_wave wMf   ;
  extern tidal_wave wMtm  ;
  extern tidal_wave wMqm  ;
  extern tidal_wave wMStm ;
  extern tidal_wave wMSqm ;

  extern tidal_wave wO1   ; 
  extern tidal_wave wQ1   ;
  extern tidal_wave w2Q1  ;
  extern tidal_wave wRo1  ;
  extern tidal_wave wSig1 ;
  extern tidal_wave wJ1   ;
  extern tidal_wave wKi1  ;
  extern tidal_wave wTta1 ;
  extern tidal_wave wS1   ;
  extern tidal_wave wPi1  ;
  extern tidal_wave wK1   ;
  extern tidal_wave wP1   ;
  extern tidal_wave wPsi1 ;
  extern tidal_wave wPhi1 ;
  extern tidal_wave wOO1  ;
  extern tidal_wave wMP1  ;/*LR: 2/04/02*/
  extern tidal_wave wSO1  ;/*LR: 2/04/02*/
  extern tidal_wave wKQ1  ;/*LR: 2/04/02*/
  extern tidal_wave wM1   ;/*LR: 2/04/02*/

  extern tidal_wave wE2   ;
  extern tidal_wave w2N2  ;
  extern tidal_wave wMu2  ;
  extern tidal_wave wN2   ;
  extern tidal_wave wNu2  ;
  extern tidal_wave wM2   ;
  extern tidal_wave wL2   ;
  extern tidal_wave wLa2  ;
  extern tidal_wave wT2   ;
  extern tidal_wave wS2   ;
  extern tidal_wave wK2   ;
  extern tidal_wave wKJ2  ;
  extern tidal_wave wR2   ;
  extern tidal_wave wOQ2  ;/*LR: 2/04/02*/
  extern tidal_wave w2MK2 ;/*LR: 2/04/02*/
  extern tidal_wave wMSK2 ;/*LR: 2/04/02*/
  extern tidal_wave wMSN2 ;/*LR: 2/04/02*/
  extern tidal_wave w2SM2 ;/*LR: 2/04/02*/
  extern tidal_wave wM_SK_2 ;/*LR: 2/04/02*/
  extern tidal_wave wM_KS_2;/*LR: 2/04/02*/
  extern tidal_wave wMKS2 ;/*LR: 2/04/02*/
  extern tidal_wave wOP2  ;/*LR: 2/04/02*/
  extern tidal_wave wMNS2 ;/*LR: 2/04/02*/

  extern tidal_wave wM3   ;
  extern tidal_wave wS3   ;
  extern tidal_wave w2MK3 ;/*LR: 2/04/02*/
  extern tidal_wave wSO3  ;/*LR: 2/04/02*/
  extern tidal_wave wMK3  ;/*LR: 2/04/02*/
  extern tidal_wave wSK3  ;/*LR: 2/04/02*/
  extern tidal_wave wMO3  ;/*LR: 2/04/02*/

  extern tidal_wave wN4   ;
  extern tidal_wave wM4   ;
  extern tidal_wave wS4   ;
  extern tidal_wave wMN4  ;
  extern tidal_wave wMS4  ;
  extern tidal_wave wMK4  ;
  extern tidal_wave wSN4  ;/*LR: 2/04/02*/
  extern tidal_wave w3MS4 ;/*LR: 2/04/02*/
  extern tidal_wave wSK4  ;/*LR: 2/04/02*/
 
  extern tidal_wave wM6   ;
  extern tidal_wave w2MN6 ;
  extern tidal_wave w2MS6 ;
  extern tidal_wave w2MK6 ;
  extern tidal_wave wMSN6 ;/*LR: 29/03/02 */
  extern tidal_wave w2SM6 ;/*LR: 2/04/02*/
  extern tidal_wave wMSK6 ;/*LR: 2/04/02*/

  extern tidal_wave w3MS8 ;/*LR: 2/04/02*/

  extern tidal_wave wRH5  ;
  extern tidal_wave wNUL  ;
  
  extern double pi;
  extern double pi2;
  extern double  deg_to_rad;      

  
#endif
