#include <netcdf.h>
#include "fes2004_lib.h"

//#define H_THREAD
//
// Uncomment this define to compile with the multithread option


/*####################################################*/
/*                                                    */
/*                                                    */
/*      running the multithread prediction            */
/*                                                    */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/*                                                    */
/*####################################################*/
void multi_t_prediction(int nocycle, int CPU, mega_struct *P)
{

#ifdef H_THREAD 
  pthread_t *tid;
#endif
		
  int cnt,i;

#ifdef H_THREAD 
  tid=malloc(CPU*sizeof(pthread_t));
  printf("starting a multithread prediction on %d processors \n",CPU);
#endif
  cnt=0;
  while(cnt<nocycle)
    {
      for(i=0;i<CPU;i++) 
	{
	  P[i].cnt=cnt;
#ifdef H_THREAD 
	  pthread_create(&tid[i], NULL,pred_coeur,(void *)&(P[i])); 
#else
	  pred_coeur( (void *)&(P[i]) );
#endif
	  cnt++;
	} 
#ifdef H_THREAD 
      for (i=0;i<CPU;i++)pthread_join(tid[i], NULL);  
#endif
      
      if(CPU>(nocycle-cnt))CPU=nocycle-cnt;
    }/* end loop on AT points*/
}


/*####################################################*/
/*                                                    */
/*                                                    */
/*      running the multithread extraction            */
/*                                                    */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/*                                                    */
/*####################################################*/
void multi_t_extraction(int nocycle, int CPU, mega_struct *P)
{

#ifdef H_THREAD 
  pthread_t *tid;
#endif
  int cnt,i;

#ifdef H_THREAD 
  tid=malloc(CPU*sizeof(pthread_t));
  printf("starting a multithread extraction on %d processors \n",CPU);
#endif
  cnt=0;
  while(cnt<nocycle)
    {
      for(i=0;i<CPU;i++) 
	{
	  P[i].cnt=cnt;
#ifdef H_THREAD 
	  pthread_create(&tid[i], NULL,extract_coeur,(void *)&(P[i])); 
#else
	  extract_coeur( (void *)&(P[i]) );
#endif
	  cnt++;
	} 
#ifdef H_THREAD 
      for (i=0;i<CPU;i++)pthread_join(tid[i], NULL); 
#endif
      if(CPU>(nocycle-cnt))CPU=nocycle-cnt;
    }/* end loop on AT points*/
}




/*####################################################*/
/*                                                    */
/*                                                    */
/*          the thread prediction function            */
/*                                                    */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/*                                                    */
/*####################################################*/

void * pred_coeur(void *input)
{
  int i,cnt;
  mega_struct *P;
  double lat,lon,decal_lat,decal_lon;

  size_t start[3],count[3];
  int indice_y,indice_x,rstatus;
  double pi;
  
  double prediction;
 
  //THE FAMOUS PI
  pi=acos(-1.0);

  //simplifie the variables name
  P=(mega_struct *) input;
  cnt=P->cnt;
  lat=P->lat[cnt];
  lon=P->lon[cnt];

  //to recale lat with the indice I must :
  decal_lat=lat+90;
  if(lon<0.0) decal_lon=360.0+lon; else decal_lon=lon;


  //find the indice in the data tab
  indice_x=(int) floor( decal_lon/(RESOLUTION) );
  indice_y=(int) floor( decal_lat /(RESOLUTION) );

  //configure the loading parameters
  start[0]=0;start[1]=indice_y;start[2]=indice_x;
  count[0]=14;count[1]=2;count[2]=2;

  //load data
  rstatus=nc_get_vara_float(P->ncid,3,start,count,P->data_amp );
  rstatus=nc_get_vara_float(P->ncid,4,start,count,P->data_phi );

  //interpolation
  rstatus=interpolation_w_mask(decal_lon,decal_lat,indice_x,indice_y,P->weight,P->data_amp);

  if(rstatus==-99)
    {
      P->Otide[cnt]=MASK;
      return input;
    }

  //initialisation of the prediction output
  P->Otide[cnt]=0;

  //prediction
  for(i=0;i<4;i++)
    {
      if(P->weight[i]!=0.0)
	{
	  prediction=predic_and_admit(P->data_amp,P->data_phi,i,P,lat,lon);
	  P->Otide[cnt]+=P->weight[i]*prediction;
	} 
    }
  return input;
}/*fin de pred_coeur*/





/*####################################################*/
/*                                                    */
/*                                                    */
/*          the thread prediction function            */
/*                                                    */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/*                                                    */
/*####################################################*/

void * extract_coeur(void *input)
{
  int i,j,cnt;
  mega_struct *P;
  double lat,lon,decal_lat,decal_lon;

  size_t start[3],count[3];
  int indice_y,indice_x,rstatus;
  double pi;
 
  //THE FAMOUS PI
  pi=acos(-1.0);

  //simplifie the variables name
  P=(mega_struct *) input;
  cnt=P->cnt;
  lat=P->lat[cnt];
  lon=P->lon[cnt];

  //to recale lat with the indice I must :
  decal_lat=lat+90;
  if(lon<0.0) decal_lon=360.0+lon; else decal_lon=lon;


  //find the indice in the data tab
  indice_x=(int) floor( decal_lon/(RESOLUTION) );
  indice_y=(int) floor( decal_lat /(RESOLUTION) );

  //configure the loading parameters
  start[0]=0;start[1]=indice_y;start[2]=indice_x; 
  count[0]=14;count[1]=2;count[2]=2;

  //load data
  rstatus=nc_get_vara_float(P->ncid,3,start,count,P->data_amp );
  rstatus=nc_get_vara_float(P->ncid,4,start,count,P->data_phi );

  //interpolation
  rstatus=interpolation_w_mask(decal_lon,decal_lat,indice_x,indice_y,P->weight,P->data_amp);

  if(rstatus==-99)
    {
      for(i=0;i<14;i++)
	{
	  P->amplitude[cnt][i]=MASK;
	  P->phase[cnt][i]=MASK;
	}
      return input;
    }

  //extraction

  for(i=0;i<14;i++)
    {
	      
      for(j=0;j<4;j++)
	{
	  if(i==4) //the M4 case//
	    {
	      if((P->data_amp[i*4+j]!=MASK)&&(P->data_amp[i*4+j]==MASK) )
		{
		  P->amplitude[cnt][i]=MASK;                                    //M4 has a different Mask
		  P->phase[cnt][i]=MASK;
		}                                    //M4 has a different Mask
	    }
	  else
	    {
	      P->amplitude[cnt][i]+=P->data_amp[i*4+j]*P->weight[j];
	      P->phase[cnt][i]+=P->data_phi[i*4+j]*P->weight[j];
	    }
	  
	}

    }
  return input;
}/*fin de extract_coeur*/




/*####################################################*/
/*                                                    */
/*                                                    */
/*      the prediction and admittance function        */
/*                                                    */
/*                                                    */
/*           Thierry LETELLIER 04Oct2005              */
/*                                                    */
/*####################################################*/

double predic_and_admit(float *data_amp,float *data_phi,int pt,mega_struct *P,double lat,double lon)
{
  int i;
  double ampli,phase;
  double prediction=0;
  double Amp,G;
  for(i=0;i<14;i++)
    {
      ampli=data_amp[i*4+pt];
      if(ampli==MASK)
	{
	  //  printf("wave %s is masked --> not inclued in prediction\n lat : %lf --- lon : %lf\n", P->spectrum[i].wave.name,lat,lon);
	  P->spectrum[i].prediction=0;
	}
      else
	{
	  phase=data_phi[i*4+pt];
	  //going the complex space
	  P->spectrum[i].Z.reel=ampli*cos(phase*-1.0*3.14/180.0);
	  P->spectrum[i].Z.imag=ampli*sin(phase*-1.0*3.14/180.0);
	  //prediction !!!!
	  P->spectrum[i].prediction=Tide_prediction(P->time[P->cnt],P->spectrum[i].wave,P->spectrum[i].Z,0,P->CTO);  
	}
    }
  /*------------------------------------------------------------------------------*/
  /* spline admittance coefficients for semi-diurnal tidal waves                  */
  /*------------------------------------------------------------------------------*/      
  // compute_admittance(P->Mat,P->Perm,P->sec_r,P->sec_i,P->res_r,P->res_i,P->spectrum, SEMI_DIURNAL,P->sindice,P->aindice); LR, change: 2008/12/17 
  compute_admittance(P->Mat,P->Perm,P->sec_r,P->sec_i,P->res_r,P->res_i,P->spectrum, SEMI_DIURNAL,P->sindice,P->aindice,P->time[P->cnt]);
  
  /*------------------------------------------------------------------------------*/
  /* spline admittance coefficients for diurnal tidal waves                       */
  /*------------------------------------------------------------------------------*/
  // compute_admittance(P->Mat,P->Perm,P->sec_r,P->sec_i,P->res_r,P->res_i,P->spectrum, DIURNAL,P->sindice,P->aindice); LR, change: 2008/12/17 
  compute_admittance(P->Mat,P->Perm,P->sec_r,P->sec_i,P->res_r,P->res_i,P->spectrum, DIURNAL,P->sindice,P->aindice,P->time[P->cnt]);
  
  /*------------------------------------------------------------------------------*/
  /* spline admittance coefficients for long period tidal waves                   */
  /*------------------------------------------------------------------------------*/  
  // compute_admittance(P->Mat,P->Perm,P->sec_r,P->sec_i,P->res_r,P->res_i,P->spectrum, LONG,P->sindice,P->aindice); LR, change: 2008/12/17 
  compute_admittance(P->Mat,P->Perm,P->sec_r,P->sec_i,P->res_r,P->res_i,P->spectrum, LONG,P->sindice,P->aindice,P->time[P->cnt]);
  
  /*------------------------------------------------------------------------------*/
  /* compute prediction for Ssa tidal wave                                        */
  /*------------------------------------------------------------------------------*/  
  tidal_potential( P->spectrum[35].wave, lat, lon, &Amp, &G,0 );
  P->spectrum[35].Z.reel=Amp*cos(G);
  P->spectrum[35].Z.imag=Amp*sin(G);
  P->spectrum[35].prediction=Tide_prediction(P->time[P->cnt],P->spectrum[34].wave,P->spectrum[34].Z,0,P->CTO);  
  
  for(i=0;i<36;i++)prediction+=P->spectrum[i].prediction;
  return(prediction);
    
}





/*####################################################*/
/*                                                    */
/*                                                    */
/*          an simple interpolation funtion           */
/*                                                    */
/*            it take in account the mask             */
/*                                                    */
/*           Thierry LETELLIER 04Oct2005              */
/*                                                    */
/*####################################################*/

int interpolation_w_mask(double lon,double lat,int indice_x,int indice_y,double *weight,float *data_amp)
{
  int i,nb_data;
  int mask_true[4],redo;
  double surface;
  //DEBUG
  double somme;

  nb_data=4;
  //init mask_true
  for(i=0;i<4;i++)mask_true[i]=1;

  for(i=0;i<4;i++){if(data_amp [i]==MASK) {mask_true[i]=0; nb_data--; } }
  
  redo=1;
  while(redo==1)
    {
      switch (nb_data)
	{
	case 0 :
	  {
	    return(-99);
	    break;
	  }
	case 1 :
	  {
	    for(i=0;i<4;i++)weight[i]=1.00*mask_true[i];
	    redo=0;
	    break;
	  }
	case 2 :
	  {
	    weight[0]=geo_d_km(lon,lat,indice_x*RESOLUTION,indice_y*RESOLUTION) * mask_true[0];
	    weight[1]=geo_d_km(lon,lat,(indice_x+1)*RESOLUTION,indice_y*RESOLUTION) * mask_true[1];
	    weight[2]=geo_d_km(lon,lat,indice_x*RESOLUTION,(indice_y+1)*RESOLUTION) * mask_true[2];
	    weight[3]=geo_d_km(lon,lat,(indice_x+1)*RESOLUTION,(indice_y+1)*RESOLUTION) * mask_true[3];
	    redo=0;
	    break;
	  }
	case 3 :
	  {
	    surface=RESOLUTION*RESOLUTION;
	    somme=( (indice_x+1)*RESOLUTION - lon )*((indice_y+1)*RESOLUTION - lat);
	    if(somme==0){mask_true[0]=0;nb_data--;break;}
	    weight[0]=surface*somme;
	
	    somme=( lon - indice_x*RESOLUTION  )*((indice_y+1)*RESOLUTION - lat);
	    if(somme==0){mask_true[0]=0;nb_data--;break;}
	    weight[1]=surface*somme;

	    somme=( (indice_x+1)*RESOLUTION - lon  )*(lat - indice_y*RESOLUTION);
	    if(somme==0){mask_true[0]=0;nb_data--;break;}
	    weight[2]=surface*somme;

	    somme=(  lon - indice_x*RESOLUTION   )*(lat - indice_y*RESOLUTION);
	    if(somme==0){mask_true[0]=0;nb_data--;break;}
	    weight[3]=surface*somme;
	    redo=0;

	    break;
	  }
	case 4 :
	  {
	    surface=RESOLUTION*RESOLUTION;
	    somme=( (indice_x+1)*RESOLUTION - lon )*((indice_y+1)*RESOLUTION - lat);
	    if(somme==0){mask_true[0]=0;nb_data--;break;}
	    weight[0]=surface*somme;
	
	    somme=( lon - indice_x*RESOLUTION  )*((indice_y+1)*RESOLUTION - lat);
	    if(somme==0){mask_true[0]=0;nb_data--;break;}
	    weight[1]=surface*somme;

	    somme=( (indice_x+1)*RESOLUTION - lon  )*(lat - indice_y*RESOLUTION);
	    if(somme==0){mask_true[0]=0;nb_data--;break;}
	    weight[2]=surface*somme;

	    somme=(  lon - indice_x*RESOLUTION   )*(lat - indice_y*RESOLUTION);
	    if(somme==0){mask_true[0]=0;nb_data--;break;}
	    weight[3]=surface*somme;
	    redo=0;
	    break;
	  }
	}
  
    }
  surface=weight[0]+weight[1]+weight[2]+weight[3];

  weight[0]/=surface;
  weight[1]/=surface;
  weight[2]/=surface;
  weight[3]/=surface;

  //DEBUG
  somme=weight[0]+weight[1]+weight[2]+weight[3];
  
  return(0);

}




/*####################################################*/
/*                                                    */
/*                                                    */
/*        compute the admittance to extend            */
/*            the prediction spectrum                 */
/*                                                    */
/*           Thierry LETELLIER 04Oct2005              */
/*                                                    */
/*####################################################*/

void compute_admittance(gsl_matrix *Mat,gsl_permutation *Perm,gsl_vector *sec_r,gsl_vector *sec_i,gsl_vector *res_r,gsl_vector *res_i,spectrum_struct *spectrum,int GROUPE,int *sindice ,int *aindice, double time)
{

  int i,n,gsl_out;
  double d2r,tau;
  int signum;
  double omega;
  int nb;

  d2r=pi/180.0;
  tau=2*24;
  n=0;

  init_admittance_coeff(sindice,aindice,&nb,GROUPE);
 
  for(i=0;i<3;i++)
    {
      omega=pulsation(spectrum[sindice[i]].wave);
      gsl_matrix_set(Mat,n,0,cos(tau*omega*d2r)*spectrum[sindice[i]].wave.Ap);
      gsl_matrix_set(Mat,n,1,sin(tau*omega*d2r)*spectrum[sindice[i]].wave.Ap);
      gsl_matrix_set(Mat,n,2,spectrum[sindice[i]].wave.Ap);
      gsl_vector_set(sec_r,n,spectrum[sindice[i]].Z.reel);
      gsl_vector_set(sec_i,n,spectrum[sindice[i]].Z.imag);
      n++;
    }	
  gsl_out=gsl_linalg_LU_decomp(Mat,Perm,&signum);
  gsl_out=gsl_linalg_LU_solve(Mat,Perm,sec_r,res_r);
  gsl_out=gsl_linalg_LU_solve(Mat,Perm,sec_i,res_i);
  
  for(i=0;i<nb;i++)
    {
      omega=pulsation(spectrum[aindice[i]].wave);
      
      spectrum[aindice[i]].Z.reel=(gsl_vector_get(res_r,0)*cos(tau*omega*d2r)+gsl_vector_get(res_r,1)*sin(tau*omega*d2r)+gsl_vector_get(res_r,2))*spectrum[aindice[i]].wave.Ap;
      spectrum[aindice[i]].Z.imag=(gsl_vector_get(res_i,0)*cos(tau*omega*d2r)+gsl_vector_get(res_i,1)*sin(tau*omega*d2r)+gsl_vector_get(res_i,2))*spectrum[aindice[i]].wave.Ap;
      spectrum[aindice[i]].prediction=Tide_prediction(time,spectrum[aindice[i]].wave,spectrum[aindice[i]].Z,0,spectrum[aindice[i]].CTO); 
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double geo_d_km(double t1,double p1,double t2,double p2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

{
  double dt,dp;
  double a1,a2,b1,b2,ux,uy,uz,vx,vy,vz,pds,angle;
  double ro;
  double pi=3.14159265358979323846;

  double dtr=pi/180.0;


  dt=t2-t1;
  dp=p2-p1;

  if ((dt == 0.0) && (dp == 0.0)) return(0.);

  a1=t1*dtr;
  b1=p1*dtr;
  a2=t2*dtr;
  b2=p2*dtr;

  ux=cos(a1)*cos(b1);
  uy=sin(a1)*cos(b1);
  uz=sin(b1);
  vx=cos(a2)*cos(b2);
  vy=sin(a2)*cos(b2);
  vz=sin(b2);

  pds=ux*vx+uy*vy+uz*vz;

  if(pds < 1.0) 
    {
    angle=acos(pds);
    /* Conversion real double de geo_mean_radius*/
    ro=(double)(6400.*angle);
    return(ro);
    }
  else return(0.);

}






//##############################################################################################################################
//##############################################################################################################################
//##############################################################################################################################
//##############################################################################################################################
//##############################################################################################################################



//##############################################################################################################################
//##############################################################################################################################
//##############################################################################################################################
//##############################################################################################################################
//##############################################################################################################################



//##############################################################################################################################
//##############################################################################################################################
//##############################################################################################################################
//##############################################################################################################################
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
//##############################################################################################################################
//##############################################################################################################################
//##############################################################################################################################
//##############################################################################################################################






double Tide_prediction(double time,tidal_wave wave,fcomplex Z,int verbose,int CTO)
/*---------------------------------------------------------------*/
/* The prediction is made from a true complex,                   */
/*      Z  =     (A cos(-G), A sin(-G) )                         */
/*---------------------------------------------------------------*/
{
  double V,f,ret;
  astro_ang_struct astro_ang;

 init_argument(time,verbose,CTO,&astro_ang);
 V=greenwhich_argument(wave,&astro_ang)+nodal_phase(wave,&astro_ang);
 f=nodal_factort(wave.formula,&astro_ang);

 ret=f*(cos(V)*Z.reel-sin(V)*Z.imag);
 return(ret); 
}

			
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void init_argument(double first_time, int verbose, int CTO,astro_ang_struct *astro_ang)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

{
  int N;

  date_t t_schureman={1,1,1900,0.0};    /* LR: 21/07/2003, -43200 -> 0 */
  date_t date58=     { 1, 1, 1958, 0.0};/* LR: 21/07/2003, -43200 -> 0 */ 
  date_t date50=     { 1, 1, 1950, 0.0};/* LR: 21/07/2003, -43200 -> 0 */
  date_t date85=     { 1, 1, 1985, 0.0};/* LR: 21/07/2003, -43200 -> 0 */
  date_t date2000=     { 1, 1, 2000, 0.0};/* LR: 21/07/2003, -43200 -> 0 */
  date_t date;
  double tj;

  
  switch (CTO)
    {
    case(0): {date =date50;break;}
    case(1): {date =date58;break;}
    case(2): {date =date85;break;}
    case(3): {date =date2000;break;}
    default: {printf("error in time reference, exit\n");exit(22);}
    }
  
  N=  julian_day(date.month,date.day,date.year)
     -julian_day(t_schureman.month,t_schureman.day,t_schureman.year);
  /* number of day elapsed between 1950 and 1900 (CNES Time and SCHUREMAN Time) */
  /*or number of day elapsed between 1958 and 1900 (CTO Time and SCHUREMAN Time) */
  /*or number of day elapsed between 1985 and 1900 (ESA Time and SCHUREMAN Time) */

  tj=((double) N + first_time/(double) 24.) /(double) 36525.0;

  astronomic_angle(tj,verbose,astro_ang);

}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void astronomic_angle(double tj, int verbose, astro_ang_struct *astro_ang)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


/*----------------------------------------------------------------------

  tj is time elapsed since 1st January 1900, 0 hour, counted in julian
  century (e.g. time in days divide by 36525)
----------------------------------------------------------------------*/

{
  double dtr,days,tt;
  date_t actual;
  double cosI,p,q;
  double t2,t4,sin2I,s2;
  double tgI2,P;
  double sh_tgn2,at1,at2;

  pi=acos(-1.0);
  dtr=pi/180.0;

/*----------------------------------------------------------------------
  Sh_n Longitude of ascending lunar node
----------------------------------------------------------------------*/

  astro_ang->sh_N=(259.1560563 -1934.1423972 *tj)*dtr;

/*----------------------------------------------------------------------
 T mean solar angle (Greenwhich time)
----------------------------------------------------------------------*/
  days=36525.0*tj;
  astro_ang->sh_T=((days - (int)days)*24.0*15.0+180.0)*dtr;

/*----------------------------------------------------------------------
 h mean solar Longitude
----------------------------------------------------------------------*/

  astro_ang->sh_h=(280.1895015 +36000.76892 *tj)*dtr;

/*----------------------------------------------------------------------
 s mean lunar Longitude
----------------------------------------------------------------------*/

  astro_ang->sh_s=(277.0256206 +481267.892 *tj)*dtr;

/*----------------------------------------------------------------------
 p1 Longitude of solar perigee
----------------------------------------------------------------------*/

  astro_ang->sh_p1=(281.2208568 +1.719175 *tj)*dtr;

/*----------------------------------------------------------------------
 p Longitude of lunar perigee
----------------------------------------------------------------------*/

  astro_ang->sh_p=(334.3837214 +4069.0322056 *tj)*dtr;

  astro_ang->sh_N =fmod(astro_ang->sh_N ,2*pi);
  astro_ang->sh_s =fmod(astro_ang->sh_s ,2*pi);
  astro_ang->sh_h =fmod(astro_ang->sh_h, 2*pi);
  astro_ang->sh_p =fmod(astro_ang->sh_p, 2*pi);
  astro_ang->sh_p1=fmod(astro_ang->sh_p1,2*pi);

  cosI=0.913694997 -0.035692561 *cos(astro_ang->sh_N);

  astro_ang->sh_I=acos(cosI);

  sin2I=sin(astro_ang->sh_I);
  sh_tgn2=tan(astro_ang->sh_N/2.0);
  
  at1=atan(1.01883*sh_tgn2);
  at2=atan(0.64412*sh_tgn2);
  
  astro_ang->sh_xi=-at1-at2+astro_ang->sh_N;

  if (astro_ang->sh_N > pi) astro_ang->sh_xi=astro_ang->sh_xi-2.0*pi;

  astro_ang->sh_nu=at1-at2;

/*----------------------------------------------------------------------
 For constituents l2 k1 k2
----------------------------------------------------------------------*/

  tgI2=tan(astro_ang->sh_I/2.0);
  P=astro_ang->sh_p-astro_ang->sh_xi;
  
  t2=tgI2*tgI2;
  t4=t2*t2;
  astro_ang->sh_x1ra=sqrt(1.0-12.0*t2*cos(2.0*P)+36.0*t4);
  
  p=sin(2.0*P);
  q=1.0/(6.0*t2)-cos(2.0*P);
  astro_ang->sh_R=atan(p/q);
  
  p=sin(2.0*astro_ang->sh_I)*sin(astro_ang->sh_nu);
  q=sin(2.0*astro_ang->sh_I)*cos(astro_ang->sh_nu)+0.3347;
  astro_ang->sh_nuprim=atan(p/q);
  
  s2=sin(astro_ang->sh_I)*sin(astro_ang->sh_I);
  p=s2*sin(2.0*astro_ang->sh_nu);
  q=s2*cos(2.0*astro_ang->sh_nu)+0.0727;
  astro_ang->sh_nusec=0.5*atan(p/q); 

  tt=(tj*36525.-18262)*24;
  getcnesdate(tt,&actual);


  if (verbose)
    {
    printf ("%d/%d/%d \n",actual.day,actual.month,actual.year);
    printf ("s: %f h: %f p: %f p1: %f \n",astro_ang->sh_s/dtr,astro_ang->sh_h/dtr,astro_ang->sh_p/dtr,astro_ang->sh_p1/dtr);
    printf ("I: %f N: %f \n",astro_ang->sh_I/dtr,astro_ang->sh_N/dtr);
    }
 
}

    
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double greenwhich_argument(tidal_wave w,astro_ang_struct *astro_ang)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------

  Returns the phase of the tidal potential relative to the Greenwhich 
  meridian (e.g. the position of the fictuous celestial body). Units are
  radian.

----------------------------------------------------------------------*/

{
  double V0;

  V0=astro_ang->sh_T*w.nT+astro_ang->sh_s*w.ns+astro_ang->sh_h*w.nh+astro_ang->sh_p*w.np+astro_ang->sh_p1*w.np1+w.shift*deg_to_rad;
  V0=fmod(V0,pi2);
  return(V0);

}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double nodal_phase(tidal_wave w,astro_ang_struct *astro_ang)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------

  Returns the phase correction u due to nodal motion. Units are radian.

----------------------------------------------------------------------*/

{
  double u;

  u=astro_ang->sh_xi*w.nksi+astro_ang->sh_nu*w.nnu0+astro_ang->sh_nuprim*w.nnu1+astro_ang->sh_nusec*w.nnu2+astro_ang->sh_R*w.R;
  return(u);

}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int julian_day(int mm,int id,int iyyy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

/*-----------------------------------------------------------------------
  tides_julianday routine returns the julian day number which begins at noon
  of the calendar date specified by month mm, day id, and year iyyy.
  positive year signifies a.d.; negative, b.c.  (remember that the
  year after 1 b.c. was 1 a.d.
  this routine has been lifted directly from the book
  press et al., numerical recipes, cambridge univ. press, 1986.
  routine handles changeover to gregorian calendar on oct. 15, 1582.
  note: to get the corresponding modified julian day,
        set mjd = julday - 2400001.

	f77 to C: 7/10/2001 (Thierry)
-----------------------------------------------------------------------*/

  {
   int igreg=15+31*(10+12*1582);
   int jy,jm,ja, tmp_iyyy;
   double temp_tides_juliandays = 0.;
   
   tmp_iyyy=iyyy;
   
   if (tmp_iyyy != 0)
/*      if (tmp_iyyy == 2000) printf("Probleme avec le passage a l annee 2000"); */
/*      else  */
       { 
       if (tmp_iyyy < 0) tmp_iyyy = tmp_iyyy + 1;
       if (mm > 2)
         {
	 jy = tmp_iyyy;
	 jm = mm + 1;
	 }
       else
         {
	 jy = tmp_iyyy - 1;
	 jm = mm + 13;
	 }
       temp_tides_juliandays=floor(365.25*jy)+floor(30.6001*jm)+id+1720995;
       if(id+31*(mm+12*tmp_iyyy) >= igreg)
         {
	 ja=floor(0.01*jy);
         temp_tides_juliandays=temp_tides_juliandays+2-ja+floor(0.25*ja);
	 }
       }/*else*/
       
    return(temp_tides_juliandays);   
    
  }/*end*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double nodal_factort(int formula, astro_ang_struct *astro_ang)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double s,f=0.;

  switch (formula)
    {

/* formule 0, solar waves */
  
    case 0:
      f=1.0;
      break;
   
/* formule 1, compound waves (78 x 78)*/
  
    case 1:
      f=nodal_factort(78,astro_ang)*nodal_factort(78,astro_ang);
      break;

/* formule 2, compound waves (78 x 0)  ===  (78)*/
  
    case 2:
      f=nodal_factort(78,astro_ang)*nodal_factort(0,astro_ang);
      break;
  
  
  
/* formule 4,  compound waves (78 x 235) */
  
    case 4:
      f=nodal_factort(78,astro_ang)*nodal_factort(235,astro_ang);
      break;

  
/* formule 5,  compound waves (78 *78 x 235) */
  
    case 5:
      f=nodal_factort(78,astro_ang)*nodal_factort(78,astro_ang)*nodal_factort(235,astro_ang);
      break;
  
/* formule 6,  compound waves (78 *78 x 0) */
  
    case 6:
      f=nodal_factort(78,astro_ang)*nodal_factort(78,astro_ang)*nodal_factort(0,astro_ang);
      break;
  
/* formule 7, compound waves (75 x 75)*/
  
    case 7:
      f=nodal_factort(75,astro_ang)*nodal_factort(75,astro_ang);
      break;
      
/* formule 8,  compound waves (78 x 0 x 235) */
  
    case 8:
      f=nodal_factort(78,astro_ang)*nodal_factort(0,astro_ang)*nodal_factort(235,astro_ang);
      break;
  
/* formule 9,  compound waves (78 x 0 x 227) */
  
    case 9:
      f=nodal_factort(78,astro_ang)*nodal_factort(0,astro_ang)*nodal_factort(227,astro_ang);
      break;
  
/* formule 10,  compound waves (78 x 227) */
  
    case 10:
      f=nodal_factort(78,astro_ang)*nodal_factort(227,astro_ang);
      break;
  
/* formule 11,  compound waves (75 x 0) */
  
    case 11:
      f=nodal_factort(75,astro_ang)*nodal_factort(0,astro_ang);
      break;
      
/* formule 12,  compound waves (78 x 78 x 78 x 0) */
  
    case 12:
      f=nodal_factort(78,astro_ang)*nodal_factort(78,astro_ang)*nodal_factort(78,astro_ang)*nodal_factort(0,astro_ang);
      break;
      
/* formule 13, compound waves (78 x 75)*/
  
    case 13:
      f=nodal_factort(78,astro_ang)*nodal_factort(75,astro_ang);
      break;
  
/* formule 14, compound waves (235 x 0)  ===  (235)*/
  
    case 14:
      f=nodal_factort(235,astro_ang)*nodal_factort(0,astro_ang);
      break;
  
/* formule 15, compound waves (235 x 75) */
  
    case 15:
      f=nodal_factort(235,astro_ang)*nodal_factort(75,astro_ang);
      break;
  
/* formule 16, compound waves (78 x 0 x 0)  ===  (78)*/
  
    case 16:
      f=nodal_factort(78,astro_ang)*nodal_factort(0,astro_ang)*nodal_factort(0,astro_ang);
      break;
      
/* formule 17,  compound waves (227 x 0) */
  
    case 17:
      f=nodal_factort(227,astro_ang)*nodal_factort(0,astro_ang);
      break;
      
/* formule 18,  compound waves (78 x 78 x 78 ) */
  
    case 18:
      f=nodal_factort(78,astro_ang)*nodal_factort(78,astro_ang)*nodal_factort(78,astro_ang);
      break;
  
/* formule 19, compound waves (78 x 0 x 0 x 0)  ===  (78)*/
  
    case 19:
      f=nodal_factort(78,astro_ang)*nodal_factort(0,astro_ang)*nodal_factort(0,astro_ang);
      break;
      

/* formule 73 */
  
    case 73:
      s=sin(astro_ang->sh_I);
      f=(2./3.-s*s)/0.5021;
      break;

/* formule 74 */
  
    case 74:
      s=sin(astro_ang->sh_I);
      f=s*s/0.1578;
      break;
  
/* formule 75 */
  
    case 75:
      s=cos (astro_ang->sh_I/2);
      f=sin (astro_ang->sh_I)*s*s/0.3800;
      break;  

/* formule 76 */
  
    case 76:
      f=sin (2*astro_ang->sh_I)/0.7214;
      break;
  
/* formule 77 */
  
    case 77:
      s=sin (astro_ang->sh_I/2);
      f=sin (astro_ang->sh_I)*s*s/0.0164;
      break;
  
/* formule 78 */
  
    case 78:
      s=cos (astro_ang->sh_I/2);
      f=s*s*s*s/0.9154;
      break;

/* formule 79 */
    
    case 79:
      s=sin(astro_ang->sh_I);
      f=s*s/0.1565;
      break;
  
/* formule 144 */
  
    case 144:
      s=sin (astro_ang->sh_I/2);
      f=(1-10*s*s+15*s*s*s*s)*cos(astro_ang->sh_I/2)/0.5873;
      break;

/* formule 149 */
  
    case 149:
      s=cos (astro_ang->sh_I/2);
      f=s*s*s*s*s*s/0.8758;
      break;

/* formule 215 */
  
    case 215:
      s=cos (astro_ang->sh_I/2);
      f=s*s*s*s/0.9154*astro_ang->sh_x1ra;
      break;
  
/* formule 227 */
  
    case 227:
      s=sin (2*astro_ang->sh_I);
      f=sqrt (0.8965*s*s+0.6001*s*cos (astro_ang->sh_nu)+0.1006);
      break;

/* formule 235 */
   
    case 235:
      s=sin (astro_ang->sh_I);
      f=sqrt (19.0444*s*s*s*s+2.7702*s*s*cos (2*astro_ang->sh_nu)+.0981);
      break;

 
    default:
      exit (1);
    }
  return(f);
}
/*-----------------------------------------------------------------------*/

void calendary(int njd,date_t *actual)

/*-----------------------------------------------------------------------*/

{
  int njul,nj,nb,nm1,nj3,m,j,ndj;
  int na,nm,nd;
  int n[12]= {31,28,31,30,31,30,31,31,30,31,30,31};

/*-----------------------------------------------------------------------
! Input = njd: nombre ecoules depuis le 1er  Janvier 1950, 0 heure
! Output= nd,nm,na: jour, mois annee calendaire
!-----------------------------------------------------------------------*/

  njul=njd+1;
  na=njul/365;
  nj=njul-na*365;
  nb=(na+1)/ 4;
  nj=nj-nb;

  if (nj >  0) goto a1;

  na=na+1949;
  nm=12;
  nd=nj+31;
  goto a9000;

  a1:
  j=na-2-nb*4;
  na=na+1950;

  if (j < 0) goto a5000;

  if (60-nj < 0)  goto  a4500;
  if (60-nj == 0) goto  a7000;
  goto a5000;

  a4500: 
  nm1=60;

  m=3;
  goto a6000;

  a5000: 
  nm1=0;

  m=1;

  a6000: 
  ndj=nm1+n[m-1];

  nj3=nj-ndj;

  if (nj3 <= 0) goto a8000;

  m=m+1;

  nm1=ndj;
  goto a6000;

  a7000: 
  nm=2;

  nd=29;

  goto a9000;

  a8000: 
  nm=m;

  nd=nj-nm1;

  a9000: 
  actual->year=na;
  actual->month=nm;
  actual->day=nd;
  actual->second=0;

}

/*-------------------------------------------*/

void getcnesdate(double t,date_t *actual)

/*-------------------------------------------*/
{
  int nday;
  float second;

/* t is time elapsed from 1/1/1950 in hours */
  
  nday=floor(t/24.0);

  calendary(nday,actual);
  second=(t-nday*24)*3600.;
  actual->second=second;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
   
double pulsation( tidal_wave wave)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* redundant with tide_wavearray_initC */

{
  double scale=36525*24.0;
  double omega;


double omega_T=  13149000.0;
double omega_s=    481267.892;
double omega_h=     36000.76892; 
double omega_p=      4069.0322056;
double omega_p1=        1.719175;

  omega = omega_T * wave.nT
        + omega_s * wave.ns
        + omega_h * wave.nh
        + omega_p * wave.np
        + omega_p1* wave.np1;
  omega /= scale;
  return(omega);
}


/*--------------------------------------------------------------------------------*/

void tidal_potential(tidal_wave wave,double lat,double lon, double *Amp, double *G,int terrestre)

/*--------------------------------------------------------------------------------*/
{
  float k2=0.3,h2=0.6;
  double C,S,a,dV,pi,dtr;

  pi=acos(-1.0);
  dtr=pi/180.0;

  if (terrestre) a=wave.Ap*h2;
  else a=wave.Ap*(1+k2-h2);
  switch (wave.nT)
    {
    case (0):

/*######################### Long period tide #########################

   potential/g = A*(1/2 -3/2 sin^2(Y)*cos(w*t+V0)

   dP/dx=  0
   dP/dy= -3*A*cos(Y)sin(Y)*cos(w*t+V0)

----------------------------------------------------------------------*/
      dV=0.0;
      C=cos(lat*dtr);
      S=sin(lat*dtr);
      *Amp = a*(0.5-1.5*S*S)/100;
      *G   = -dV*dtr;
    break;

    case (1):

/*########################### Diurnal tide ###########################

   potential/g = A*sin(2Y)*cos(w*t+V0+X)

----------------------------------------------------------------------*/

      dV=lon;
      C=cos(lat*dtr);
      S=sin(lat*dtr);
      *Amp = 2*a*S*C/100;
      *G   = -dV*dtr;
    break;

    case (2):

/*######################### Semi-diurnal tide #########################

   potential/g = A*cos^2(Y)*cos(w*t+V0+2*X)

----------------------------------------------------------------------*/

      dV=2*lon;
      C=cos(lat*dtr);
      S=sin(lat*dtr);
      *Amp = a*C*C/100;
      *G   = -dV*dtr;
    break;

/*####################### non-astronomical tide #######################
   potential/g = 0

----------------------------------------------------------------------*/
    default:
    break;

  }/*end switch */
}/* end*/

