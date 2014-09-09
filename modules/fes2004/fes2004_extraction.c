#include "fes2004_lib.h"




int fes2004_extraction (char *netcdf_filename,int nb_position,double *lat,double *lon,double **amplitude, double **phase,int nb_CPU)
{

  /*####################################################*/
  /* variable*/
  /*####################################################*/
  //classical variables
  int rstatus;
  
  //prediction typedef variables
  spectrum_struct *spectrum=NULL;
  mega_struct *P=NULL;
  
 /*####################################################*/
  /* allocation*/
  /*####################################################*/

  alloc_tide_spectrum(&spectrum);
  P=calloc(nb_CPU,sizeof(mega_struct));
  alloc_extraction_threads(P,nb_CPU);
 
  /*####################################################*/
  /*init prediction spectrum*/
  /*####################################################*/

  rstatus=init_spectrum(spectrum,99);
  if (rstatus != 0 ) print_error_2();  

  /*####################################################*/
  /* load data files*/
  /*####################################################*/

  load_netcdf_fes2004_data(netcdf_filename,P,nb_CPU);
   

  /*####################################################*/
  /*init thread struct*/
  /*####################################################*/

  init_thread_struct(nb_CPU, P, 99, lat, lon, NULL, NULL, spectrum,amplitude,phase);

  /*####################################################*/
  /*Multithreaded extraction*/
  /*####################################################*/
  
   multi_t_extraction( nb_position, nb_CPU, P);

  /*####################################################*/
  /*free memory and exit*/
  /*####################################################*/
   free_threads(P,nb_CPU);
   free(P);

   return 0;
}/*end*/		  
