#include "fes2004_lib.h"


int fes2004_prediction (char *netcdf_filename,int time_reference,int nb_position,double *lat,double *lon,double *time,double *prediction,int nb_CPU)
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
   alloc_prediction_threads(P,nb_CPU);
      
   /*####################################################*/
   /*init prediction spectrum*/
   /*####################################################*/
   
   rstatus=init_spectrum(spectrum,time_reference);
   if (rstatus != 0 )  print_error_2();  
   
  /*####################################################*/
  /* load data files*/
  /*####################################################*/

   load_netcdf_fes2004_data(netcdf_filename,P,nb_CPU);
   
  /*####################################################*/
  /*init thread struct*/
  /*####################################################*/

   init_thread_struct(nb_CPU, P, time_reference, lat, lon, time, prediction, spectrum,NULL,NULL);

  /*####################################################*/
  /*Multithreaded prediction*/
  /*####################################################*/

   multi_t_prediction(nb_position , nb_CPU, P);

  /*####################################################*/
  /*free memory and exit*/
  /*####################################################*/
   free_threads(P,nb_CPU);
   free(P);


   printf("\n------------- prediction completed -------------\n");

   return 0;
}/*end*/		  
