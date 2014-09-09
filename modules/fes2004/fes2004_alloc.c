#include "fes2004_lib.h"


/*####################################################*/
/*                                                    */
/*     spectrum variable allocation                   */
/*                                                    */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/*                                                    */
/*####################################################*/

void alloc_tide_spectrum(spectrum_struct **spectrum )
{
  int nb_wave;
  nb_wave=NB_WAVE_IN_SPECTRUM;
  *spectrum=(spectrum_struct *)calloc(nb_wave,sizeof(spectrum_struct));
  if(*spectrum==NULL) print_error_3("error in spectrum allocation, you may use a larger memory computer -->exit");

}



/*####################################################*/
/*                                                    */
/*     threads variable allocation                    */
/*                                                    */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/*                                                    */
/*####################################################*/

void alloc_prediction_threads(mega_struct *P,int CPU)
{
  int i;
  
  for(i=0;i<CPU;i++)
    {
      P[i].Mat=gsl_matrix_calloc(coef_spline,coef_spline);
      if(P[i].Mat==NULL) print_error_3("error in threads allocation, you may use a larger memory computer or reduce the number of CPU -->exit");

      P[i].sec_r=gsl_vector_calloc(coef_spline);
      if(P[i].sec_r==NULL) print_error_3("error in threads allocation, you may use a larger memory computer or reduce the number of CPU -->exit");

      P[i].sec_i=gsl_vector_calloc(coef_spline);
      if(P[i].sec_i==NULL) print_error_3("error in threads allocation, you may use a larger memory computer or reduce the number of CPU -->exit");

      P[i].res_r=gsl_vector_calloc(coef_spline);
      if(P[i].res_r==NULL) print_error_3("error in threads allocation, you may use a larger memory computer or reduce the number of CPU -->exit");

      P[i].res_i=gsl_vector_calloc(coef_spline);
      if(P[i].res_i==NULL) print_error_3("error in threads allocation, you may use a larger memory computer or reduce the number of CPU -->exit");

      P[i].Perm=gsl_permutation_calloc(coef_spline);
      if(P[i].Perm==NULL) print_error_3("error in threads allocation, you may use a larger memory computer or reduce the number of CPU -->exit");

      P[i].spectrum=calloc(36,sizeof(spectrum_struct));
      if(P[i].spectrum==NULL) print_error_3("error in threads allocation, you may use a larger memory computer or reduce the number of CPU -->exit");

      P[i].data_amp=malloc(56*sizeof(float));
      P[i].data_phi=malloc(56*sizeof(float));

      P[i].weight=malloc(4*sizeof(double));
      P[i].sindice=malloc(3*sizeof(int));
      P[i].aindice=malloc(10*sizeof(int));
     

    }
}

/*####################################################*/
/*                                                    */
/*     threads variable allocation                    */
/*                                                    */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/*                                                    */
/*####################################################*/

void alloc_extraction_threads(mega_struct *P,int CPU)
{
  int i;
  for(i=0;i<CPU;i++)
    {
      P[i].spectrum=calloc(36,sizeof(spectrum_struct));
      if(P[i].spectrum==NULL) print_error_3("error in threads allocation, you may use a larger memory computer or reduce the number of CPU -->exit");
      P[i].data_amp=malloc(56*sizeof(float));
      P[i].data_phi=malloc(56*sizeof(float)); 
      P[i].weight=malloc(4*sizeof(double));
    }

}

/*####################################################*/
/*                                                    */
/*              threads variable free                 */
/*                                                    */
/*                                                    */
/*           Thierry LETELLIER 06Oct2005              */
/*                                                    */
/*####################################################*/


void free_threads(mega_struct *P,int CPU)
{
  int i;  
  for(i=0;i<CPU;i++)
    {
      if(P[i].Mat!=NULL)gsl_matrix_free( P[i].Mat);
      if(P[i].sec_r!=NULL)gsl_vector_free(P[i].sec_r);
      if(P[i].sec_i!=NULL)gsl_vector_free(P[i].sec_i);
      if(P[i].res_r!=NULL)gsl_vector_free(P[i].res_r);
      if(P[i].res_i!=NULL)gsl_vector_free(P[i].res_i);
      if(P[i].Perm!=NULL)gsl_permutation_free(P[i].Perm);

      if(P[i].spectrum!=NULL)free(P[i].spectrum);
      if(P[i].data_amp!=NULL)free(P[i].data_amp);
      if( P[i].data_phi!=NULL)free( P[i].data_phi);
      if(P[i].weight!=NULL) free(P[i].weight);
      if(P[i].sindice!=NULL)free(P[i].sindice);
      if(P[i].aindice!=NULL)free(P[i].aindice);
    }
}
