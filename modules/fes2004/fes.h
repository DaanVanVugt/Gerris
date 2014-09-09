
// this file is a hearder file that prototype the libfes.a functions
// The inputs are :
//    char an string that contain the ath to load the NETCDF SOLUTION
//    int time reference 0->01/01/1950   1->01/01/1958   2->01/01/1985   3->01/0/2000
//    int nb_position the number of point the predict and or to extract
//    double lat lon time 3 tabulars(nb_position) that contain the position and time for the prediction extraction
// The returns are :
//    double prediction  tabular(nb_position) that contains the prediction in meters 
//    double amplitude phase 2 dimension tabular[nb_position][nb_wave=14] that contain the amplitude and the phase_lag extracted for each point and all waves

// the int nb_cpu can be use if chachc in the source code src/fes2004_kernel.c
// and activate the multithread mode ... read the end of the README file

//extern int fes2004_prediction (char *netcdf_filename,int time_reference,int nb_position,double *lat,double *lon,double *time,double *prediction,int nb_CPU);
//extern int fes2004_extraction (char *netcdf_filename,int nb_position,double *lat,double *lon,double **amplitude,double **phase,int nb_CPU);

//  prototypes



extern int fes2004_prediction (char *,int ,int ,double *,double *,double *,double *,int );
extern int fes2004_extraction (char *,int ,double *,double *,double **,double **,int );
