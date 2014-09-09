#include <netcdf.h>
#include "fes2004_lib.h"


/*####################################################*/
/*                                                    */
/*        load the FES2004 data from FES2004.nc       */
/*                                                    */
/*           Thierry LETELLIER 03Oct2005              */
/*                                                    */
/*####################################################*/
void  load_netcdf_fes2004_data(char *filename,mega_struct *P,int CPU)
{


  int i,rstatus;
 
  for(i=0;i<CPU;i++) 
    {
      rstatus=nc_open( filename, NC_NOWRITE,&P[i].ncid );
      //  if (rstatus != NC_NOERR) handle_error(rstatus);
    }
}
