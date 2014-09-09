/* Gerris - The GNU Flow Solver
 * Copyright (C) 2011 Daniel Fuster
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.  
 */

#include <fftw3.h>
#include <stdlib.h>

#include "output.h"
#include "cartesian.h"

/* GfsOutputSpectra: header */

typedef struct _GfsOutputSpectra                     GfsOutputSpectra;

struct _GfsOutputSpectra {
  /*< private >*/
  GfsOutput parent;
  GfsCartesianGrid * cgd;
  guint dir[3];

  /*< public >*/
  GfsVariable * v;
  FttVector L, pos;
  gint level, Ndim;
};


#define GFS_OUTPUT_SPECTRA(obj)            GTS_OBJECT_CAST (obj,\
							    GfsOutputSpectra, \
							    gfs_output_spectra_class ())
#define GFS_IS_OUTPUT_SPECTRA(obj)         (gts_object_is_from_class (obj,\
								      gfs_output_spectra_class ()))

GfsOutputClass * gfs_output_spectra_class  (void);

/** \beginobject{GfsOutputSpectra} */

typedef struct {
  FILE * fp;
  fftw_complex *out;
  FttVector L,kmax;
  guint n1,n2,n3;
  guint dir1,dir2; 
} Datawrite;

static void fill_cartesian_matrix ( GfsCartesianGrid * cgd, GfsVariable * v, GfsDomain * domain )
{
  guint i,j,k;
  FttVector pos;
  FttCell * cell;

  for (i = 0; i < cgd->n[0]; i++)
    for (j = 0; j < cgd->n[1]; j++)
      for (k = 0; k < cgd->n[2]; k++) {
        pos.x = cgd->x[0][i];
        pos.y = cgd->x[1][j];
        pos.z = cgd->x[2][k];

        cell = gfs_domain_locate (domain, pos, -1, NULL);
        cgd->v[k+cgd->n[2]*(i*cgd->n[1]+j)] = gfs_interpolate (cell, pos, v);
      }
}

static FttVector init_kmax (FttVector L)
{
  FttVector kmax;
  guint i;

  for (i = 0; i < 3; i++) {
    if ((&(L.x))[i] != 0) 
      (&(kmax.x))[i] = 2.*M_PI/(&(L.x))[i]; 
    else
      (&(kmax.x))[i] = 0;
  }

  return kmax;
}

static void write_spectra_1D ( Datawrite * data )
{
  guint i; 
  FttVector k;
  fputs ("# 1:kx 2:ky 3:kz 4:real 5:img\n", data->fp);

  for ( i = 0; i < data->n1; i++ ) {
    k.x = data->kmax.x*i;
    k.y = data->kmax.y*i;
    k.z = data->kmax.z*i;
    fprintf (data->fp, "%g %g %g %g %g\n",
	     k.x, k.y, k.z , data->out[i][0], data->out[i][1] );
  }
}

static void write_spectra_2D ( Datawrite * data )
{
  guint i,j; 
  gint aux;
  FttVector k;
  k.x = k.y = k.z = 0;
  fputs ("# 1:kx 2:ky 3:kz 4:real 5:img\n", data->fp);

  for ( i = 0; i < data->n1; i++ ) {
    if ( i < data->n1/2. +1 ) (&(k.x))[data->dir1] = (&(data->kmax.x))[data->dir1]*i;
    else {
      aux = i-data->n1;
      (&(k.x))[data->dir1] = (&(data->kmax.x))[data->dir1]*aux;
    }
    for ( j = 0; j < data->n2; j++ ) {
      (&(k.x))[data->dir2] = (&(data->kmax.x))[data->dir2]*j;
      fprintf (data->fp, "%g %g %g %g %g\n",
	       k.x, k.y, k.z , data->out[i*data->n2+j][0], data->out[i*data->n2+j][1] );
    }
  }
}

static void write_spectra_3D ( Datawrite * data )
{
  guint i,j,l; 
  gint aux;
  FttVector k;
  fputs ("# 1:kx 2:ky 3:kz 4:real 5:img\n", data->fp);

  for ( i = 0; i < data->n1; i++ ) {
    if ( i < data->n1/2. +1 ) k.x = data->kmax.x*i;
    else { 
      aux = i-data->n1;
      k.x = data->kmax.x*aux;
    }
    for ( j = 0; j < data->n2; j++ ) {
      if ( j < data->n2/2. +1 ) k.y = data->kmax.y*j;
      else {
        aux = j-data->n2;
        k.y = data->kmax.y*aux;
      }
      for ( l = 0; l < data->n3; l++ ) {
        k.z = data->kmax.z*(gdouble) l;
        fprintf (data->fp, "%g %g %g %g %g\n", 
		 k.x, k.y, k.z , 
		 data->out[l+data->n3*(i*data->n2+j)][0], 
		 data->out[l+data->n3*(i*data->n2+j)][1]);
      }
    }
  }
}

static gboolean output_spectra_event (GfsEvent * event, 
				      GfsSimulation * sim) 
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_spectra_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsOutputSpectra * v = GFS_OUTPUT_SPECTRA (event);
    fftw_plan p;
    Datawrite data;

    data.fp  = GFS_OUTPUT (event)->file->fp;
    data.L   = v->L;
    data.kmax = init_kmax(v->L);
    data.dir1 = v->dir[0];
    data.dir2 = v->dir[1];

    fill_cartesian_matrix( v->cgd, v->v, domain);
    switch (v->Ndim) {
      case 1: {
	data.n1 = ( v->cgd->n[v->dir[0]] / 2 ) + 1;
	data.out = fftw_malloc( sizeof(fftw_complex)*data.n1 );
	p = fftw_plan_dft_r2c_1d( v->cgd->n[v->dir[0]], v->cgd->v, data.out, FFTW_ESTIMATE);
	fftw_execute(p);
	write_spectra_1D ( &data );
	break;
      }
      case 2: {
	data.n1 = v->cgd->n[v->dir[0]];
	data.n2 = ( v->cgd->n[v->dir[1]] / 2 ) + 1;
	data.out = fftw_malloc( sizeof(fftw_complex)*v->cgd->n[v->dir[0]]*data.n2 );
	p = fftw_plan_dft_r2c_2d( v->cgd->n[v->dir[0]], v->cgd->n[v->dir[1]], 
				  v->cgd->v, data.out, FFTW_ESTIMATE);
	fftw_execute(p); 
	write_spectra_2D ( &data );
	break;
      }
    case 3: {
      data.n1 = v->cgd->n[0];
      data.n2 = v->cgd->n[1];
      data.n3 = ( v->cgd->n[2] / 2 ) + 1;
      data.out = fftw_malloc( sizeof(fftw_complex)*v->cgd->n[0]*v->cgd->n[1]*data.n3 );
      p = fftw_plan_dft_r2c_3d( v->cgd->n[0], v->cgd->n[1], v->cgd->n[2], 
				v->cgd->v, data.out, FFTW_ESTIMATE);
      fftw_execute(p); 
      write_spectra_3D ( &data );
      break;
    }
    default:
      g_assert_not_reached ();
    }
    
    fftw_destroy_plan(p);
    fftw_free ( data.out );

    return TRUE;
  }
  return FALSE;
}

static void output_spectra_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_spectra_class ())->parent_class->read) (o, fp); 
  if (fp->type == GTS_ERROR)
    return;
  GfsOutputSpectra * v = GFS_OUTPUT_SPECTRA (*o);

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (v)");
    return;
  }
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(v->v = gfs_variable_from_name (domain->variables, fp->token->str))) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  GtsFileVariable var[] = {
    {GTS_DOUBLE, "x",  TRUE, &v->pos.x},
    {GTS_DOUBLE, "y",  TRUE, &v->pos.y},
    {GTS_DOUBLE, "z",  TRUE, &v->pos.z},
    {GTS_DOUBLE, "Lx", TRUE, &v->L.x},
    {GTS_DOUBLE, "Ly", TRUE, &v->L.y},
    {GTS_DOUBLE, "Lz", TRUE, &v->L.z},
    {GTS_NONE}
  };
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type == GTS_INT) {
    v->level = atoi (fp->token->str);
    gts_file_next_token (fp);
  }
  else
    v->level = gfs_domain_depth (domain);

  guint i, j, k, size = 1;

  v->cgd = gfs_cartesian_grid_new (gfs_cartesian_grid_class ()); 

  /* number of dims of the fft */
  v->cgd->N = 3;
  v->Ndim = 0;
  k = 0;
  for (i = 0; i < 3; i++) {
    if ((&(v->L.x))[i] != 0) {
      v->Ndim++;
      v->dir[k] = i;
      k++;
    }
  }

  if (v->Ndim == 0) {
    gts_file_error (fp, "There must be at least one L component larger than 0");
    return;
  }

  /* number of points in each direction */
  v->cgd->n = g_malloc (3*sizeof (guint));  
  for (i = 0; i < 3; i++) {
    if ((&(v->L.x))[i] == 0 )
      v->cgd->n[i] = 1;
    else
      v->cgd->n[i] = pow(2,v->level);
    size *= v->cgd->n[i];
  }

  /* mesh coordinates */
  v->cgd->x = g_malloc0 (3*sizeof (gdouble *));
  for (i = 0; i < 3; i++) {
    v->cgd->x[i] = g_malloc (v->cgd->n[i]*sizeof (gdouble));
    if (v->cgd->n[i] != 1)
      for (j = 0; j < v->cgd->n[i]; j++)
        v->cgd->x[i][j] = 
	  (&(v->pos.x))[i] + (&(v->L.x))[i]*(gdouble)j/((gdouble)(v->cgd->n[i]-1)) - 0.5;
    else
      v->cgd->x[i][0] = (&(v->pos.x))[i];
  }

  /* memory data allocation */
  v->cgd->v = g_malloc0( sizeof ( gdouble ) * 2*(size/2+1)  );
}

static void output_spectra_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_output_spectra_class ())->parent_class->write) (o, fp); 
  
  GfsOutputSpectra * v = GFS_OUTPUT_SPECTRA (o);
  fprintf (fp, " %s { x = %g y = %g z = %g Lx = %g Ly = %g Lz = %g } %d",
	   v->v->name, v->pos.x, v->pos.y, v->pos.z, v->L.x, v->L.y, v->L.z, v->level);
}

static void output_spectra_destroy ( GtsObject * o )
{
  if (GFS_OUTPUT_SPECTRA (o)->cgd)
    gts_object_destroy (GTS_OBJECT (GFS_OUTPUT_SPECTRA (o)->cgd));

  (* GTS_OBJECT_CLASS (gfs_output_spectra_class ())->parent_class->destroy) (o);
}

static void output_spectra_init ( GtsObject * o )
{
  GfsOutputSpectra * v = GFS_OUTPUT_SPECTRA (o);
  v->L.x = v->L.y = v->L.z = 0.;
  v->pos.x = v->pos.y = v->pos.z = 0.;
}

static void output_spectra_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = output_spectra_event;
  klass->read =  output_spectra_read;
  klass->write = output_spectra_write;
  klass->destroy = output_spectra_destroy;
}

GfsOutputClass * gfs_output_spectra_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsOutputSpectra",
      sizeof (GfsOutputSpectra),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) output_spectra_class_init,
      (GtsObjectInitFunc) output_spectra_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsOutputSpectra} */

/* Initialize module */

const gchar gfs_module_name[] = "fft";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{ 
  gfs_output_spectra_class ();
  return NULL; 
}
