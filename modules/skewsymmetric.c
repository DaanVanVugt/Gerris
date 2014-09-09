/* Gerris - The GNU Flow Solver
 * Copyright (C) 2011-2012 Daniel Fuster
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

#include "simulation.h"
#include "source.h"
#include "adaptive.h"
#include "output.h"
#include "init.h"

/* GfsSkewSymmetric: Header */

typedef struct _GfsSkewSymmetric              GfsSkewSymmetric;

struct _GfsSkewSymmetric {
  /*< private >*/
  GfsSimulation parent;
  GfsVariable * velold[FTT_NEIGHBORS];

  /*< public >*/
  gdouble beta;         /* parameter to define the position of the intermediate step */
  GfsVariable * velfaces[FTT_NEIGHBORS];
};

#define GFS_SKEW_SYMMETRIC(obj)            GTS_OBJECT_CAST (obj,		\
							 GfsSkewSymmetric,	\
							 gfs_skew_symmetric_class ())
#define GFS_IS_SKEW_SYMMETRIC(obj)         (gts_object_is_from_class (obj,	\
								   gfs_skew_symmetric_class ()))

GfsSimulationClass * gfs_skew_symmetric_class  (void);

typedef struct {
  GfsVariable **velfaces , **velold , **u; 
  GfsVariable *p;
  gdouble * dt, beta; 
} FaceData;

typedef struct {
  GfsSourceDiffusion * d; 
  GfsFunction * alpha;
  FaceData * fd;
} DataDif;

/** \beginobject{GfsSkewSymmetric} */

static void gfs_skew_symmetric_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_skew_symmetric_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != '{') 
    return;

  GtsFileVariable var[] = {
    {GTS_DOUBLE, "beta", TRUE, &GFS_SKEW_SYMMETRIC (*o)->beta},
    {GTS_NONE}
  };
  gts_file_assign_variables (fp, var);
}

static void gfs_skew_symmetric_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_skew_symmetric_class ())->parent_class->write) (o, fp);

  fprintf (fp, " { beta = %g }", GFS_SKEW_SYMMETRIC (o)->beta);
}

static void gfs_skew_symmetric_run (GfsSimulation * sim);

static void gfs_skew_symmetric_class_init (GfsSimulationClass * klass) 
{
  GTS_OBJECT_CLASS (klass)->read  = gfs_skew_symmetric_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_skew_symmetric_write;
  klass->run  = gfs_skew_symmetric_run;
}

static void gfs_skew_symmetric_init (GfsSkewSymmetric * object)
{
  object->beta = 0.05;

  GfsDomain * domain = GFS_DOMAIN (object);
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++) {
    gchar * name = g_strdup_printf ("Uface%d", d);
    gchar * descr = g_strdup_printf ("%d-component of face velocity", d);
    object->velfaces[d] = gfs_domain_add_variable (domain, name, descr);
    object->velfaces[d]->units = 1.;
    g_free(name); g_free(descr);

    name = g_strdup_printf ("Ufaceold%d", d);
    descr = g_strdup_printf ("%d-component of old face velocity", d);
    object->velold[d]   = gfs_domain_add_variable (domain, name, descr);
    object->velold[d]->units = 1.;
    g_free(name); g_free(descr);
  }

//  gfs_variable_set_vector (object->velfaces, FTT_NEIGHBORS);
//  gfs_variable_set_vector (object->velold, FTT_NEIGHBORS);
}

GfsSimulationClass * gfs_skew_symmetric_class (void)
{
  static GfsSimulationClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_skew_symmetric_info = {
      "GfsSkewSymmetric",
      sizeof (GfsSkewSymmetric),
      sizeof (GfsSimulationClass),
      (GtsObjectClassInitFunc) gfs_skew_symmetric_class_init,
      (GtsObjectInitFunc) gfs_skew_symmetric_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_simulation_class ()),
        &gfs_skew_symmetric_info);
  }

  return klass;
}

static void get_face_values (FttCell * cell, FaceData * fd)
{
  GfsStateVector * s = GFS_STATE (cell);

  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++) {
    FttComponent c = d/2;
    s->f[d].un = GFS_VALUE (cell, fd->u[c])/2.;
    if (ftt_cell_neighbor (cell, d))
      s->f[d].un += GFS_VALUE (ftt_cell_neighbor (cell, d), fd->u[c])/2.;
    else
      s->f[d].un  = 0;
  }
}

static void reset_unold (FttCell * cell, FaceData * fd)
{
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++) 
    GFS_VALUE (cell, fd->velold[d]) = 0.;
  
}

static void initialize_unold (FttCell * cell, FaceData * fd)
{
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++) 
    GFS_VALUE (cell, fd->velold[d]) = GFS_VALUE (cell, fd->velfaces[d]);
}

static void get_velfaces (FttCell * cell, FaceData * fd)
{
  GfsStateVector * s = GFS_STATE (cell);
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++) {
    GFS_VALUE (cell, fd->velfaces[d]) = s->f[d].un;
    s->f[d].un = ( fd->beta + 0.5 ) * s->f[d].un - ( fd->beta - 0.5 ) * GFS_VALUE (cell, fd->velold[d]);
  }
}

static void get_cell_values (FttCell * cell, FaceData * fd)
{
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    GFS_VALUE (cell, fd->u[c]) = (GFS_VALUE (cell, fd->velfaces[2*c]) +
				  GFS_VALUE (cell, fd->velfaces[2*c + 1]))/2.;
}

static void advance_face_values (FttCell * cell, FaceData * fd)
{
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++)
    GFS_VALUE (cell, fd->velfaces[d]) = ((1.0 + fd->beta)*GFS_VALUE (cell, fd->velfaces[d]) -
					 fd->beta*GFS_VALUE (cell, fd->velold[d]));

}

/* d: direction of the face required */
/* d2: cell direction with respect to cellref */
static gdouble interpolate_value_skew (FttCell * cellref,
				       FttDirection d, 
				       FttDirection * d2, 
				       FaceData * fd)
{
  FttCell * cell;
  FttComponent c = d/2;
  if (d2)
    cell = ftt_cell_neighbor (cellref, *d2);
  else
    cell = cellref;

  if (!cell) {
    /* Symmetric BC */
    if ( d == (*d2) ) 
      return -GFS_VALUE (cellref,fd->velfaces[FTT_OPPOSITE_DIRECTION(d)]);
    else
      return GFS_VALUE (cellref,fd->velfaces[d]);
  } 

  if (!FTT_CELL_IS_LEAF (cell)) { 
    FttDirection corner[FTT_DIMENSION];
    FttCell * cell2;
    gdouble val;
#if FTT_2D
    if ( d == (*d2) ) {
      FttComponent c1 = FTT_ORTHOGONAL_COMPONENT (c);
      corner[0]=2*c1;
      corner[1]=FTT_OPPOSITE_DIRECTION(*d2);
      cell2 = ftt_cell_child_corner(cell, corner);  
      val = GFS_VALUE (cell2,fd->velfaces[d]);
      corner[0]=2*c1+1;
      cell2 = ftt_cell_child_corner(cell, corner);
      return ( val + GFS_VALUE (cell2,fd->velfaces[d]) ) / 2. ;
    }   
    else {
      corner[0]=d;
      corner[1]=FTT_OPPOSITE_DIRECTION(*d2);
      cell2 = ftt_cell_child_corner(cell, corner);
      return GFS_VALUE (cell2,fd->velfaces[d]);
    }
#else
  static FttComponent orthogonal[FTT_DIMENSION][2] = {
    {FTT_Y, FTT_Z}, {FTT_X, FTT_Z}, {FTT_X, FTT_Y}
  };
    val = 0.;
    gint i,j;
    if ( d == (*d2) ) {
      FttVector pc;
      ftt_cell_pos (cell, &pc);
      corner[0]=FTT_OPPOSITE_DIRECTION(*d2);
      for ( i = 0; i < 2; i++ ) {
        for ( j = 0; i < 2; i++ ) {
          corner[1]=2*orthogonal[c][0]+i;     
          corner[2]=2*orthogonal[c][1]+j;     
          cell2 = ftt_cell_child_corner(cell, corner);
          val += GFS_VALUE (cell2,fd->velfaces[d]);
        }
      }
      return val / 4.;
    }
    else {
      corner[0]=FTT_OPPOSITE_DIRECTION(*d2);
      corner[1]=d;
      FttComponent c2 = (*d2) / 2;
      if ( c != orthogonal[c2][0] )
        c2 = orthogonal[c2][0];
      else
        c2 = orthogonal[c2][1];
      for ( i = 0; i < 2; i++ ) {
        corner[2]=2*c2+i;
        cell2 = ftt_cell_child_corner(cell, corner);
        val += GFS_VALUE (cell2,fd->velfaces[d]);
      }
      return val / 2.;
    }
#endif 
  }
  else {
    if ( ftt_cell_level(cell) == ftt_cell_level(cellref) || d == (*d2) )
      return GFS_VALUE (cell,fd->velfaces[d]);
    else {
      FttVector pos_next, pos_ref;
      ftt_cell_pos (cellref, &pos_ref);
      ftt_cell_pos (cell, &pos_next);
      if ( ( (&(pos_ref.x))[c] < (&(pos_next.x))[c] && (d % 2) != 0 ) || 
           ( (&(pos_ref.x))[c] > (&(pos_next.x))[c] && (d % 2) == 0 )  )
        return GFS_VALUE (cell,fd->velfaces[d]);
      else 
        return ( GFS_VALUE (cell,fd->velfaces[d]) + GFS_VALUE (cell,fd->velfaces[FTT_OPPOSITE_DIRECTION(d)]) ) / 2;
    }
  }
}

/* b Adaptative boolean */
static gdouble transverse_advection (FttCell * cell, 
				     FttComponent oc,
				     FttDirection d,
				     gdouble un,
				     FaceData * fd,
				     gboolean b)
{
  gdouble uauxbot, uauxtop,size_ratio;
  gdouble vn, vntop, vnbot, vndiag;
  FttDirection daux;
  FttCell * cellnext = ftt_cell_neighbor (cell, d);
  if (!cellnext) cellnext = cell;
  size_ratio = ftt_cell_size (cell);

  if (!b) {
    size_ratio = ftt_cell_size (cellnext)/size_ratio;
    if (!FTT_CELL_IS_LEAF (cellnext))
      size_ratio /= 2.;
    vn      = interpolate_value_skew (cell,2*oc,NULL,fd);
    vntop   = interpolate_value_skew (cell,2*oc,&d  ,fd);
    vndiag  = interpolate_value_skew (cell,2*oc+1,&d ,fd);
    vnbot   = interpolate_value_skew (cell,2*oc+1,NULL,fd);
    daux    = 2*oc;
    uauxtop = interpolate_value_skew (cell, d, &daux, fd);
    daux    = 2*oc+1;
    uauxbot = interpolate_value_skew (cell, d, &daux, fd);
  } else {
    size_ratio = size_ratio/ftt_cell_size (cellnext);
    if (!FTT_CELL_IS_LEAF (cellnext))
      size_ratio *= 2.;
    daux    = FTT_OPPOSITE_DIRECTION(d);
    vn      = interpolate_value_skew (cell,2*oc,&daux, fd);
    vntop   = interpolate_value_skew (cell,2*oc,&daux, fd);
    vndiag  = interpolate_value_skew (cell,2*oc+1,NULL,fd);
    vnbot   = interpolate_value_skew (cell,2*oc,&daux ,fd);
    daux    = 2*oc;
    uauxtop = interpolate_value_skew (cell, FTT_OPPOSITE_DIRECTION(d), &daux, fd);
    daux    = 2*oc+1;
    uauxbot = interpolate_value_skew (cell, FTT_OPPOSITE_DIRECTION(d), &daux, fd);
  }

  return (uauxtop*(vn + vntop*size_ratio) - uauxbot*(vnbot + vndiag*size_ratio)) / 4.;
}

static void advection_term (FttCell * cell, FaceData * fd)
{
  gdouble un, unext, unprev;

  FttDirection d0;
  for (d0 = 0; d0 < FTT_NEIGHBORS; d0++) {

  GfsStateVector * s = GFS_STATE (cell);
  FttComponent c = d0/2;
  FttDirection d;
  gboolean cond;

  un = GFS_VALUE (cell,fd->velfaces[d0]);
  if ((d0 % 2 ) != 0 ) {
    cond = TRUE;
    d = FTT_OPPOSITE_DIRECTION (d0);
    unext     = interpolate_value_skew (cell, d,    NULL , fd);
    unprev    = interpolate_value_skew (cell, d0, &d0, fd); 
  }
  else { 
    cond = FALSE;
    d = d0;
    unext     = interpolate_value_skew (cell, d, &d, fd);
    unprev    = interpolate_value_skew (cell, FTT_OPPOSITE_DIRECTION(d), NULL, fd);
  }

  s->f[d0].v = ((un + unext)*unext - (un + unprev)*unprev) / 4.;
#if FTT_2D
  s->f[d0].v += transverse_advection (cell, 
				      FTT_ORTHOGONAL_COMPONENT (c), d, un, fd, cond);
#else  /* FTT_3D */
  static FttComponent orthogonal[FTT_DIMENSION][2] = {
    {FTT_Y, FTT_Z}, {FTT_X, FTT_Z}, {FTT_X, FTT_Y}
  };
  s->f[d0].v += transverse_advection (cell, orthogonal[c][0], d, un, fd, cond);
  s->f[d0].v += transverse_advection (cell, orthogonal[c][1], d, un, fd, cond); 
#endif
  }
}


static gdouble get_size_next ( FttCell * cell, FttDirection d )
{
  FttCell * cellnext = ftt_cell_neighbor (cell, d);
  if (!cellnext) return ftt_cell_size (cell);
  gdouble size;

  if (!FTT_CELL_IS_LEAF (cellnext))
    size = ftt_cell_size (cell) / 2.;
  else 
    size = ftt_cell_size (cellnext);

  return size;
}

static gdouble transverse_diffusion (FttCell * cell, 
				     FttComponent oc,
				     FttDirection d,
				     gdouble un,
				     FaceData * fd)
{
  gdouble uaux, size, flux = 0;
  gint i;  

  for ( i = 0; i < 2; i++ ) {
    FttDirection daux = 2*oc + i;
    uaux = interpolate_value_skew (cell, d, &daux, fd);
    size = ( ftt_cell_size (cell) + get_size_next (cell, daux) ) / 2;
    flux += (uaux - un) / size;
  }

  return flux;
}

static void diffusion_term (FttCell * cell, DataDif * data)
{
  /* fixme: I need to account for the metric */
  gdouble size, sizenext, size_ratio;
  gdouble un, unext, unprev;

  FttDirection d0;
  for (d0 = 0; d0 < FTT_NEIGHBORS; d0++) {

  FttCellFace face = gfs_cell_face(cell, d0);
  gdouble flux = 0.;  
  gdouble invdens = data->alpha ? gfs_function_face_value (data->alpha, &face) : 1.;
  gdouble visc = gfs_diffusion_cell (data->d->D, cell);

  GfsStateVector * s = GFS_STATE (cell);

  FttDirection od = FTT_OPPOSITE_DIRECTION(d0);

  un = interpolate_value_skew (cell, d0, NULL, data->fd);

  if ((d0 % 2) != 0) {
    unext    = interpolate_value_skew (cell, od , NULL, data->fd);
    unprev   = interpolate_value_skew (cell, d0 , &(d0) , data->fd); 
    sizenext = ftt_cell_size (cell); 
    size     = get_size_next (cell, d0);
  }
  else {
    unext    = interpolate_value_skew (cell, d0, &(d0), data->fd);
    unprev   = interpolate_value_skew (cell, od,      NULL,    data->fd);
    size     = ftt_cell_size (cell); 
    sizenext = get_size_next (cell, d0);
  } 
  size_ratio = ( 1. + sizenext / size ) / 2;
  flux = ( (unext - un)/sizenext - (un - unprev)/size );

  FttComponent c = d0/2;
#if FTT_2D
  FttComponent oc = FTT_ORTHOGONAL_COMPONENT (c);
  flux += size_ratio * transverse_diffusion(cell, oc, d0, un, data->fd);
#else
  static FttComponent orthogonal[FTT_DIMENSION][2] = {
    {FTT_Y, FTT_Z}, {FTT_X, FTT_Z}, {FTT_X, FTT_Y}
  };
  flux += size_ratio * transverse_diffusion(cell, orthogonal[c][0], d0, un, data->fd);
  flux += size_ratio * transverse_diffusion(cell, orthogonal[c][1], d0, un, data->fd);
#endif 

  s->f[d0].v -= invdens*visc*flux;
  }
}

static void update_vel (FttCell * cell, FaceData * fd)
{
  GfsStateVector * s = GFS_STATE (cell);
  gdouble size;

  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++) {
    size = ( ftt_cell_size (cell) + get_size_next (cell, d) ) / 2;
    GFS_VALUE (cell, fd->velfaces[d]) = (GFS_VALUE (cell, fd->velfaces[d]) + 
					 fd->beta*GFS_VALUE (cell, fd->velold[d]))/(1.+fd->beta); 
    s->f[d].un = (2*fd->beta*GFS_VALUE (cell, fd->velfaces[d]) + 
		  (0.5-fd->beta)*GFS_VALUE (cell, fd->velold[d]) - 
		  s->f[d].v*(*fd->dt)/size)/(0.5+fd->beta);
    GFS_VALUE (cell, fd->velold[d]) = GFS_VALUE (cell, fd->velfaces[d]);
    s->f[d].v = s->f[d].un;
  }
}

/* Same as in source.c used here to obtain viscosity (make it more general?) */
static GfsSourceDiffusion * source_diffusion_viscosity (GfsVariable * v)
{
  if (v->sources) {
    GSList * i = GTS_SLIST_CONTAINER (v->sources)->items;

    while (i) {
      GtsObject * o = i->data;

      if (GFS_IS_SOURCE_DIFFUSION (o))
        return GFS_SOURCE_DIFFUSION (o);
      i = i->next;
    }
  }
  return NULL;
}

static void correct_face_velocity (FttCell * cell)
{                               
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++) {
    FttCellFace face = gfs_cell_face(cell, d);
    if (GFS_FACE_FRACTION_RIGHT (&face) != 0. && face.neighbor) {

    switch (ftt_face_type (&face)) {
      case FTT_FINE_FINE:
        GFS_STATE (face.neighbor)->f[FTT_OPPOSITE_DIRECTION(face.d)].un = GFS_STATE (cell)->f[face.d].un;
        break;
      case FTT_FINE_COARSE:
        GFS_STATE (cell)->f[face.d].un = GFS_STATE (face.neighbor)->f[FTT_OPPOSITE_DIRECTION(face.d)].un;
        break;
      default:
        g_assert_not_reached ();
    }
    }
  }
}

static void obtain_face_fluxes (const FttCell * cell)
{                               
  FttCellChildren child;
  GfsStateVector * s = GFS_STATE (cell);

  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++) {
    FttCell * neighbor = ftt_cell_neighbor (cell, d);
    if (neighbor) {
      if (!FTT_CELL_IS_LEAF (neighbor)) {
        gint i, n = ftt_cell_children_direction (neighbor, FTT_OPPOSITE_DIRECTION(d), &child);
        s->f[d].v = 0;
        for (i = 0; i < n; i++)
          if (child.c[i])  
            s->f[d].v += GFS_STATE (child.c[i])->f[FTT_OPPOSITE_DIRECTION(d)].v;
        s->f[d].v /= n;
      }
      else if ((d % 2) > 0 && ftt_cell_level(cell) == ftt_cell_level(neighbor))
        s->f[d].v = GFS_STATE (neighbor)->f[FTT_OPPOSITE_DIRECTION(d)].v;
    }
    else
      s->f[d].v = 0;
  }
}

static void gfs_skew_symmetric_momentum (GfsSimulation * sim, FaceData * fd, GfsVariable **gmac)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  FttComponent c;
  FttDirection d;
  /* it is used for implementation of viscosity (improve?) */
  GfsSourceDiffusion * dif = source_diffusion_viscosity (fd->u[0]);

  gfs_domain_cell_traverse (domain,
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) advance_face_values, fd);
  
  /* boundary conditions */
  for (d = 0; d <  FTT_NEIGHBORS; d++)
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, fd->velfaces[d]);

  gfs_domain_cell_traverse (domain, 
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) advection_term, fd);

  if (dif) { 
    DataDif dd = { dif , sim->physical_params.alpha, fd };
    gfs_domain_cell_traverse (domain, 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) diffusion_term, &dd); 
  }

  /* regularize flux at faces */
  for (c = 0; c <  FTT_DIMENSION; c++)
    gfs_domain_face_bc (domain, c, fd->u[c]);

  gfs_domain_cell_traverse (domain, 
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) obtain_face_fluxes, NULL);

  gfs_domain_cell_traverse (domain, 
			    FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) update_vel, fd);

  gfs_velocity_face_sources (domain, fd->u, (*fd->dt), sim->physical_params.alpha, gmac);

  gfs_domain_cell_traverse (domain, 
                            FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1, 
                            (FttCellTraverseFunc) correct_face_velocity, NULL);
}

static void gfs_skew_symmetric_run (GfsSimulation * sim)
{
  GfsVariable * p,  * res = NULL, * gmac[FTT_DIMENSION]; 
  GfsDomain * domain;
  GSList * i;

  domain = GFS_DOMAIN (sim);

  p = gfs_variable_from_name (domain->variables, "P");

  g_assert (p);
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++) 
    gmac[c] = gfs_temporary_variable (domain);

  gfs_variable_set_vector (gmac, FTT_DIMENSION);

  gfs_simulation_refine (sim);
  gfs_simulation_init (sim);

  i = domain->variables;
  while (i) {
    if (GFS_IS_VARIABLE_RESIDUAL (i->data))
      res = i->data;
    i = i->next;
  }

  gfs_simulation_set_timestep (sim);

  GfsVariable ** u = gfs_domain_velocity (domain);
  GfsVariable ** velfaces = GFS_SKEW_SYMMETRIC(sim)->velfaces;
  GfsVariable ** velold   = GFS_SKEW_SYMMETRIC(sim)->velold;

  FaceData fd = { velfaces, velold, u, p, &sim->advection_params.dt, GFS_SKEW_SYMMETRIC(sim)->beta};

  if (sim->time.i == 0) {

    gfs_domain_cell_traverse (domain, 
                              FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) reset_unold, &fd);
    
    


    gfs_domain_cell_traverse (domain, 
        FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
        (FttCellTraverseFunc) get_face_values, &fd);
  
    gfs_mac_projection (domain,
			&sim->projection_params, 
			sim->advection_params.dt/2.,
			p, sim->physical_params.alpha, gmac, NULL);
 
    gfs_domain_cell_traverse (domain, 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) get_velfaces, &fd);

    gfs_domain_cell_traverse (domain, 
                              FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) initialize_unold, &fd);
  
  }

  while (sim->time.t < sim->time.end && sim->time.i < sim->time.iend) {
    
    gdouble tstart = gfs_clock_elapsed (domain->timer);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);

    gfs_skew_symmetric_momentum (sim, &fd, gmac);

    gfs_mac_projection (domain,
			&sim->projection_params, 
			sim->advection_params.dt/2.,
			p, sim->physical_params.alpha, gmac, NULL);

    gfs_domain_cell_traverse (domain, 
                              FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1, 
                              (FttCellTraverseFunc) correct_face_velocity, NULL);

    gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_half_do, sim); 
    gfs_domain_cell_traverse (domain, 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) get_velfaces, &fd);

    gfs_domain_cell_traverse (domain, 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) get_cell_values, &fd);

    gfs_domain_cell_traverse (domain,
			      FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_coarse_init, domain);
    gfs_simulation_adapt (sim);

    sim->time.t = sim->tnext;
    sim->time.i++;

    gfs_simulation_set_timestep (sim);
    gfs_advance_tracers (sim, sim->advection_params.dt);

    gts_range_add_value (&domain->timestep, gfs_clock_elapsed (domain->timer) - tstart);
    gts_range_update (&domain->timestep);
    gts_range_add_value (&domain->size, gfs_domain_size (domain, FTT_TRAVERSE_LEAFS, -1));
    gts_range_update (&domain->size);
  }
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gfs_event_do, sim);  
  gts_container_foreach (GTS_CONTAINER (sim->events), (GtsFunc) gts_object_destroy, NULL);

  for (c = 0; c < FTT_DIMENSION; c++) 
    gts_object_destroy (GTS_OBJECT (gmac[c]));

}

/** \endobject{GfsSkewSymmetric} */

/** \beginobject{GfsInitFaceValues} */

#define GFS_IS_INIT_FACE_VALUES(obj)         (gts_object_is_from_class (obj, \
                                              gfs_init_face_values_class ()))

GfsGenericInitClass * gfs_init_face_values_class (void);

typedef struct {
  GfsVariable * v[FTT_DIMENSION];
  GfsFunction * f[FTT_DIMENSION];
  guint n;
} VarFunc;

typedef struct {
  GfsVariable * v1, * v2;
  GfsFunction * f;
} FaceInitData;

static void init_fd (FttCellFace * face, FaceInitData * fd)
{
  if (face->d % 2 != 0) 
    GFS_VALUE (face->cell, fd->v2) = gfs_function_face_value (fd->f, face);
  else
    GFS_VALUE (face->cell, fd->v1) = gfs_function_face_value (fd->f, face);
}

static gboolean gfs_init_face_values_event (GfsEvent * event, GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_face_values_class ())->parent_class)->event) 
      (event, sim)) {
    GSList * i = GFS_INIT (event)->f;
    while (i) {
      VarFunc * vf = i->data;
      FaceInitData data;
      FttComponent c = FTT_DIMENSION;
      if (!strcmp (vf->v[0]->name, "U")) {
	if (vf->n > 1)
	  g_assert_not_implemented ();
        data.v1 = GFS_SKEW_SYMMETRIC(sim)->velfaces[0];
        data.v2 = GFS_SKEW_SYMMETRIC(sim)->velfaces[1];
        data.f  = vf->f[0];
        c  = FTT_X;
      }
      else if (!strcmp (vf->v[0]->name, "V")) {
        data.v1 = GFS_SKEW_SYMMETRIC(sim)->velfaces[2];
        data.v2 = GFS_SKEW_SYMMETRIC(sim)->velfaces[3];
        data.f  = vf->f[0];
        c  = FTT_Y;
      }
#if (!FTT_2D)
      else if (!strcmp (vf->v[0]->name, "W")) {
        data.v1 = GFS_SKEW_SYMMETRIC(sim)->velfaces[4];
        data.v2 = GFS_SKEW_SYMMETRIC(sim)->velfaces[5];
        data.f  = vf->f[0];
        c  = FTT_Z;
      }
#endif
      if (c < FTT_DIMENSION) {
	gfs_catch_floating_point_exceptions ();
	gfs_domain_face_traverse (GFS_DOMAIN (sim), c,
				  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1, 
				  (FttFaceTraverseFunc) init_fd, &data);
	gfs_restore_fpe_for_function (vf->f[0]);
      }
      i = i->next;
    }
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_face_values_class_init (GfsGenericInitClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_init_face_values_event;
}

GfsGenericInitClass * gfs_init_face_values_class (void)
{
  static GfsGenericInitClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsInitFaceValues",
      sizeof (GfsInit),
      sizeof (GfsGenericInitClass),
      (GtsObjectClassInitFunc) gfs_init_face_values_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_init_class ()), &info);
  }

  return klass;
}

/** \endobject{GfsInitFaceValues} */

/* Initialize module */

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "skewsymmetric";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  gfs_skew_symmetric_class ();
  gfs_init_face_values_class ();
  return NULL;
} 
