/* Gerris - The GNU Flow Solver
 * Copyright (C) 2010-2012 Daniel Fuster/CNRS
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

#include "particulatecommon.h"
#include "ftt.h"

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

/* GfsBubble: Header */

typedef struct _GfsBubble GfsBubble;

struct _GfsBubble {
  /*< private >*/
  GfsParticulate parent;
  gdouble rliq;
  
  /*< public >*/
  gdouble velR, p0, R0, vol_liq;
};

#define GFS_BUBBLE(obj)            GTS_OBJECT_CAST (obj, GfsBubble, gfs_bubble_class ())
#define GFS_IS_BUBBLE(obj)         (gts_object_is_from_class (obj, gfs_bubble_class ()))

static GfsEventClass * gfs_bubble_class  (void);

/* GfsBubble: Object */
/* The radius of each bubble varies according to the Rayleigh-Plesset equation */

typedef struct {
  gdouble liqpres, liqdens;
  GfsBubble * bubble;
} RPData;

static gdouble p_state_ec (GfsBubble * bubble, gdouble rb)
{
  return bubble->p0*pow (bubble->R0/rb, 3.*1.4);
}

static int func (double t, const double y[], double f[], void * params)
{
  RPData * rp = (RPData *) params;
  f[0] = y[1];
  /* interface acceleration- incompressible RP equation */
  gdouble pbubble = p_state_ec (rp->bubble, y[0]);
  f[1] = ((pbubble - rp->liqpres)/rp->liqdens - 3./2.*y[1]*y[1])/y[0];
  return GSL_SUCCESS;
}

/* jacobian matrix */
int static jac (double t, const double y[], double *dfdy, 
		double dfdt[], void *params)
{
  RPData * rp = (RPData *) params;
  gdouble pbubble = p_state_ec (rp->bubble, y[0]);
  gdouble dddRdR  = 2.*rp->liqpres-2.*(1. + 3.*1.4)*pbubble + 3.*rp->liqdens*y[1]*y[1];
  dddRdR  /= 2.*y[0]*y[0]*rp->liqdens;
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
  gsl_matrix * m = &dfdy_mat.matrix; 
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 1, 0, dddRdR);
  gsl_matrix_set (m, 1, 1, - 3.*y[1]/y[0]);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

static gboolean gfs_bubble_event (GfsEvent * event, 
				  GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_bubble_class ())->parent_class)->event) 
      (event, sim)) {
    GfsParticle * p = GFS_PARTICLE (event);
    GfsParticulate * particulate = GFS_PARTICULATE (event);
    GfsBubble * bubble = GFS_BUBBLE (event);
    GfsDomain * domain = GFS_DOMAIN (sim);

    GfsVariable * liqpres = gfs_variable_from_name (domain->variables, "P");
  
    FttCell * cell = gfs_domain_locate (domain, p->pos, -1, NULL);
    if (cell == NULL) 
      return TRUE;
    gdouble liq_rho = sim->physical_params.alpha ? 1./
      gfs_function_value (sim->physical_params.alpha, cell) : 1.;

    FttVector pos = p->pos;
    gfs_simulation_map (sim, &pos);

    gdouble point_pres = gfs_interpolate (cell, p->pos, liqpres);

    RPData rp = { point_pres, liq_rho, bubble };

    const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;

    gsl_odeiv_step * s    = gsl_odeiv_step_alloc (T, 2);
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (1e-6, 0.0);
    gsl_odeiv_evolve * e  = gsl_odeiv_evolve_alloc (2);

    gsl_odeiv_system sys = {func, jac, 2, &rp};

    gdouble t = sim->time.t;
    gdouble t1 = t + sim->advection_params.dt;
    gdouble h = 1e-6; /* better criterion?? */
    /* variables R, dot{R} */
    gdouble y[2] = { pow(3./(4.*M_PI)*particulate->volume,1./3.) , bubble->velR };

    while (t < t1) {
      int status = gsl_odeiv_evolve_apply (e, c, s,
					   &sys, &t, t1, &h, y);
      if (status != GSL_SUCCESS) 
	g_error ("Error in the RK method");
    }

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);

    bubble->velR = y[1];
    particulate->volume = 4./3.*M_PI*y[0]*y[0]*y[0];

    return TRUE;
  }
  return FALSE;
} 

static void gfs_bubble_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_bubble_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_bubble_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  GfsBubble * p = GFS_BUBBLE (*o);
  GfsParticulate * part = GFS_PARTICULATE (*o);
  gdouble L = gfs_object_simulation (*o)->physical_params.L;
    
  p->vol_liq = 0;
  p->R0 = pow (part->volume*3./(4.*M_PI), 1./3.);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (radial velocity)");
    return;
  }
  p->velR = atof (fp->token->str)/L;
  gts_file_next_token (fp);
  
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (reference pressure)");
    return;
  }
  p->p0 = atof (fp->token->str)*L;
  gts_file_next_token (fp);
}

static void gfs_bubble_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_bubble_class ())->parent_class->write) (o, fp); 
  GfsBubble * p = GFS_BUBBLE (o);
  gdouble L = gfs_object_simulation (o)->physical_params.L;
  fprintf (fp, " %g %g", p->velR*L, p->p0/L);
}

static void gfs_bubble_class_init (GfsEventClass * klass)
{
  klass->event = gfs_bubble_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_bubble_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_bubble_write;
}


static GfsEventClass * gfs_bubble_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bubble_info = {
      "GfsBubble",
      sizeof (GfsBubble),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_bubble_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particulate_class ()),
				  &gfs_bubble_info);
  }
  return klass;
}

/* GfsBubbleFraction: header */

typedef struct _GfsBubbleFraction                GfsBubbleFraction;

struct _GfsBubbleFraction {
  /*< private >*/
  GfsParticulateField parent;

  /*< public >*/
  gdouble rliq;
  GfsFunction * kernel_function;
};

#define GFS_BUBBLE_FRACTION(obj)            GTS_OBJECT_CAST (obj,	\
							     GfsBubbleFraction, \
							     gfs_bubble_fraction_class ())
#define GFS_IS_BUBBLE_FRACTION(obj)         (gts_object_is_from_class (obj, \
								       gfs_bubble_fraction_class ()))

GfsVariableClass * gfs_bubble_fraction_class  (void);

typedef struct {
  gdouble correction;
  GfsBubble * bubble;
  GfsVariable * v;
  GfsBubbleFraction * bf;
} BubbleData;

typedef struct {
  FttVector * pos;
  gdouble distance; 
} CondData;

/** \beginobject{GfsBubbleFraction} */
/* do it more general? */

static void distance_normalization (FttVector * pos1, GfsParticulate * p)
{
  gdouble rb = pow (3.*p->volume/(4.*M_PI), 1./3.);
  FttVector * pos2 = &(GFS_PARTICLE (p)->pos);
  pos1->x = (pos1->x - pos2->x)/rb;
  pos1->y = (pos1->y - pos2->y)/rb;
  pos1->z = 0.;
#if !FTT_2D
  pos1->z = (pos1->z - pos2->z)/rb;
#endif
}

static void voidfraction_from_bubbles (FttCell * cell, BubbleData * p)
{
  FttVector pos;
  ftt_cell_pos (cell, &pos);
  distance_normalization (&pos, GFS_PARTICULATE (p->bubble));
  GFS_VALUE (cell, p->v) += GFS_PARTICULATE (p->bubble)->volume*
    gfs_function_spatial_value (p->bf->kernel_function, &pos)/p->correction;
}

static void kernel_volume (FttCell * cell, BubbleData * p)
{
  gdouble cellvol = gfs_cell_volume (cell, p->v->domain);

  p->bubble->vol_liq += cellvol;

  /* correction term to make a discretely conservative kernel */
  FttVector pos;
  ftt_cell_pos (cell, &pos);
  distance_normalization (&pos, GFS_PARTICULATE (p->bubble));
  p->correction += gfs_function_spatial_value (p->bf->kernel_function, &pos)*cellvol;
}

static gboolean cond_bubble (FttCell * cell, gpointer data)
{
  CondData * p = data;
  FttVector pos;
  ftt_cell_pos (cell, &pos);
  gdouble radeq;
  gdouble size = ftt_cell_size (cell)/2.;

#if FTT_2D
  radeq = size*sqrt(2.);
#else
  radeq = size*sqrt(3.);
#endif /* 3D */

  if (ftt_vector_distance (&pos, p->pos) - radeq <= p->distance) 
    return TRUE;
    
  /* Check also if the bubble is inside the cell*/
  if (p->pos->x > pos.x + size || p->pos->x < pos.x - size ||
      p->pos->y > pos.y + size || p->pos->y < pos.y - size 
#if !FTT_2D
      || p->pos->z > pos.z + size || p->pos->z < pos.z - size 
#endif
      )
    return FALSE;

  return TRUE;
}

static gboolean bubble_fraction_event (GfsEvent * event, 
				       GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_particulate_field_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsVariable * v = GFS_VARIABLE (event);
    GfsParticulateField * pfield = GFS_PARTICULATE_FIELD (v);
    GfsBubbleFraction * bf = GFS_BUBBLE_FRACTION (event);
    
    /* Reset variable */
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_reset, v);
    /* Loop over the list of particles in the selected object */
    GSList * i = GFS_EVENT_LIST (pfield->plist)->list->items;
    while (i) {
      GfsBubble * bubble = GFS_BUBBLE (i->data);
      bubble->vol_liq = 0;
      bubble->rliq = pow (GFS_PARTICULATE(i->data)->volume*3./(4.*M_PI), 1./3.)*bf->rliq;
      BubbleData p = { 0, bubble, v, bf };
      CondData cd = { &GFS_PARTICLE (i->data)->pos, bubble->rliq };
      gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                                          (FttCellTraverseFunc) kernel_volume, &p,
                                          cond_bubble, &cd);
      gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
                                          (FttCellTraverseFunc) pfield->voidfraction_func, &p,
                                          cond_bubble, &cd);
      i = i->next;
    }
    return TRUE;
  }
  return FALSE;
}

static void bubble_fraction_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_BUBBLE_FRACTION (o)->kernel_function));

  (* GTS_OBJECT_CLASS (gfs_bubble_fraction_class ())->parent_class->destroy) (o); 
}

static void bubble_fraction_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_bubble_fraction_class ())->parent_class->read) (o, fp); 
  if (fp->type == GTS_ERROR)
    return;

  GfsBubbleFraction * b = GFS_BUBBLE_FRACTION (*o);

  if (fp->type != '{') {
    gts_file_error (fp, "expecting an opening brace");
    return;
  }
  fp->scope_max++;
  gts_file_next_token (fp);

  while (fp->type != GTS_ERROR && fp->type != '}') {
    if (fp->type == '\n') {
      gts_file_next_token (fp);
      continue;
    }
    if (fp->type != GTS_STRING) {
      gts_file_error (fp, "expecting a keyword");
      return;
    }  
    else if (!strcmp (fp->token->str, "rkernel")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting '='");
	return;
      }    
      gts_file_next_token (fp);
      b->rliq = atof (fp->token->str);
      gts_file_next_token (fp);
    }
    else if (!strcmp (fp->token->str, "kernel")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
	gts_file_error (fp, "expecting '='");
	return;
      }
      gts_file_next_token (fp);
      gfs_function_read (b->kernel_function, gfs_object_simulation (*o), fp);
    }
    else {
      gts_file_error (fp, "unknown keyword `%s'", fp->token->str);
      return;
    }
  }
  if (fp->type == GTS_ERROR)
    return;
  if (fp->type != '}') {
    gts_file_error (fp, "expecting a closing brace");
    return;
  }
  fp->scope_max--;
  gts_file_next_token (fp);
}

static void bubble_fraction_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_bubble_fraction_class ())->parent_class->write) (o, fp);
  fprintf (fp, " { rkernel = %g ", GFS_BUBBLE_FRACTION (o)->rliq);
  fputs (" kernel =", fp);
  gfs_function_write (GFS_BUBBLE_FRACTION (o)->kernel_function, fp);
  fputc ('}', fp);
}

static void bubble_fraction_init (GfsVariable * v)
{
  v->units = 0.;
  GFS_PARTICULATE_FIELD (v)->voidfraction_func = voidfraction_from_bubbles;
  GFS_BUBBLE_FRACTION (v)->kernel_function = gfs_function_new (gfs_function_spatial_class (), 0.);
}

static void bubble_fraction_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = bubble_fraction_event;
  klass->destroy = bubble_fraction_destroy;
  klass->read =  bubble_fraction_read;
  klass->write = bubble_fraction_write;
}

GfsVariableClass * gfs_bubble_fraction_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bubble_fraction_info = {
      "GfsBubbleFraction",
      sizeof (GfsBubbleFraction),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) bubble_fraction_class_init,
      (GtsObjectInitFunc) bubble_fraction_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particulate_field_class ()),
				  &gfs_bubble_fraction_info);
  }

  return klass;
}

/** \endobject{GfsBubbleFraction} */

/* GfsBubbleFractionDt: header */

#define GFS_IS_BUBBLE_FRACTION_DT(obj)         (gts_object_is_from_class (obj, \
						     gfs_bubble_fraction_dt_class ()))

GfsVariableClass * gfs_bubble_fraction_dt_class (void);

/** \beginobject{GfsBubbleFractionDt} */

static void dVpdt_from_particles (FttCell * cell, BubbleData * p)
{
  gdouble rad = pow (3.0*GFS_PARTICULATE (p->bubble)->volume/(4.0*M_PI), 1./3.);
  GFS_VALUE (cell, p->v) += 3.0*GFS_PARTICULATE (p->bubble)->volume
    *p->bubble->velR/(p->bubble->vol_liq*rad);
}

static void bubble_fraction_dt_init (GtsObject * o)
{
  GFS_PARTICULATE_FIELD (o)->voidfraction_func = dVpdt_from_particles;
}

GfsVariableClass * gfs_bubble_fraction_dt_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_bubble_fraction_dt_info = {
      "GfsBubbleFractionDt",
      sizeof (GfsBubbleFraction),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) bubble_fraction_dt_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS ( gfs_bubble_fraction_class ()),
				  &gfs_bubble_fraction_dt_info);
  }

  return klass;
}

/** \endobject{GfsBubbleFractionDt} */

/* Initialize module */

const gchar gfs_module_name[] = "bubbles";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{ 
  gfs_bubble_class ();
  gfs_bubble_fraction_class ();
  gfs_bubble_fraction_dt_class ();
  return NULL; 
}
