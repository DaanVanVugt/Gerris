/* Gerris - The GNU Flow Solver
 * Copyright (C) 2011 Sebastien Delaux, National Institute of Water
 * and Atmospheric Research
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

#include "solid.h"
#include "output.h"
#include "vof.h"

/* This implementation of submarine landslide initial condition is directly
   inpsired from the software TOPICS developped by P. Watts.
   Sources are that used in the version 1.2 of TOPICS, last modified in
   August 2009 P. Watts.
*/

/* GfsInitSubmarineLandslide: Header */

typedef struct _GfsInitSubmarineLandslide         GfsInitSubmarineLandslide;

struct _GfsInitSubmarineLandslide {
  /*< private >*/
  GfsGenericInit parent;

  /*< public >*/
  GfsVariable * v;
  GfsVariable * vu, * vv;

  gdouble xo, yo; /* initial x/y-axis mass failure center */
  gdouble alpha;  /* CCW angle of north in degrees from grid top */

  gdouble depth;  /* d initial depth of the middle of slide (m) */
  gdouble theta;  /* mean slope along failure plane (degrees) */
  gdouble length; /* b  initial slide length during failure (m) */
  gdouble thick;  /* t Enter maximum initial slide thickness (m) */
  gdouble width;  /* w Enter maximum initial slide width (m) */
  gdouble gamma;  /* specific density */
  gdouble vol;    /* Volume of the slide */

  gdouble eta;
  gdouble lambda;
  gdouble so;
};

#define GFS_INIT_SUBMARINE_LANDSLIDE(obj) GTS_OBJECT_CAST (obj, \
				          GfsInitSubmarineLandslide,\
			     		  gfs_init_submarine_landslide_class ())
#define GFS_IS_INIT_SUBMARINE_LANDSLIDE(obj)(gts_object_is_from_class (obj,\
					 gfs_init_submarine_landlside_class ()))

static GfsEventClass * gfs_init_submarine_landslide_class  (void);


/* GfsInitSubmarineSlump: Header */

typedef struct _GfsInitSubmarineSlump         GfsInitSubmarineSlump;

struct _GfsInitSubmarineSlump {
  /*< private >*/
  GfsGenericInit parent;

  /*< public >*/
  GfsVariable * v;
  GfsVariable * vu, * vv;

  gdouble xo, yo; /* initial x/y-axis mass failure center */
  gdouble alpha;  /* CCW angle of north in degrees from grid top */

  gdouble depth;  /* d initial depth of the middle of slide (m) */
  gdouble theta;  /* mean slope along failure plane (degrees) */
  gdouble length; /* b  initial slide length during failure (m) */
  gdouble thick;  /* t Enter maximum initial slide thickness (m) */
  gdouble width;  /* w Enter maximum initial slide width (m) */
  gdouble gamma;  /* specific density */

  gdouble dist;   /* distance travelled by the center of gravity */

  gdouble lambda;
  gdouble eta;
  gdouble so;
};

#define GFS_INIT_SUBMARINE_SLUMP(obj)    GTS_OBJECT_CAST (obj, \
					 GfsInitSubmarineSlump, \
					 gfs_init_submarine_slump_class ())
#define GFS_IS_INIT_SUBMARINE_SLUMP(obj) (gts_object_is_from_class (obj,\
					 gfs_init_submarine_landlside_class ()))


static GfsEventClass * gfs_init_submarine_slump_class  (void);


/* GfsInitSubaerialLandslide: Header */

typedef struct _GfsInitSubaerialLandslide     GfsInitSubaerialLandslide;

struct _GfsInitSubaerialLandslide {
  /*< private >*/
  GfsGenericInit parent;

  /*< public >*/
  GfsVariable * v;
  GfsVariable * vu, * vv;

  gdouble xo, yo; /* initial x/y-axis mass failure center */
  gdouble alpha;  /* CCW angle of north in degrees from grid top */

  gdouble vol;    /* Landslide volume within water */
  gdouble width;  /* Landslide width at shoreline */
  gdouble depth;  /* Water depth near end of landslide (m)*/
  
  gdouble ut;     /* Landslide velocity at shoreline */
  gdouble so;     /* Landslide length in water */
  gdouble to;     /* Landslide runout time in water */

  gdouble lambda;
  gdouble eta;
};

#define GFS_INIT_SUBAERIAL_LANDSLIDE(obj) GTS_OBJECT_CAST (obj, \
					  GfsInitSubaerialLandslide,\
					  gfs_init_subaerial_landslide_class ())
#define GFS_IS_INIT_SUBAERIAL_LANDSLIDE(obj)(gts_object_is_from_class (obj,\
					 gfs_init_subaerial_landlside_class ()))

static GfsEventClass * gfs_init_subaerial_landslide_class  (void);


/* GfsInitPyroclastic: Header */

typedef struct _GfsInitPyroclastic         GfsInitPyroclastic;

struct _GfsInitPyroclastic {
  /*< private >*/
  GfsGenericInit parent;

  /*< public >*/
  GfsVariable * v;
  GfsVariable * vu, * vv;

  gdouble xo, yo; /* xo,yo position controls source location */
  gdouble alpha;  /* CCW angle of north in degrees from grid top */

  gdouble depth;  /* Water depth near end of flow */
  gdouble vol;    /* Submerged volume (m^3) */
  gdouble ut;     /* Flow velocity at shoreline (m/s) - control */
  gdouble so;     /* Flow runout length under water (m) */
  gdouble to;     /* Flow runout time under water (s) */
  gdouble width;      /* Flow width at the shoreline (m) */

  gdouble lambda;
  gdouble eta;
};

#define GFS_INIT_PYROCLASTIC(obj)       GTS_OBJECT_CAST (obj, \
					GfsInitPyroclastic, \
					gfs_init_pyroclastic_class ())
#define GFS_IS_INIT_PYROCLASTIC(obj)    (gts_object_is_from_class (obj,	\
				        gfs_init_pyroclastic_class ()))

static GfsEventClass * gfs_init_pyroclastic_class  (void);


/* GfsInitSubmarineLandslide: Object */

static void gfs_init_submarine_landslide_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_init_submarine_landslide_class ())
   ->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a variable name");
    return;
  }
  GfsInitSubmarineLandslide * l = GFS_INIT_SUBMARINE_LANDSLIDE (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(l->v = gfs_domain_get_or_add_variable (domain, fp->token->str, NULL))) {
    gts_file_error (fp, "`%s' is a reserved keyword", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
  gchar * vname = g_strconcat (l->v->name, "_U", NULL);
  l->vu = gfs_domain_get_or_add_variable (domain, vname, NULL);
  vname = g_strconcat (l->v->name, "_V", NULL);
  l->vv = gfs_domain_get_or_add_variable (domain, vname, NULL);

  GtsFileVariable var[] = {
    {GTS_DOUBLE, "x",         TRUE, &l->xo},
    {GTS_DOUBLE, "y",         TRUE, &l->yo},
    {GTS_DOUBLE, "alpha",     TRUE, &l->alpha},
    {GTS_DOUBLE, "depth",     TRUE, &l->depth},
    {GTS_DOUBLE, "theta",     TRUE, &l->theta},
    {GTS_DOUBLE, "length",    TRUE, &l->length},
    {GTS_DOUBLE, "thickness", TRUE, &l->thick},
    {GTS_DOUBLE, "width",     TRUE, &l->width},
    {GTS_DOUBLE, "volume",    TRUE, &l->vol},
    {GTS_DOUBLE, "gamma",     TRUE, &l->gamma},
    {GTS_NONE}
  };
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  /* Checks the validity of the model */
  if ( l->theta > 30 )
    g_warning ("The incline angle theta >30 degrees: %g. This will make the amplitude inaccurate!",
	       l->theta);
  if ( l->depth/l->length < 0.12 )
    g_warning ("The ratio depth/length <0.12: %g. This will make the amplitude inaccurate!",
	       l->depth/l->length);
  if ( l->thick/l->length > 0.2 )
    g_warning ("The ratio thickness/length > 0.2: %g. This will make the amplitude inaccurate!",
	       l->thick/l->length);
  if ( l->thick/l->depth > 3.33 )
    g_warning ("The ratio thickness/depth > 3.33: %g. This will make the amplitude inaccurate!",
	       l->thick/l->depth);
  if ( l->width/l->length < 0.06 )
    g_warning ("The ratio width/length < 0.06: %g. This will make the amplitude inaccurate!", 
	       l->width/l->length);
  if ( l->width/l->length > 1.0 )
    g_warning ("The ratio width/length > 1.0: %g. This will make the amplitude inaccurate!", 
	       l->width/l->length);

  gdouble g = GFS_SIMULATION(domain)->physical_params.g;
  gdouble sint = sin(l->theta*M_PI/180.);
  gdouble gmo = l->gamma - 1.;
  gdouble ao = g*sint*gmo / (l->gamma + 1.); /* Characteristic acceleration */
  gdouble ut = sqrt(0.5*g*l->length*M_PI*sint*gmo); /* Characteristic velocity*/
  l->so = pow(ut, 2.) / ao; /* Characteristic displacement */
  gdouble to = ut / ao; /* Characteristic time-scale */

  l->lambda = to*sqrt(g*l->depth); /* Initial wavelength of tsunami */

  gdouble hao = l->lambda / l->length;
  if ( hao < 1.0 )
    g_warning ("The Hammack number is < 1: %g. This will make the amplitude inaccurate\n", hao);

  gdouble sg = l->so*sint / l->depth;
  if ( sg > 0.35 )
    g_warning ("The Submergence number is > 0.35: %g. This will make the amplitude inaccurate\n", 
	       sg);

  /* Calculate tsunami initial amplitude */
  l->eta = 0.723*l->so * (4.772e-02 - 3.559e-02*sint + 8.13e-03*sint*sint) *
    (l->thick / l->length) * pow( l->length*sint / l->depth, 1.25) *
    1.18*(1.0 - exp(-2.2027*gmo));

  gdouble term = l->eta / (l->so * pow(sint, 1.5));
  if ( term > 0.2 )
    g_warning ("The term (eta/so*sinq^1.5) >0.2 : %g. This will make the amplitude inaccurate!", 
	       term);
}

static void gfs_init_submarine_landslide_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_init_submarine_landslide_class ())
   ->parent_class->write) (o, fp);
  GfsInitSubmarineLandslide * l = GFS_INIT_SUBMARINE_LANDSLIDE (o);
  fprintf (fp, " %s {\n", l->v->name);
  fprintf (fp,
	   "  x = %g y = %g alpha = %g\n"
	   "  depth = %g theta = %g\n"
	   "  length = %g width = %g\n"
	   "  thickness = %g volume = %g\n"
	   "  gamma = %g\n"
	   "}",
	   l->xo, l->yo, l->alpha,
	   l->depth, l->theta,
	   l->length, l->width,
	   l->thick, l->vol,
	   l->gamma);
}

static void init_submarine_landslide (FttCell * cell, GfsInitSubmarineLandslide * l)
{
  GfsSimulation * sim = gfs_object_simulation (l);

  gdouble cosa = cos(l->alpha*M_PI/180.);
  gdouble sina = sin(l->alpha*M_PI/180.);
  gdouble cost = cos(l->theta*M_PI/180.);
  gdouble sint = sin(l->theta*M_PI/180.);
  gdouble tant = tan(l->theta*M_PI/180.);

  FttVector p, o;
  gfs_cell_cm (cell, &p);
  o.x = l->xo; o.y = l->yo;
  gfs_simulation_map (sim, &o);
  gdouble L = sim->physical_params.L;
  p.x = (p.x - o.x)*L; p.y = (p.y - o.y)*L;
  FttVector q;
  q.x = -sina*p.x + cosa*p.y;
  q.y =- cosa*p.x - sina*p.y;

  gdouble cut = 200. * l->width;

  gdouble g = sim->physical_params.g; /* Gravity */

  gdouble xg = (l->depth + l->thick/cost) / tant;
  gdouble xmin = 0.95*( (xg + 0.4338*l->so*cost) - xg);
  gdouble nmin = -1.2*2.1*l->eta;
  gdouble nmax = 0.64*l->eta*(0.8 + 0.2*l->depth / (l->length*sint));
  
  gdouble term = 1.0 - exp( -2.0906*(l->width / l->lambda) *
			    (1.0 + 1.0903*(l->width / l->lambda)) );

  /* Initialize tsunami shape */
  GFS_VALUE (cell, l->v) = term * 
    (nmin*exp(-pow(nmin*(q.x - xmin) / (l->lambda * nmax), 2.)) +
     nmax * exp(-pow((q.x - xmin - 0.5*l->lambda) / (l->lambda), 2.))) *
    pow(2. / ( exp(3.*term*q.y / l->width) + exp(-3.*term*q.y / l->width)), 2.);
  if ( fabs(q.y) > cut) {
    GFS_VALUE (cell, l->v) *= exp(-pow(5.*(fabs(q.y) - cut) / cut, 2.));
  }

  /* Initialize velocity  from estimated linear wave quantities */
  gdouble kappa = 2.*M_PI / l->lambda;
  term = kappa*l->depth;
  gdouble omega = sqrt(g*kappa*tanh(term));

  if (GFS_VALUE (cell, l->v) > 0.) {
    gdouble utot = GFS_VALUE (cell, l->v)*g*kappa*cosh(0.469*term) / 
      (omega*cosh(term));
    GFS_VALUE (cell, l->vu) = -utot*sina;
    GFS_VALUE (cell, l->vv) = utot*cosa;
  }
  else {
    GFS_VALUE (cell, l->vu) = 0.;
    GFS_VALUE (cell, l->vv) = 0.;
  }
}

static gboolean gfs_init_submarine_landslide_event (GfsEvent * event,
						    GfsSimulation * sim)
{
  if ((*GFS_EVENT_CLASS (GTS_OBJECT_CLASS(gfs_init_submarine_landslide_class ())
			 ->parent_class)->event) (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim),
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) init_submarine_landslide,
			      event);
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_submarine_landslide_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_init_submarine_landslide_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_submarine_landslide_write;
  klass->event = gfs_init_submarine_landslide_event;
}

static GfsEventClass * gfs_init_submarine_landslide_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_submarine_landslide_info = {
      "GfsInitSubmarineLandslide",
      sizeof (GfsInitSubmarineLandslide),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_init_submarine_landslide_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_submarine_landslide_info);
  }

  return klass;
}

/* GfsInitSubmarineSlump: Object */

static void gfs_init_submarine_slump_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_init_submarine_slump_class ())
   ->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a variable name");
    return;
  }
  GfsInitSubmarineSlump * l = GFS_INIT_SUBMARINE_SLUMP (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(l->v = gfs_domain_get_or_add_variable (domain, fp->token->str, NULL))) {
    gts_file_error (fp, "`%s' is a reserved keyword", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
  gchar * vname = g_strconcat (l->v->name, "_U", NULL);
  l->vu = gfs_domain_get_or_add_variable (domain, vname, NULL);
  vname = g_strconcat (l->v->name, "_V", NULL);
  l->vv = gfs_domain_get_or_add_variable (domain, vname, NULL);

  GtsFileVariable var[] = {
    {GTS_DOUBLE, "x",         TRUE, &l->xo},
    {GTS_DOUBLE, "y",         TRUE, &l->yo},
    {GTS_DOUBLE, "alpha",     TRUE, &l->alpha},
    {GTS_DOUBLE, "depth",     TRUE, &l->depth},
    {GTS_DOUBLE, "theta",     TRUE, &l->theta},
    {GTS_DOUBLE, "length",    TRUE, &l->length},
    {GTS_DOUBLE, "thickness", TRUE, &l->thick},
    {GTS_DOUBLE, "width",     TRUE, &l->width},
    {GTS_DOUBLE, "distance",  TRUE, &l->dist},
    {GTS_DOUBLE, "gamma",     TRUE, &l->gamma},
    {GTS_NONE}
  };
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  if ( l->depth/l->length < 0.02 )
    g_warning ("The initial depth is too small for the physics!");
  if ( l->theta > 30 )
    g_warning ("The incline angle theta >30 degrees: %g. This will make the amplitude inaccurate!",
	       l->theta);
  if ( l->depth/l->length < 0.12 )
    g_warning ("The ratio depth/length <0.12: %g. This will make the amplitude inaccurate!",
	       l->depth/l->length);
  if ( l->thick/l->length > 0.2 )
    g_warning ("The ratio thickness/length > 0.2: %g. This will make the amplitude inaccurate!",
	       l->thick/l->length);
  if ( l->thick/l->depth > 3.33 )
    g_warning ("The ratio thickness/depth > 3.33: %g. This will make the amplitude inaccurate!", 
	       l->thick/l->depth);
  if ( l->width/l->length < 0.25 )
    g_warning ("The ratio width/length < 0.25: %g. This will make the amplitude inaccurate!", 
	       l->width/l->length);
  if ( l->width/l->length > 2.0 )
    g_warning ("The ratio width/length > 2.0: %g. This will make the amplitude inaccurate!", 
	       l->width/l->length);

  gdouble g = GFS_SIMULATION(domain)->physical_params.g; /* Gravity */
  gdouble gmo = l->gamma - 1.;
  gdouble sint = sin(l->theta*M_PI/180.);
  gdouble r = 0.125*pow(l->length, 2.) / l->thick + l->thick / 2.0;
  gdouble dphi = l->dist / r;

  if ( dphi > 0.53 )
    g_warning ("The angular motion dphi > 0.53 radians: %g. This will make the amplitude inaccurate!", dphi);
  if ( r/l->length > 2.0 || r/l->length < 1.0)
    g_warning ("The ratio r/length >2 or r/length <1: %g. This will make the amplitude inaccurate!", r/l->length);

  l->so = l->dist / 2.0;
  gdouble to = sqrt((r*(l->gamma + 1.0)) / (g*gmo));
  l->lambda = 2.0*to*sqrt(g*l->depth);

  gdouble hao = 0.5*l->lambda / l->length;
  if ( hao < 1.0 )
    g_warning ("The Hammack number is < 1: %g. This will make the amplitude inaccurate\n", hao);

  gdouble sg = l->so*sint / l->depth;
  if ( sg > 0.35 )
    g_warning ("The Submergence number is > 0.35: %g. This will make the amplitude inaccurate\n", 
	       sg);

  l->eta = 0.723*l->so*(1.4662*gmo - 0.3454*pow(gmo, 2.))*pow(sint, 0.22) *
    (l->thick / l->length) * pow(l->length / l->depth, 1.25) *
    pow(dphi, 0.39) * pow(l->length / r, 0.63) * 0.1309;

  gdouble term = l->eta / (l->so*pow(sint, 1.5));
  if ( term > 0.2 )
    g_warning ("The term (eta/so*sinq^1.5) > 0.2: %g. This will make the amplitude inaccurate\n",
	       term);
}

static void gfs_init_submarine_slump_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_init_submarine_slump_class ())
   ->parent_class->write) (o, fp);
  GfsInitSubmarineSlump * l = GFS_INIT_SUBMARINE_SLUMP (o);
  fprintf (fp, " %s {\n", l->v->name);
  fprintf (fp,
	   "  x = %g y = %g alpha = %g\n"
	   "  depth = %g theta = %g\n"
	   "  length = %g width = %g\n"
	   "  thickness = %g distance = %g\n"
	   "  gamma = %g\n"
	   "}",
	   l->xo, l->yo, l->alpha,
	   l->depth, l->theta,
	   l->length, l->width,
	   l->thick, l->dist,
	   l->gamma);
}

static void init_submarine_slump (FttCell * cell, GfsInitSubmarineSlump * l)
{
  GfsSimulation * sim = gfs_object_simulation (l);

  gdouble cosa = cos(l->alpha*M_PI/180.);
  gdouble sina = sin(l->alpha*M_PI/180.);
  gdouble cost = cos(l->theta*M_PI/180.);
  gdouble sint = sin(l->theta*M_PI/180.);
  gdouble tant = tan(l->theta*M_PI/180.);

  FttVector p, o;
  gfs_cell_cm (cell, &p);
  o.x = l->xo; o.y = l->yo;
  gfs_simulation_map (sim, &o);
  gdouble L = sim->physical_params.L;
  p.x = (p.x - o.x)*L; p.y = (p.y - o.y)*L;
  FttVector q;
  q.x = -sina*p.x + cosa*p.y;
  q.y =- cosa*p.x - sina*p.y;

  gdouble g = sim->physical_params.g; /* Gravity */
  gdouble cut = 200.*l->width;

  gdouble sg = l->so*sint / l->depth;
  gdouble xg = (l->depth + l->thick / cost) / tant;
  gdouble xmin = 0.565*(xg + 0.4597*l->so*cost) - xg;
  gdouble delx = 0.5 * l->lambda;
  gdouble zmin = -l->eta * ((2.480*0.2892 - 0.7904*sg +  1.3376*pow(sg, 2.)) /
			    (0.2892 + 0.9163 * sg));
  gdouble zmax = l->eta * ((1.686*0.3498 - 0.3531*sg + 0.6466*pow(sg, 2.)) /
		       (0.3498 + 1.0257*sg));
  gdouble nmin = 1.22*1.15*zmin;
  gdouble nmax = 1.22*zmax;
  gdouble denom = 0.5*l->lambda;
  gdouble shift = 0.8*delx;
  gdouble wid = 0.5*l->lambda;

  gdouble term = 1.0 - 
    exp(-2.0906* (l->width / wid) * (1.0 + 1.0903*(l->width / wid)));

  /* Initialize tsunami shape */
  GFS_VALUE (cell, l->v) = term * 
    (nmin * exp(-pow(nmin*(q.x - xmin) / (denom*nmax), 2.)) +
     nmax * exp(-pow((q.x - xmin - shift) / (denom), 2.))) *
    pow(2. / (exp(3.*term*q.y / l->width) + exp(-3.*term*q.y / l->width)), 2.);
  if ( fabs(q.y) > cut) {
    GFS_VALUE (cell, l->v) *= exp(-pow(5.0*(fabs(q.y) - cut) / cut, 2.));
  }

  /* Initialize velocity  from estimated linear wave quantities */
  gdouble kappa = 2.0*M_PI / l->lambda;
  term = kappa*l->depth;
  gdouble omega = sqrt(g*kappa*tanh(term));

  if (GFS_VALUE (cell, l->v) > 0.) {
    gdouble utot = GFS_VALUE (cell, l->v)*g*kappa*cosh(0.469*term) /
      (omega*cosh(term));
    GFS_VALUE (cell, l->vu) = -utot*sina;
    GFS_VALUE (cell, l->vv) = utot*cosa;
  }
  else {
    GFS_VALUE (cell, l->vu) = 0.;
    GFS_VALUE (cell, l->vv) = 0.;
  }
}

static gboolean gfs_init_submarine_slump_event (GfsEvent * event,
						GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_submarine_slump_class ())
			  ->parent_class)->event) (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim),
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) init_submarine_slump,
			      event);
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_submarine_slump_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_init_submarine_slump_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_submarine_slump_write;
  klass->event = gfs_init_submarine_slump_event;
}

static GfsEventClass * gfs_init_submarine_slump_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_submarine_slump_info = {
      "GfsInitSubmarineSlump",
      sizeof (GfsInitSubmarineSlump),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_init_submarine_slump_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_submarine_slump_info);
  }

  return klass;
}

/* GfsInitSubaerialLandslide: Object */

static void gfs_init_subaerial_landslide_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_init_subaerial_landslide_class ())
   ->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a variable name");
    return;
  }
  GfsInitSubaerialLandslide * l = GFS_INIT_SUBAERIAL_LANDSLIDE (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(l->v = gfs_domain_get_or_add_variable (domain, fp->token->str, NULL))) {
    gts_file_error (fp, "`%s' is a reserved keyword", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
  gchar * vname = g_strconcat (l->v->name, "_U", NULL);
  l->vu = gfs_domain_get_or_add_variable (domain, vname, NULL);
  vname = g_strconcat (l->v->name, "_V", NULL);
  l->vv = gfs_domain_get_or_add_variable (domain, vname, NULL);

  GtsFileVariable var[] = {
    {GTS_DOUBLE, "x",         TRUE, &l->xo},
    {GTS_DOUBLE, "y",         TRUE, &l->yo},
    {GTS_DOUBLE, "alpha",     TRUE, &l->alpha},
    {GTS_DOUBLE, "volume",    TRUE, &l->vol},
    {GTS_DOUBLE, "width",     TRUE, &l->width},
    {GTS_DOUBLE, "depth",     TRUE, &l->depth},
    {GTS_DOUBLE, "ut",        TRUE, &l->ut},
    {GTS_DOUBLE, "so",        TRUE, &l->so},
    {GTS_DOUBLE, "to",        TRUE, &l->to},
    {GTS_NONE}
  };
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  gdouble g = GFS_SIMULATION(domain)->physical_params.g;
  gdouble fr = l->ut / sqrt( g * l->depth);
  if (fr > 4.0)
    g_warning ("The Froude number is > 4. This will make the amplitude inaccurate!\n");
  if (fr < 1.0)
    g_warning ("The Froude number is < 1. This will make the amplitude inaccurate!\n");

  l->lambda = 0.27*l->to*sqrt(g*l->depth);

  gdouble etal = 1.32*l->depth*
    pow(l->vol*2.0*l->ut / (M_PI*l->width*l->so*l->depth*sqrt(l->depth*g)), 
	0.68);
  gdouble etat = 1.32*l->depth*pow(l->vol /
				   (l->width*l->to*l->depth*sqrt(l->depth*g)),
				   0.68);
  
  /* Choose the smallest amplitude for two measures of time. */
  if (etal < etat)
    l->eta = etal;
  else
    l->eta = etat;

  if ( l->eta > 0.86*l->depth) {
    g_warning ("The wave amplitude is excessively large. The amplitude will be decreased to 0.86*depth!\n");
      l->eta = 0.86*l->depth;
  }

  gdouble test = 2.0*fabs(etal - etat)/(etal + etat);
  if ( test > 4.0 )
    g_warning ("The two wave amplitudes differ more than 40 percents. This will make the results questionable\n");

  test = 2.0*l->to*l->ut / (l->so*M_PI);
  if ( test < 0.4 || test > 2.5 )
    g_warning ("Runout length and time differ. This will make the results questionable!\n");

  test = l->to / (4.5*sqrt(10.0*sqrt(l->vol/l->width) / g));
  if ( test < 0.3 || test > 3.3 )
    g_warning ("Runout length and time differ. This will make the results questionable!\n");

  test = 3.4*pow(l->vol, 1./3.)/l->width;
  if ( test < 0.25 || test > 4.0 )
    g_warning ("Landslide width may be unusual. This will make the results questionable!\n");
}

static void gfs_init_subaerial_landslide_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_init_subaerial_landslide_class ())
   ->parent_class->write) (o, fp);
  GfsInitSubaerialLandslide * l = GFS_INIT_SUBAERIAL_LANDSLIDE (o);
  fprintf (fp, " %s {\n", l->v->name);
  fprintf (fp,
	   "  x = %g y = %g alpha = %g\n"
	   "  volume = %g width = %g\n"
	   "  depth = %g ut = %g\n"
	   "  so = %g to = %g\n"
	   "}",
	   l->xo, l->yo, l->alpha,
	   l->vol, l->width, l->depth,
	   l->ut, l->so, l->to);
}

static void init_subaerial_landslide (FttCell * cell, GfsInitSubaerialLandslide * l)
{
  GfsSimulation * sim = gfs_object_simulation (l);

  gdouble cosa = cos(l->alpha*M_PI/180.);
  gdouble sina = sin(l->alpha*M_PI/180.);

  FttVector p, o;
  gfs_cell_cm (cell, &p);
  o.x = l->xo; o.y = l->yo;
  gfs_simulation_map (sim, &o);
  gdouble L = sim->physical_params.L;
  p.x = (p.x - o.x)*L; p.y = (p.y - o.y)*L;
  FttVector q;
  q.x = -sina*p.x + cosa*p.y;
  q.y =- cosa*p.x - sina*p.y;

  gdouble g = sim->physical_params.g; /* Gravity */

  gdouble xmin = l->so + l->lambda;

  /* Evaluate tsunami shape in local x,y coordinates. */
  gdouble term = 4.0*l->eta*(l->width + l->lambda)*l->lambda / l->vol;
  if (((l->width + l->lambda) / term) < l->lambda) {
    term = (l->width + l->lambda) / l->lambda;
    l->eta = term*l->vol / (4.0*(l->width + l->lambda)*l->lambda);
  }
  
  GFS_VALUE (cell, l->v) = l->eta*
    pow(2.0 / ( exp(term*q.y / (l->width + l->lambda)) +
		exp(-term*q.y / (l->width + l->lambda)) ), 2.) *
    pow(2.0 / ( exp((q.x - xmin) / l->lambda) +
		exp(-(q.x - xmin) / l->lambda)), 2.);

  /* Initialize tsunami velocities */
  term = l->eta / l->depth;
  if (GFS_VALUE (cell, l->v) > 0.) {
    gdouble utot = sqrt (g*l->depth) * (1.0 + term / 2.0) *
      ((1.0 + 0.17006*term) * GFS_VALUE (cell, l->v) / l->depth
       - 1.25509*pow(GFS_VALUE (cell, l->v) / l->depth, 2.));
    
    GFS_VALUE (cell, l->vu) = -utot*sina;
    GFS_VALUE (cell, l->vv) = utot*cosa;
  }
  else {
    GFS_VALUE (cell, l->vu) = 0.;
    GFS_VALUE (cell, l->vv) = 0.;
  }
}

static gboolean gfs_init_subaerial_landslide_event (GfsEvent * event,
						    GfsSimulation * sim)
{
  if ((*GFS_EVENT_CLASS(GTS_OBJECT_CLASS (gfs_init_subaerial_landslide_class ())
			->parent_class)->event)
      (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim),
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) init_subaerial_landslide,
			      event);
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_subaerial_landslide_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_init_subaerial_landslide_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_subaerial_landslide_write;
  klass->event = gfs_init_subaerial_landslide_event;
}

static GfsEventClass * gfs_init_subaerial_landslide_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_subaerial_landslide_info = {
      "GfsInitSubaerialLandslide",
      sizeof (GfsInitSubaerialLandslide),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_init_subaerial_landslide_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_subaerial_landslide_info);
  }

  return klass;
}

/* GfsInitPyroclastic: Object */

static void gfs_init_pyroclastic_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_init_pyroclastic_class ())
   ->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a variable name");
    return;
  }
  GfsInitPyroclastic * l = GFS_INIT_PYROCLASTIC (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(l->v = gfs_domain_get_or_add_variable (domain, fp->token->str, NULL))) {
    gts_file_error (fp, "`%s' is a reserved keyword", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
  gchar * vname = g_strconcat (l->v->name, "_U", NULL);
  l->vu = gfs_domain_get_or_add_variable (domain, vname, NULL);
  vname = g_strconcat (l->v->name, "_V", NULL);
  l->vv = gfs_domain_get_or_add_variable (domain, vname, NULL);

  GtsFileVariable var[] = {
    {GTS_DOUBLE, "x",         TRUE, &l->xo},
    {GTS_DOUBLE, "y",         TRUE, &l->yo},
    {GTS_DOUBLE, "alpha",     TRUE, &l->alpha},
    {GTS_DOUBLE, "volume",    TRUE, &l->vol},
    {GTS_DOUBLE, "width",     TRUE, &l->width},
    {GTS_DOUBLE, "depth",     TRUE, &l->depth},
    {GTS_DOUBLE, "ut",        TRUE, &l->ut},
    {GTS_DOUBLE, "so",        TRUE, &l->so},
    {GTS_DOUBLE, "to",        TRUE, &l->to},
    {GTS_NONE}
  };
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  gdouble g = GFS_SIMULATION(domain)->physical_params.g; /* Gravity */

  gdouble fr = l->ut / sqrt(g*l->depth);
  if (fr > 4.0)
    g_warning ("The Froude number is > 4. This will make the amplitude inaccurate!\n");
  if (fr < 1.0)
    g_warning ("The Froude number is < 1. This will make the amplitude inaccurate!\n");

  l->lambda = 0.27*l->to*sqrt(g*l->depth);
  gdouble etal = 1.32*l->depth *
    pow(l->vol*2.0*l->ut / (M_PI*l->width*l->so*l->depth*sqrt(l->depth*g)),
	0.68);
  gdouble etat = 1.32*l->depth * 
    pow(l->vol / (l->width*l->to*l->depth*sqrt(l->depth*g)), 0.68);

  gdouble test = 2.*fabs(etal - etat) / (etal + etat);
  if (test > 0.4)
    g_warning ("The two wave amplitudes differ more than 40 percent. This will make the results questionable!\n");

  /* Choose the smallest amplitude for two measures of time */
  if (etal < etat)
    l->eta = etal;
  else
    l->eta = etat;

  test = 0.86*l->depth;
  if (l->eta > test) {
    g_warning ("The wave amplitude is excessively large. The amplitude will be decreased to 0.86*depth\n");
    l->eta = test;
  }

  test = 2.0*l->to*l->ut / (l->so*M_PI);
  if ( test < 0.4 || test > 2.5)
    g_warning ("Runout length and time differ: %g. This will make the results questionable!\n", 
	       test);

  test = l->to / (4.5*sqrt(10.0*sqrt(l->vol / l->width) / g));
  if ( test < 0.3 || test > 3.3)
    g_warning ("Runout length and time differ: %g. This will make the results questionable!\n", 
	       test);

  test = 30.0*pow(l->vol, 1.0 / 3.0) / l->width;
  if ( test < 0.25 || test > 4.0)
    g_warning ("Flow width may be unusual: %g. This will make the results questionable!\n", test);
}

static void gfs_init_pyroclastic_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_init_pyroclastic_class ())->parent_class->write)
    (o, fp);
  GfsInitPyroclastic * l = GFS_INIT_PYROCLASTIC (o);
  fprintf (fp, " %s {\n", l->v->name);
  fprintf (fp, 
	   "  x = %g y = %g alpha = %g\n"
	   "  volume = %g width = %g\n"
	   "  depth = %g ut = %g\n"
	   "  so = %g to = %g\n"
	   "}",
	   l->xo, l->yo, l->alpha,
	   l->vol, l->width, l->depth,
	   l->ut, l->so, l->to);
}

static void init_pyroclastic (FttCell * cell, GfsInitPyroclastic * l)
{
  GfsSimulation * sim = gfs_object_simulation (l);

  gdouble cosa = cos(l->alpha*M_PI/180.);
  gdouble sina = sin(l->alpha*M_PI/180.);

  FttVector p, o;
  gfs_cell_cm (cell, &p);
  o.x = l->xo; o.y = l->yo;
  gfs_simulation_map (sim, &o);
  gdouble L = sim->physical_params.L;
  p.x = (p.x - o.x)*L; p.y = (p.y - o.y)*L;
  FttVector q;
  q.x = -sina*p.x + cosa*p.y;
  q.y =- cosa*p.x - sina*p.y;

  gdouble g = sim->physical_params.g; /* Gravity */

  gdouble xmin = l->so + l->lambda;

  /* Evaluate tsunami shape in local x,y coordinates. */
  gdouble term = 4.0*l->eta*(l->width + l->lambda)*l->lambda / l->vol;
  if (((l->width + l->lambda) / term) < l->lambda) {
    term = (l->width + l->lambda) / l->lambda;
    l->eta = term*l->vol / (4.0*(l->width + l->lambda)*l->lambda);
  }

  GFS_VALUE (cell, l->v) = l->eta *
    pow(2.0 / (exp(term*q.y / (l->width + l->lambda)) +
	       exp(-term*q.y / (l->width + l->lambda))), 2.) *
    pow(2.0 / (exp((q.x - xmin) / l->lambda) +
	       exp(-(q.x - xmin) / l->lambda)), 2.);
    
  /* Initialize tsunami velocities */
  term = l->eta / l->depth;
  
  if (GFS_VALUE (cell, l->v) > 0.) {
    gdouble utot = sqrt (g*l->depth)*(1.0 + term / 2.0) * 
      ((1.0 + 0.17006*term)*GFS_VALUE (cell, l->v) / l->depth
       - 1.25509*pow(GFS_VALUE (cell, l->v) / l->depth, 2.));
    
    GFS_VALUE (cell, l->vu) = -utot*sina;
    GFS_VALUE (cell, l->vv) = utot*cosa;
  }
  else {
    GFS_VALUE (cell, l->vu) = 0.;
    GFS_VALUE (cell, l->vv) = 0.;
  }
}

static gboolean gfs_init_pyroclastic_event (GfsEvent * event, 
					    GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_pyroclastic_class ())
			  ->parent_class)->event) (event, sim)) {
    gfs_domain_cell_traverse (GFS_DOMAIN (sim), 
			      FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) init_pyroclastic, event);
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_pyroclastic_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_init_pyroclastic_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_pyroclastic_write;
  klass->event = gfs_init_pyroclastic_event;
}

static GfsEventClass * gfs_init_pyroclastic_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_pyroclastic_info = {
      "GfsInitPyroclastic",
      sizeof (GfsInitPyroclastic),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_init_pyroclastic_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_pyroclastic_info);
  }

  return klass;
}

/* GfsOutputInterfaceGrid: Header */

typedef struct _GfsOutputInterfaceGrid         GfsOutputInterfaceGrid;

struct _GfsOutputInterfaceGrid {
  /*< private >*/
  GfsOutput parent;
  
  /*< public >*/
  GfsVariable * v, * t;

  gdouble cost, sint;
  gdouble cosa, sina;
  gdouble theta, alpha;
  gdouble x, y;
  gdouble sx, sy, sz;
};

#define GFS_OUTPUT_INTERFACE_GRID(obj)      GTS_OBJECT_CAST (obj,\
					    GfsOutputInterfaceGrid,	\
					    gfs_output_interface_grid_class ())
#define GFS_IS_OUTPUT_INTERFACE_GRID(obj)   (gts_object_is_from_class (obj,\
					    gfs_output_interface_grid_class ()))

static GfsOutputClass * gfs_output_interface_grid_class  (void);


static void gfs_output_interface_grid_read (GtsObject ** o, GtsFile * fp)
{
  GfsOutputInterfaceGrid * output = GFS_OUTPUT_INTERFACE_GRID (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));

  if (GTS_OBJECT_CLASS (gfs_output_interface_grid_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_output_interface_grid_class ())
     ->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a variable (T)");
    return;
  }
  if ((output->t = gfs_variable_from_name (domain->variables, fp->token->str))
      == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_DOUBLE, "theta", TRUE, &output->theta},
      {GTS_DOUBLE, "x",     TRUE, &output->x},
      {GTS_DOUBLE, "y",     TRUE, &output->y},
      {GTS_DOUBLE, "alpha", TRUE, &output->alpha},
      {GTS_DOUBLE, "sx",    TRUE, &output->sx},
      {GTS_DOUBLE, "sy",    TRUE, &output->sy},
      {GTS_DOUBLE, "sz",    TRUE, &output->sz},
      {GTS_NONE}
    }; 
    gts_file_assign_variables (fp, var);
    if (fp->type == GTS_ERROR)
      return;
    output->cost = cos(output->theta*M_PI/180.);
    output->sint = sin(output->theta*M_PI/180.);
    output->cosa = cos(output->alpha*M_PI/180.);
    output->sina = sin(output->alpha*M_PI/180.);
  }

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a variable (V)");
    return;
  }
  if ((output->v = gfs_variable_from_name (domain->variables, fp->token->str))
      == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
}

static void gfs_output_interface_grid_write (GtsObject * o, FILE * fp)
{
  GfsOutputInterfaceGrid * output = GFS_OUTPUT_INTERFACE_GRID (o);

  if (GTS_OBJECT_CLASS (gfs_output_interface_grid_class ())
      ->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_output_interface_grid_class ())
     ->parent_class->write) (o, fp);

  fprintf (fp, " %s ",output->t->name);
  fputs(" {\n", fp);
  fprintf(fp, "   theta = %g\n", output->theta);
  fprintf(fp, "   x = %g\n", output->x);
  fprintf(fp, "   y = %g\n", output->y);
  fprintf(fp, "   alpha = %g\n", output->alpha);
  fprintf(fp, "   sx = %g\n", output->sx);
  fprintf(fp, "   sy = %g\n", output->sy);
  fprintf(fp, "   sz = %g\n", output->sz);
  fputs(" }", fp);
  fprintf (fp, " %s ",output->v->name);
}

#define NODATA 0.

typedef struct {
  double xmin, xmax;
  double ymin, ymax;
  double size;
  gint xn, yn;
  GfsOutputInterfaceGrid * output;
  gdouble * buf, ** data;
} GridData;

static gboolean interface_condition (FttCell * cell, gpointer data)
{
  GfsOutputInterfaceGrid * out = (GfsOutputInterfaceGrid *) data;

  if (GFS_VALUE (cell, out->t) < 1. && GFS_VALUE (cell, out->t) > 0)
    return TRUE;
  return FALSE;
}

static void print_interface (FttCell * cell, gpointer * data)
{
  GridData * d = (GridData *) data;
  GfsOutputInterfaceGrid * out = d->output;
  
  if (FTT_CELL_IS_LEAF (cell) && GFS_VALUE (cell, out->t) < 1. &&
      GFS_VALUE (cell, out->t) > 0) {
    GfsVariableTracerVOF * tv = GFS_VARIABLE_TRACER_VOF (out->t);

    double nx = GFS_VALUE(cell, tv->m[0]);
    double ny = GFS_VALUE(cell, tv->m[1]);
    double nz = GFS_VALUE(cell, tv->m[2]);
    double alpha = GFS_VALUE(cell, tv->alpha);
    
    FttVector pos;
    ftt_cell_pos (cell, &pos);
    double h = ftt_cell_size(cell);
    
    gint i, j;
    FttVector p;
    for ( i = 0; i < d->xn; i ++)
      for ( j = 0; j < d->yn; j ++) {
	double lon = d->xmin + i*d->size;
	double lat = d->ymin + j*d->size;
	double X = -lon*out->sina - lat*out->cosa;
	double Y = lon*out->cosa - lat*out->sina;
	p.x = X*out->cost;
	p.z = Y;

	if (p.x <= pos.x + h/2 && p.x > pos.x - h/2)
	  if (p.z <= pos.z + h/2 && p.z > pos.z - h/2) {
	    p.y = pos.y - h/2. + (alpha - nx*(p.x-pos.x+h/2.)/h -
				  nz*(p.z-pos.z+h/2.)/h) / ny*h;
	    d->data[i][j] = gfs_interpolate (cell, p, out->v);
	  }
      }
  }
}

static void extent (FttCell * cell, gpointer * data)
{
  GridData * d = (GridData *) data;
  GfsOutputInterfaceGrid * out = d->output;

  if (FTT_CELL_IS_LEAF (cell) && GFS_VALUE (cell, out->t) < 1. &&
      GFS_VALUE (cell, out->t) > 0) {
    FttVector pos;
    ftt_cell_pos (cell, &pos);
    
    /* Projection */
    double X = out->cost*pos.x + out->sint*pos.y;
    double Y = pos.z;

    /* Rotation */
    double x = -out->sina*X + out->cosa*Y;
    double y = -out->sina*Y - out->cosa*X;

    if (d->xmin > x )
      d->xmin = x;
    if (d->ymin > y)
      d->ymin = y;
    if (d->xmax < x)
      d->xmax = x;
    if (d->ymax < y)
      d->ymax = y;
  }
}

static gboolean gfs_output_interface_grid_event (GfsEvent * event,
						 GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_output_interface_grid_class ())
			  ->parent_class)->event) (event, sim)) {
     GfsOutputInterfaceGrid * out = GFS_OUTPUT_INTERFACE_GRID (event);
     GfsDomain * domain = GFS_DOMAIN (sim);

     /* New GridData structure */
     GridData * data = g_malloc0 (sizeof (GridData));
     data->output = out;
     data->size = 1./pow(2.,gfs_domain_depth (domain));
     data->xmin = data->ymin = G_MAXDOUBLE;
     data->xmax = data->ymax = -G_MAXDOUBLE;

     /* Get the size of the interface */
     gfs_domain_cell_traverse_condition (GFS_DOMAIN(sim), FTT_PRE_ORDER,
					 FTT_TRAVERSE_LEAFS, -1,
					 (FttCellTraverseFunc) extent, data,
  					 interface_condition, out);

     /* Number of points of the grid */
     data->xn = (data->xmax - data->xmin)/data->size+1;
     data->yn = (data->ymax - data->ymin)/data->size+1;
     
     /* Memory allocation for the grid */
     gint i,j;
     data->buf = g_malloc (sizeof (gdouble)*data->xn*data->yn);
     for ( i = 0; i < data->xn*data->yn; i++)
       data->buf[i] = NODATA;     
     data->data = g_malloc (sizeof (gdouble *)*data->xn);     
     for ( i = 0; i < data->xn; i++)
       data->data[i] = &data->buf[i*data->yn];

      
     /* Data collection */
     gfs_domain_cell_traverse_condition (GFS_DOMAIN(sim), FTT_PRE_ORDER,
					 FTT_TRAVERSE_LEAFS, -1,
					 (FttCellTraverseFunc) print_interface,
					 data, interface_condition, out);

     /* Print the data in lon/lat coordinate system in cgd format*/
     double lat_meter_to_degrees = 180./(M_PI*6371220.);
     double lon_meter_to_degrees = 180./(M_PI*6371220.*cos(out->y*M_PI/180.));
     
     fprintf(GFS_OUTPUT(out)->file->fp,"2 x y\n");
     fprintf(GFS_OUTPUT(out)->file->fp,"%i %i\n",data->xn,data->yn);
     for (i=0; i < data->xn; i++) {
       double x0 = (data->xmin + i*data->size)*sim->physical_params.L*
	 lon_meter_to_degrees*out->sx + out->x;
       fprintf(GFS_OUTPUT(out)->file->fp,"%f ", x0);
     }
     fprintf(GFS_OUTPUT(out)->file->fp,"\n");
     
     for ( i = 0; i < data->yn; i++) {
       double y0 = (data->ymin + i*data->size)*sim->physical_params.L*
	 lat_meter_to_degrees*out->sy + out->y;
       fprintf(GFS_OUTPUT(out)->file->fp,"%f ", y0);
     }
     fprintf(GFS_OUTPUT(out)->file->fp,"\n");
     
     for ( i = 0; i < data->xn; i++) {
       for ( j = 0; j < data->yn; j++) {
     	 if (data->data[i][j] == NODATA)
     	   fprintf(GFS_OUTPUT(out)->file->fp,"%f ", 0.);
     	 else
     	   fprintf(GFS_OUTPUT(out)->file->fp,"%f ", data->data[i][j]*
		   sim->physical_params.L*out->sz);
       }
       fprintf(GFS_OUTPUT(out)->file->fp,"\n");
     }

     /* Frees GridData */
     g_free (data->data);
     g_free (data->buf);
     g_free (data);
     return TRUE;
  }
  return FALSE;
}

static void gfs_output_interface_grid_class_init (GfsOutputClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_output_interface_grid_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_output_interface_grid_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_output_interface_grid_write;
}

static void gfs_output_interface_grid_init (GfsOutputInterfaceGrid * object)
{
  object->x = 0.;
  object->y = 0.;
  object->alpha = 0.;
  object->sx = object->sy = object->sz = 1.;
}

GfsOutputClass * gfs_output_interface_grid_class (void)
{
  static GfsOutputClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_output_interface_grid_info = {
      "GfsOutputInterfaceGrid",
      sizeof (GfsOutputInterfaceGrid),
      sizeof (GfsOutputClass),
      (GtsObjectClassInitFunc) gfs_output_interface_grid_class_init,
      (GtsObjectInitFunc) gfs_output_interface_grid_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_output_class ()),
				  &gfs_output_interface_grid_info);
  }

  return klass;
}

/* Initialize module */

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "topics";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  gfs_init_submarine_landslide_class ();
  gfs_init_submarine_slump_class ();
  gfs_init_subaerial_landslide_class ();
  gfs_init_pyroclastic_class ();
  gfs_output_interface_grid_class ();
  return NULL;
}
