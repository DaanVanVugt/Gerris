/* Gerris - The GNU Flow Solver
 * Copyright (C) 2009-2012 National Institute of Water and Atmospheric Research
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

#include "event.h"
#include "solid.h"

/* Implementation of the formulae of Okada, 1985, "Surface deformation
   due to shear and tensile faults in a half-space", Bulletin of the
   Seismological Society of America, 75:4, 1135-1154, */

/* formulae (25)-(30) */
static void rectangular_source (const double U[3], double cosd, double sind,
				double mulambda, double d,
				double psi, double eta, double q,
				double u[3])
{
  double R = sqrt (psi*psi + eta*eta + q*q);
  double X = sqrt (psi*psi + q*q);
  double dtilde = eta*sind - q*cosd;
  double ytilde = eta*cosd + q*sind;
  double atanp = fabs (q) > 1e-6 ? atan (psi*eta/(q*R)) : 0.;

  mulambda = mulambda/(1. + mulambda);
  double logReta = R + eta > 1e-6 ? log (R + eta) : - log (R - eta);
  double Reta = fabs (R + eta) > 1e-6 ? R + eta : 1e30;
  double I1, I2, I3, I4, I5;
  if (fabs (cosd) > 1e-6) {
    /* formula (28) */
    I5 = fabs (psi) < 1e-6 ? 0. :
      mulambda*2./cosd*atan ((eta*(X + q*cosd) + X*(R + X)*sind)/(psi*(R + X)*cosd));
    I4 = mulambda/cosd*(log (R + dtilde) - sind*logReta);
    I3 = mulambda*(1./cosd*ytilde/(R + dtilde) - logReta) + sind/cosd*I4;
    I2 = mulambda*(- logReta) - I3;
    I1 = mulambda*(-1./cosd*psi/(R + dtilde)) - sind/cosd*I5;
  }
  else {
    /* formula (29) */
    double R1 = R + dtilde;
    I1 = - mulambda/2.*psi*q/(R1*R1);
    I3 = mulambda/2.*(eta/R1 + ytilde*q/(R1*R1) - logReta);
    I2 = mulambda*(- logReta) - I3;
    I4 = - mulambda*q/R1;
    I5 = - mulambda*psi*sind/R1;
  }
    
  /* strike-slip, formula (25) */  
  if (U[0] != 0.) {
    double U1pi = U[0]/(2.*M_PI);
    u[0] -= U1pi*(psi*q/(R*Reta) + atanp + I1*sind);
    u[1] -= U1pi*(ytilde*q/(R*Reta) + q*cosd/Reta + I2*sind);
    u[2] -= U1pi*(dtilde*q/(R*Reta) + q*sind/Reta + I4*sind);
  }

  /* dip-slip, formula (26) */  
  if (U[1] != 0.) {
    double U2pi = U[1]/(2.*M_PI);
    u[0] -= U2pi*(q/R - I3*sind*cosd);
    u[1] -= U2pi*(ytilde*q/(R*(R + psi)) + cosd*atanp - I1*sind*cosd);
    u[2] -= U2pi*(dtilde*q/(R*(R + psi)) + sind*atanp - I5*sind*cosd);
  }

  /* tensile, formula (27) */  
  if (U[2] != 0.) {
    double U3pi = U[2]/(2.*M_PI);
    u[0] += U3pi*(q*q/(R*Reta) - I3*sind*sind);
    u[1] += U3pi*(-dtilde*q/(R*(R + psi)) - sind*(psi*q/(R*Reta) - atanp) - I1*sind*sind);
    u[2] += U3pi*(ytilde*q/(R*(R + psi)) + cosd*(psi*q/(R*Reta) - atanp) - I5*sind*sind);
  }
}

/* formula (24) */
static void okada_rectangular_source (const double U[3], 
				      double L, double W, double d, double delta, double mulambda,
				      double x, double y,
				      double u[3])
{
  double cosd = cos (delta), sind = sin (delta);
  double p = y*cosd + d*sind;
  double q = y*sind - d*cosd;

  u[0] = u[1] = u[2] = 0.;
  rectangular_source (U, cosd, sind, mulambda, d,
		      x, p, q,
		      u);
  rectangular_source (U, cosd, sind, mulambda, d,
		      x - L, p - W, q,
		      u);

  double u1[3] = {0., 0., 0.};
  rectangular_source (U, cosd, sind, mulambda, d,
		      x, p - W, q,
		      u1);
  rectangular_source (U, cosd, sind, mulambda, d,
		      x - L, p, q,
		      u1);
  u[0] -= u1[0];
  u[1] -= u1[1];
  u[2] -= u1[2];
}

/* GfsInitOkada: Header */

typedef struct _GfsInitOkada         GfsInitOkada;

struct _GfsInitOkada {
  /*< private >*/
  GfsGenericInit parent;
  double sina, cosa;

  /*< public >*/
  GfsVariable * v;
  gdouble x, y, depth;
  gdouble strike, dip;
  gdouble mu, lambda;
  gdouble length, width, U[3];
  gdouble R;
};

#define GFS_INIT_OKADA(obj)            GTS_OBJECT_CAST (obj,\
					         GfsInitOkada,\
					         gfs_init_okada_class ())
#define GFS_IS_INIT_OKADA(obj)         (gts_object_is_from_class (obj,	\
						 gfs_init_okada_class ()))

static GfsEventClass * gfs_init_okada_class  (void);

/* GfsInitOkada: Object */

static void gfs_init_okada_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_init_okada_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a variable name");
    return;
  }
  GfsInitOkada * okada = GFS_INIT_OKADA (*o);
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
  if (!(okada->v = gfs_domain_get_or_add_variable (domain, fp->token->str, NULL))) {
    gts_file_error (fp, "`%s' is a reserved keyword", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  gdouble U = 0., rake = 90.;
  GtsFileVariable var[] = {
    {GTS_DOUBLE, "x",      TRUE, &okada->x},      /* 0 */
    {GTS_DOUBLE, "y",      TRUE, &okada->y},      /* 1 */
    {GTS_DOUBLE, "depth",  TRUE, &okada->depth},  /* 2 */
    {GTS_DOUBLE, "strike", TRUE, &okada->strike}, /* 3 */
    {GTS_DOUBLE, "dip",    TRUE, &okada->dip},    /* 4 */
    {GTS_DOUBLE, "rake",   TRUE, &rake},          /* 5 */
    {GTS_DOUBLE, "mu",     TRUE, &okada->mu},     /* 6 */
    {GTS_DOUBLE, "lambda", TRUE, &okada->lambda}, /* 7 */
    {GTS_DOUBLE, "length", TRUE, &okada->length}, /* 8 */
    {GTS_DOUBLE, "width",  TRUE, &okada->width},  /* 9 */
    {GTS_DOUBLE, "U1",     TRUE, &okada->U[0]},   /* 10 */
    {GTS_DOUBLE, "U2",     TRUE, &okada->U[1]},   /* 11 */
    {GTS_DOUBLE, "U3",     TRUE, &okada->U[2]},   /* 12 */ 
    {GTS_DOUBLE, "U",      TRUE, &U},             /* 13 */ 
    {GTS_DOUBLE, "R",      TRUE, &okada->R},      /* 14 */ 
    {GTS_NONE}
  };
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR)
    return;

  if (var[5].set) {
    if (var[10].set || var[11].set) {
      gts_file_error (fp, "set rake and U or U1 and U2 (not both)");
      return;
    }
    okada->U[0] = U*cos (rake*M_PI/180.);
    okada->U[1] = U*sin (rake*M_PI/180.);
  }

  okada->sina = sin ((90. - okada->strike)*M_PI/180.);
  okada->cosa = cos ((90. - okada->strike)*M_PI/180.);
}

static void gfs_init_okada_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_init_okada_class ())->parent_class->write) (o, fp);
  GfsInitOkada * okada = GFS_INIT_OKADA (o);
  fprintf (fp, " %s {\n", okada->v->name);
  fprintf (fp, 
	   "  x = %g y = %g depth = %g\n"
	   "  strike = %g dip = %g\n"
	   "  length = %g width = %g\n"
	   "  U1 = %g U2 = %g U3 = %g\n"
	   "  mu = %g lambda = %g\n"
	   "}",
	   okada->x, okada->y, okada->depth,
	   okada->strike, okada->dip,
	   okada->length, okada->width, 
	   okada->U[0],	okada->U[1], okada->U[2],
	   okada->mu, okada->lambda);
}

static double delta (double theta1, double theta2)
{
  double d = theta1 - theta2;
  if (d > 180.) d -= 360.;
  if (d < -180.) d += 360.;
  return d;
}

static void init_okada (FttCell * cell, GfsInitOkada * okada)
{
  FttVector p;
  gfs_cell_cm (cell, &p);
  GfsSimulation * sim = gfs_object_simulation (okada);
  gfs_simulation_map_inverse (sim, &p);
  p.x = okada->R*cos(p.y*M_PI/180.)*delta (p.x, okada->x)*M_PI/180.;
  p.y = okada->R*delta (p.y, okada->y)*M_PI/180.;
  FttVector q;
  q.x =   okada->cosa*p.x + okada->sina*p.y;
  q.y = - okada->sina*p.x + okada->cosa*p.y;
  double u[3];
  double sind = sin (okada->dip*M_PI/180.);
  /* depth of the bottom edge */
  double depth = sind > 0. ? okada->depth + okada->width*sind : okada->depth;
  /* shift origin to the centroid */
  q.x += okada->length/2.;
  q.y += okada->width/2.*cos (okada->dip*M_PI/180.);
  okada_rectangular_source (okada->U, okada->length, okada->width, depth, 
			    okada->dip*M_PI/180.,
			    okada->mu/okada->lambda,
			    q.x, q.y,
			    u);
  GFS_VALUE (cell, okada->v) += u[2];
}

static gboolean gfs_init_okada_event (GfsEvent * event, 
				      GfsSimulation * sim)
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_init_okada_class ())->parent_class)->event) 
      (event, sim)) {
    gfs_domain_traverse_leaves (GFS_DOMAIN (sim), (FttCellTraverseFunc) init_okada, event);
    return TRUE;
  }
  return FALSE;
}

static void gfs_init_okada_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = gfs_init_okada_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_init_okada_write;
  klass->event = gfs_init_okada_event;
}

static void gfs_init_okada_init (GfsInitOkada * okada)
{
  okada->mu = okada->lambda = 1.;
  okada->R = 6371220.; /* Earth radius (metres) */
}

static GfsEventClass * gfs_init_okada_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_init_okada_info = {
      "GfsInitOkada",
      sizeof (GfsInitOkada),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_init_okada_class_init,
      (GtsObjectInitFunc) gfs_init_okada_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_generic_init_class ()),
				  &gfs_init_okada_info);
  }

  return klass;
}

/* Initialize module */

/* only define gfs_module_name for "official" modules (i.e. those installed in
   GFS_MODULES_DIR) */
const gchar gfs_module_name[] = "okada";
const gchar * g_module_check_init (void);

const gchar * g_module_check_init (void)
{
  gfs_init_okada_class ();
  return NULL;
}
