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

#include <stdlib.h>

#include "particulatecommon.h"
#include "source.h"

/* Forces on the Particle */

static FttVector subs_fttvectors (FttVector *a, FttVector *b)
{
  FttVector result;
  FttComponent c;
  for(c = 0; c< FTT_DIMENSION; c++)    
    (&result.x)[c]  = (&a->x)[c] - (&b->x)[c];  
  return result;
}

/* Same as in source.c used here to obtained viscosity */
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

/* Similar to gfs_vorticity which returns norm of the vorticity */
static void vorticity_vector (FttCell *cell, GfsVariable **v, 
			      FttVector *vort) 
{
  gdouble size;

  if (cell == NULL) return;
  if (v == NULL) return;

  size = ftt_cell_size (cell);
#if FTT_2D
  vort->x = 0.;
  vort->y = 0.;
  vort->z = (gfs_center_gradient (cell, FTT_X, v[1]->i) -
	     gfs_center_gradient (cell, FTT_Y, v[0]->i))/size;
#else  /* FTT_3D */
  vort->x = (gfs_center_gradient (cell, FTT_Y, v[2]->i) -
	     gfs_center_gradient (cell, FTT_Z, v[1]->i))/size;
  vort->y = (gfs_center_gradient (cell, FTT_Z, v[0]->i) -
	     gfs_center_gradient (cell, FTT_X, v[2]->i))/size;
  vort->z = (gfs_center_gradient (cell, FTT_X, v[1]->i) -
	     gfs_center_gradient (cell, FTT_Y, v[0]->i))/size;
#endif
}

/* GfsForceCoeff: object */

static void gfs_force_coeff_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_force_coeff_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_force_coeff_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != '\n') {
    GfsForceCoeff * force = FORCE_COEFF (*o);
    force->coefficient = gfs_function_new (gfs_function_class (), 0.);
    gfs_function_read (force->coefficient, gfs_object_simulation (*o), fp);
    GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (*o));
    
    /* fixme: "Rep", "Urelp" etc... should be derived variables not
       straight variables (i.e. there is no need to allocate memory
       for these as they are only used temporarily to compute the
       coefficient) */
    force->re_p = gfs_domain_get_or_add_variable (domain, "Rep", 
						  "Particle Reynolds number");  
    force->u_rel = gfs_domain_get_or_add_variable (domain, "Urelp", 
						   "Particle x - relative velocity");
    force->v_rel = gfs_domain_get_or_add_variable (domain, "Vrelp", 
						   "Particle y - relative velocity");
#if !FTT_2D
    force->w_rel = gfs_domain_get_or_add_variable (domain, "Wrelp", 
						   "Particle z - relative velocity");
#endif
    force->pdia = gfs_domain_get_or_add_variable (domain, "Pdia", 
						  "Particle radii");
  }
}

static void gfs_force_coeff_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_force_coeff_class ())->parent_class->write) (o, fp);
  GfsForceCoeff * force = FORCE_COEFF (o);
  if (force->coefficient)
    gfs_function_write (force->coefficient, fp);
}

static void gfs_force_coeff_destroy (GtsObject * o)
{
  if (FORCE_COEFF (o)->coefficient)
    gts_object_destroy (GTS_OBJECT (FORCE_COEFF (o)->coefficient));

  (* GTS_OBJECT_CLASS (gfs_force_coeff_class ())->parent_class->destroy) (o);
}

static void gfs_force_coeff_class_init (GtsObjectClass * klass)
{
  klass->read = gfs_force_coeff_read;
  klass->write = gfs_force_coeff_write;
  klass->destroy = gfs_force_coeff_destroy;
}
 
GtsSListContaineeClass * gfs_force_coeff_class (void)
{
  static GtsSListContaineeClass * klass = NULL;
  
  if (klass == NULL) {
    GtsObjectClassInfo gfs_force_coeff_info = {
      "GfsForceCoeff",
      sizeof (GfsForceCoeff),
      sizeof (GtsSListContaineeClass),
      (GtsObjectClassInitFunc) gfs_force_coeff_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particle_force_class ()),
				  &gfs_force_coeff_info);
  }
  return klass;
}

/* GfsForceLift: object */

static FttVector compute_lift_force (GfsParticle * p, GfsParticleForce * liftforce)
{
  GfsParticulate * particulate = GFS_PARTICULATE (p);
  GfsForceCoeff * coeff = FORCE_COEFF (liftforce);

  GfsSimulation * sim = gfs_object_simulation (particulate);
  GfsDomain * domain = GFS_DOMAIN (sim);
  
  FttVector force;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] = 0;

  FttCell * cell = gfs_domain_locate (domain, p->pos, -1, NULL);
  if (cell == NULL) return force;

  gdouble fluid_rho = sim->physical_params.alpha ? 1./
    gfs_function_value (sim->physical_params.alpha, cell) : 1.;
  GfsVariable ** u = gfs_domain_velocity (domain);
 
  gdouble viscosity = 0.;
  GfsSourceDiffusion * d = source_diffusion_viscosity (u[0]); 
  if (d) viscosity = gfs_diffusion_cell (d->D, cell);
  
  FttVector fluid_vel;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&fluid_vel.x)[c] = gfs_interpolate (cell, p->pos, u[c]);

  FttVector relative_vel = subs_fttvectors (&fluid_vel, &particulate->vel);
  FttVector vorticity;
  vorticity_vector (cell, u, &vorticity);

  gdouble cl = 0.5;
  if (coeff->coefficient) {
    gdouble norm_relative_vel = sqrt (relative_vel.x*relative_vel.x + 
				      relative_vel.y*relative_vel.y +
				      relative_vel.z*relative_vel.z);
    gdouble dia =  2.*pow(3.0*(particulate->volume)/4.0/M_PI, 1./3.);
    if (viscosity == 0) {
      g_warning ("Viscosity is 0. cannot compute lift force on particulate\n");
      g_assert_not_reached ();
    }
    gdouble Re = norm_relative_vel*dia*fluid_rho/viscosity;

    GFS_VALUE (cell, coeff->re_p) = Re;
    GFS_VALUE (cell, coeff->pdia) = dia;
    GFS_VALUE (cell, coeff->u_rel) = relative_vel.x;
    GFS_VALUE (cell, coeff->v_rel) = relative_vel.y;
#if !FTT_2D
    GFS_VALUE (cell, coeff->w_rel) = relative_vel.z;
#endif
    cl = gfs_function_value (coeff->coefficient, cell); 
  }
 
#if FTT_2D
  force.x = fluid_rho*cl*relative_vel.y*vorticity.z;
  force.y = -fluid_rho*cl*relative_vel.x*vorticity.z;
#else
  force.x = fluid_rho*cl*(relative_vel.y*vorticity.z
			  -relative_vel.z*vorticity.y);
  force.y = fluid_rho*cl*(relative_vel.z*vorticity.x
			  -relative_vel.x*vorticity.z);
  force.z = fluid_rho*cl*(relative_vel.x*vorticity.y
			  -relative_vel.y*vorticity.x);
#endif

  return force; 
}

static void gfs_force_lift_init (GfsParticleForce * force)
{
  force->force = compute_lift_force;
}

GtsSListContaineeClass * gfs_force_lift_class (void)
{
  static GtsSListContaineeClass * klass = NULL;
  
  if (klass == NULL) {
    GtsObjectClassInfo gfs_force_lift_info = {
      "GfsForceLift",
      sizeof (GfsForceCoeff),
      sizeof (GtsSListContaineeClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_force_lift_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_force_coeff_class ()),
				  &gfs_force_lift_info);
  }
  return klass;
}

/* GfsForceDrag: object */

static FttVector compute_drag_force (GfsParticle * p, GfsParticleForce * dragforce)
{
  GfsParticulate * particulate = GFS_PARTICULATE (p);
  GfsForceCoeff * coeff = FORCE_COEFF (dragforce);
  GfsSimulation * sim = gfs_object_simulation (particulate);
  GfsDomain * domain = GFS_DOMAIN (sim);

  FttVector force;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] = 0;

  FttCell * cell = gfs_domain_locate (domain, p->pos, -1, NULL);
  if (cell == NULL) return force;

  gdouble fluid_rho = sim->physical_params.alpha ? 1./
    gfs_function_value (sim->physical_params.alpha,cell) : 1.;
  GfsVariable ** u = gfs_domain_velocity (domain);  

  gdouble viscosity = 0.;  

  GfsSourceDiffusion * d = source_diffusion_viscosity (u[0]); 
  if (d) viscosity = gfs_diffusion_cell (d->D, cell);
  
  FttVector fluid_vel;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&fluid_vel.x)[c] = gfs_interpolate (cell, p->pos, u[c]);

  FttVector relative_vel = subs_fttvectors (&fluid_vel, &particulate->vel);

  gdouble dia = 2.*pow(3.0*(particulate->volume)/4.0/M_PI, 1./3.);
#if !FTT_2D
  gdouble norm_relative_vel = sqrt (relative_vel.x*relative_vel.x + 
				    relative_vel.y*relative_vel.y +
				    relative_vel.z*relative_vel.z);
#else
  gdouble norm_relative_vel = sqrt (relative_vel.x*relative_vel.x + 
				    relative_vel.y*relative_vel.y);
#endif

  gdouble cd = 0.;
  gdouble Re;
  if (viscosity == 0)    
    return force;
  else
    Re = norm_relative_vel*dia*fluid_rho/viscosity;
  
  if (coeff->coefficient) {
    GFS_VALUE (cell, coeff->re_p) = Re;
    GFS_VALUE (cell, coeff->u_rel) = relative_vel.x;
    GFS_VALUE (cell, coeff->v_rel) = relative_vel.y;
#if !FTT_2D
    GFS_VALUE (cell, coeff->w_rel) = relative_vel.z;
#endif
    GFS_VALUE (cell, coeff->pdia) = dia;
    cd = gfs_function_value (coeff->coefficient, cell); 
  }
  else {
    if (Re < 1e-8)
      return force;
    else if (Re < 50.0)
      cd = 16.*(1. + 0.15*pow(Re,0.5))/Re;
    else
      cd = 48.*(1. - 2.21/pow(Re,0.5))/Re;
  }
  for (c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] += 3./(4.*dia)*cd*norm_relative_vel*(&relative_vel.x)[c]*fluid_rho;
  
  return force;
}

static void gfs_force_drag_init (GfsParticleForce * force)
{
  force->force = compute_drag_force;
}

GtsSListContaineeClass * gfs_force_drag_class (void)
{
  static GtsSListContaineeClass * klass = NULL;
  
  if (klass == NULL) {
    GtsObjectClassInfo gfs_force_drag_info = {
      "GfsForceDrag",
      sizeof (GfsForceCoeff),
      sizeof (GtsSListContaineeClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_force_drag_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_force_coeff_class ()),
				  &gfs_force_drag_info);
  }
  return klass;
}

/* GfsForceBuoy: object */

static FttVector compute_buoyancy_force (GfsParticle * p, GfsParticleForce * buoyforce)
{
  GfsParticulate * particulate = GFS_PARTICULATE (p); 
  GfsSimulation * sim = gfs_object_simulation (particulate);
  GfsDomain * domain = GFS_DOMAIN (sim);

  FttVector force;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] = 0;

  FttCell * cell = gfs_domain_locate (domain, p->pos, -1, NULL);
  if (cell == NULL) return force;

  gdouble fluid_rho = sim->physical_params.alpha ? 1./
    gfs_function_value (sim->physical_params.alpha, cell) : 1.;
  GfsVariable ** u = gfs_domain_velocity (domain);
 
  gdouble g[3];
  for (c = 0; c < FTT_DIMENSION; c++) {
    g[c] = 0.;
    if (u[c]->sources) {
      GSList * i = GTS_SLIST_CONTAINER (u[c]->sources)->items;
      
      while (i) {
	if (GFS_IS_SOURCE (i->data)) {
	  g[c] += gfs_function_value (GFS_SOURCE ((GfsSourceGeneric *) i->data)->intensity, 
				      cell);
	}
	i = i->next;
      }
    }
  }

  for (c = 0; c < FTT_DIMENSION; c++)
    (&force.x)[c] += (particulate->mass/particulate->volume-fluid_rho)*g[c];

  return force;
}

static void gfs_force_buoy_init (GfsParticleForce * force)
{
  force->force = compute_buoyancy_force;
}

GtsSListContaineeClass * gfs_force_buoy_class (void)
{
  static GtsSListContaineeClass * klass = NULL;
  
  if (klass == NULL) {
    GtsObjectClassInfo gfs_force_buoy_info = {
      "GfsForceBuoy",
      sizeof (GfsParticleForce),
      sizeof (GtsSListContaineeClass), 
      (GtsObjectClassInitFunc) NULL, 
      (GtsObjectInitFunc) gfs_force_buoy_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particle_force_class ()),
				  &gfs_force_buoy_info);
  }
  return klass;
}

/* GfsParticleForce: object */

static void gfs_particle_force_read (GtsObject ** o, GtsFile * fp)
{ 
  GtsObjectClass *klass;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsParticleClass)");
    return;
  }

  klass = gfs_object_class_from_name (fp->token->str);
  if (klass == NULL) {
    gts_file_error (fp, "unknown class `%s'", fp->token->str);
    return;
  }
  if (!gts_object_class_is_from_class (klass, gfs_particle_force_class ())) {
    gts_file_error (fp, "`%s' is not a GfsParticleForce", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
}

static void gfs_particle_force_write (GtsObject * o, FILE * fp)
{
  fprintf (fp, "%s", o->klass->info.name);
}

static void gfs_particle_force_class_init (GtsObjectClass * klass)
{
  GTS_OBJECT_CLASS(klass)->read = gfs_particle_force_read;
  GTS_OBJECT_CLASS(klass)->write = gfs_particle_force_write;
}

GtsSListContaineeClass * gfs_particle_force_class (void)
{
  static GtsSListContaineeClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_particle_force_info = {
      "GfsParticleForce",
      sizeof (GfsParticleForce),
      sizeof (GtsSListContaineeClass),
      (GtsObjectClassInitFunc) gfs_particle_force_class_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_slist_containee_class ()),
				  &gfs_particle_force_info);
  }
  return klass;
}

/* GfsParticulate: Object */

static void compute_forces (GfsParticleForce * event, GfsParticulate * p)
{ 
  FttComponent c;
  FttVector new_force = (event->force) (GFS_PARTICLE (p), event);
  FttVector total_force;
     
  for ( c = 0 ; c < FTT_DIMENSION; c++)
    (&total_force.x)[c] = (&new_force.x)[c]*p->volume + (&p->force.x)[c];
    
  p->force = total_force;
}

static gboolean gfs_particulate_event (GfsEvent * event, 
				       GfsSimulation * sim)
{
  GfsParticle * p = GFS_PARTICLE (event);
  GfsParticulate * particulate = GFS_PARTICULATE (event);

  if (particulate->forces == NULL)
    (* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_particulate_class ())->parent_class)->event)
      (event, sim);
  else {
    FttVector pos = p->pos;
    gfs_simulation_map (sim, &pos);

    FttComponent c;
    /* Velocity Verlet Algorithm */
    for (c = 0; c < FTT_DIMENSION; c++) {
      (&pos.x)[c] += (&particulate->force.x)[c]*sim->advection_params.dt*sim->advection_params.dt
	/particulate->mass/2.+ (&particulate->vel.x)[c]*sim->advection_params.dt;
      (&particulate->vel.x)[c] += (&particulate->force.x)[c]*sim->advection_params.dt
	/(2.*particulate->mass);
    }

    /* Compute forces */
    for (c = 0; c < FTT_DIMENSION; c++)
      (&particulate->force.x)[c] = 0.;      
    gts_container_foreach (GTS_CONTAINER (particulate->forces), 
			   (GtsFunc) compute_forces, particulate);

    for (c = 0; c < FTT_DIMENSION; c++)
      (&particulate->vel.x)[c] += 
	(&particulate->force.x)[c]*sim->advection_params.dt/(2.*particulate->mass);

    gfs_simulation_map_inverse (sim, &pos);
    p->pos = pos;   
  }
  return TRUE;
} 

static void gfs_particulate_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_particulate_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_particulate_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;
  GfsParticulate * p = GFS_PARTICULATE (*o);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (mass)");
    return;
  }
  p->mass = atof (fp->token->str);
  gts_file_next_token (fp);
  
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (volume)");
    return;
  }
  gdouble L = gfs_object_simulation (*o)->physical_params.L;
// L is not correctly initialized. WHY???
//  printf("L %g \n", L);
  p->volume = atof (fp->token->str);
//  p->volume /= pow(L, FTT_DIMENSION);
  gts_file_next_token (fp);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (v.x)");
    return;
  }
  p->vel.x = atof (fp->token->str)/L;
  gts_file_next_token (fp);
  
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (v.y)");
    return;
  }
  p->vel.y = atof (fp->token->str)/L;
  gts_file_next_token (fp);
  
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (v.z)");
    return;
  }
  p->vel.z = atof (fp->token->str)/L;
  gts_file_next_token (fp);

  if (fp->type == GTS_INT || fp->type == GTS_FLOAT) {
    p->force.x = atof (fp->token->str)/L;
    gts_file_next_token (fp);
  }
  
  if (fp->type == GTS_INT || fp->type == GTS_FLOAT) {
    p->force.y = atof (fp->token->str)/L;
    gts_file_next_token (fp);
  }
  
  if (fp->type == GTS_INT || fp->type == GTS_FLOAT) {
    p->force.z = atof (fp->token->str)/L;
    gts_file_next_token (fp);
  }
}

static void gfs_particulate_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_particulate_class ())->parent_class->write) (o, fp);
 
 GfsParticulate * p = GFS_PARTICULATE (o);
  gdouble L = gfs_object_simulation (o)->physical_params.L;
  //fprintf (fp, " %g %g %g %g %g", p->mass, p->volume*pow(L, FTT_DIMENSION), 
  //             p->vel.x*L, p->vel.y*L, p->vel.z*L);
  fprintf (fp, " %g %g %g %g %g", p->mass, p->volume, 
               p->vel.x*L, p->vel.y*L, p->vel.z*L);
  fprintf (fp, " %g %g %g", p->force.x*L, p->force.y*L, p->force.z*L);
}

static void gfs_particulate_class_init (GfsEventClass * klass)
{
  klass->event = gfs_particulate_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_particulate_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_particulate_write;
}

GfsEventClass * gfs_particulate_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_particulate_info = {
      "GfsParticulate",
      sizeof (GfsParticulate),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_particulate_class_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particle_class ()),
				  &gfs_particulate_info);
  }
  return klass;
}

/* GfsParticleList: Object */

static void assign_forces(GfsParticulate *particulate, GtsSListContainer *forces)
{
  particulate->forces = forces;
}

static void gfs_particle_list_read (GtsObject ** o, GtsFile * fp)
{
  if (GTS_OBJECT_CLASS (gfs_particle_list_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_particle_list_class ())->parent_class->read) 
      (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsParticleList * p = GFS_PARTICLE_LIST (*o);  
  GfsEventList * l = GFS_EVENT_LIST (p);  
  if (fp->type == '{') {
    fp->scope_max++;
    gts_file_next_token (fp);
    
    while (fp->type == '\n')
      gts_file_next_token (fp);
  
    GfsSimulation * sim = gfs_object_simulation (*o);
    GtsObjectClass * klass;
    while (fp->type != '}') {      

      if (fp->type != GTS_STRING) {
	gts_file_error (fp, "expecting a keyword (GfsParticleForce)");
	break;
      }
      klass = gfs_object_class_from_name (fp->token->str);
 
      if (klass == NULL) {
	gts_file_error (fp, "unknown class `%s'", fp->token->str);
	break;
      }
      if (!gts_object_class_is_from_class (klass, gfs_particle_force_class ())) {
	gts_file_error (fp, "'%s' is not a GfsParticleForce", fp->token->str);
	break;
      }
  
      GtsObject * object = gts_object_new (klass);
      gfs_object_simulation_set (object, sim);
  
      (* klass->read) (&object, fp);
    
      if (fp->type == GTS_ERROR) {
	gts_object_destroy (object);
	break;
      }
  
      while (fp->type == '\n') 
	gts_file_next_token (fp);
      
      gts_container_add (GTS_CONTAINER (p->forces), GTS_CONTAINEE (object));   
    }
    
    if (fp->type != '}') {
      gts_file_error (fp, "expecting a closing brace");
      return;
    }
    fp->scope_max--;
    gts_file_next_token (fp); 
  }

  if (p->forces->items != NULL) {
    p->forces->items = g_slist_reverse (p->forces->items);
    gts_container_foreach (GTS_CONTAINER (l->list), (GtsFunc) assign_forces, p->forces);
  }

  if(fp->type == GTS_INT){
    p->idlast = atoi (fp->token->str);
    gts_file_next_token (fp);
  }    
}

static void gfs_particle_list_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_particle_list_class ())->parent_class->write) (o, fp);

  GfsParticleList * p = GFS_PARTICLE_LIST (o);
  fputs (" {\n", fp);
  GSList * i = p->forces->items;
  while (i) {
    fputs ("    ", fp);
    (* GTS_OBJECT (i->data)->klass->write) (i->data, fp);
    fputc ('\n', fp);
    i = i->next; 
  }
  fputc ('}', fp);

  fprintf (fp, " %d", p->idlast);
}

static void gfs_particle_list_init (GtsObject *o){

  GfsParticleList * plist = GFS_PARTICLE_LIST(o);

  plist->forces = 
    GTS_SLIST_CONTAINER (gts_container_new (GTS_CONTAINER_CLASS (gts_slist_container_class ())));
 
}

static void gfs_particle_list_destroy (GtsObject * o)
{
  GfsParticleList * plist = GFS_PARTICLE_LIST(o);
 
  gts_container_foreach (GTS_CONTAINER (plist->forces), (GtsFunc) gts_object_destroy, NULL);
  gts_object_destroy (GTS_OBJECT (plist->forces));

  (* GTS_OBJECT_CLASS (gfs_particle_list_class ())->parent_class->destroy) (o);
}

static void gfs_particle_list_class_init (GtsObjectClass * klass)
{
  klass->read = gfs_particle_list_read;
  klass->write = gfs_particle_list_write;  
  klass->destroy = gfs_particle_list_destroy;  
}

GfsEventClass * gfs_particle_list_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_particle_list_info = {
      "GfsParticleList",
      sizeof (GfsParticleList),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_particle_list_class_init,
      (GtsObjectInitFunc) gfs_particle_list_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_event_list_class ()),
				  &gfs_particle_list_info);
  }
  return klass;
}

typedef struct {
  FttVector pos, vel;
  gdouble volume;
} Droplets;

typedef struct {
  GfsVariable * tag, * c, *t;
  Droplets * drops;
  GfsVariable **u;
  guint * sizes;
  guint n, min;
  gdouble resetval;
  gdouble density;
  GfsFunction *fc;
} DropletsPar;

static int greater (const void * a, const void * b)
{
  return *((guint *)a) > *((guint *)b) ? -1 : 1;
}

static void reset_small_fraction (FttCell * cell, DropletsPar * p)
{
  gint i = GFS_VALUE (cell, p->tag);
  if (i > 0 && p->sizes[i - 1] < p->min)
    GFS_VALUE (cell, p->c) = p->resetval;
}

static void compute_droplet_properties (FttCell * cell, DropletsPar * p)
{
  gint i = GFS_VALUE (cell, p->tag);
  gdouble h = ftt_cell_size (cell), vol;
  FttVector pos; 
  ftt_cell_pos (cell, &pos);
  GfsVariable ** u = p->u;

  if (i > 0) {
    p->sizes[i - 1]++;
    vol = pow (h, FTT_DIMENSION);
    p->drops[i-1].volume += vol*GFS_VALUE (cell, p->c);
    FttComponent c;
    for(c = 0; c < FTT_DIMENSION; c++){
      (&(p->drops[i-1].pos.x))[c] +=  (&pos.x)[c];
      (&(p->drops[i-1].vel.x))[c] += GFS_VALUE (cell,u[c]);
    }
  }  
}

static void convert_droplets (GfsDomain * domain, 
			      DropletsPar * pars, GfsParticleList * plist)
{
  GfsSimulation * sim = gfs_object_simulation (plist); 
  guint i;
  
  GfsDropletToParticle * d = DROPLET_TO_PARTICLE (plist);
  GfsEventList * l = GFS_EVENT_LIST (plist); 

  pars->sizes = g_malloc0 (pars->n*sizeof (guint));  
  pars->drops = g_malloc0 (pars->n*sizeof (Droplets));

  FttComponent c;
  /* Initialize drops */
  for (i = 0; i < pars->n; i++){
    pars->drops[i].volume = 0.;
    pars->sizes[i] = 0;
    for(c = 0; c < FTT_DIMENSION; c++) {
      (&(pars->drops[i].pos.x))[c] = 0.;
      (&(pars->drops[i].vel.x))[c] = 0.;
    }
  }
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) compute_droplet_properties, pars);

#ifdef HAVE_MPI
  if (domain->pid >= 0) {
    guint * sizes = g_malloc0 (pars->n*sizeof (guint));
    MPI_Allreduce (pars->sizes, sizes, pars->n, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    g_free (pars->sizes);
    pars->sizes = sizes;
  }
#endif
  if (d->min >= 0)
    pars->min = d->min;
  else {
    guint * tmp = g_malloc (pars->n*sizeof (guint));
    memcpy (tmp, pars->sizes, pars->n*sizeof (guint));
    qsort (tmp, pars->n, sizeof (guint), greater);
    g_assert (-1 - d->min < pars->n);
    /* fixme: this won't work for parallel jobs */
    pars->min = tmp[-1 - d->min];
    g_free (tmp);
  }
  
  for (i = 0; i < pars->n; i++) {
    if (pars->sizes[i] < pars->min){
      for (c = 0; c < FTT_DIMENSION; c++) {
      	(&pars->drops[i].pos.x)[c] = (&pars->drops[i].pos.x)[c]/pars->sizes[i];
	(&pars->drops[i].vel.x)[c] = (&pars->drops[i].vel.x)[c]/pars->sizes[i];
      }
      FttCell * cell = gfs_domain_locate (domain, pars->drops[i].pos, -1, NULL);    
      if (cell) {
	/* Construct an Object */
	GtsObjectClass * klass = l->klass;
	if (klass == NULL) {
	  gfs_error (0, "Unknown particle class\n");
	  return;
	}
	GtsObject * object = gts_object_new (klass);
	gfs_object_simulation_set (object, sim);
	l->list->items = g_slist_reverse (l->list->items);	
	gts_container_add (GTS_CONTAINER (l->list), GTS_CONTAINEE (object));
	l->list->items = g_slist_reverse (l->list->items);
	GfsEvent * list = GFS_EVENT (l);	
	gfs_event_set (GFS_EVENT (object), 
		       list->start, list->end, list->step, list->istart, list->iend, list->istep);
	GfsParticulate * drop = GFS_PARTICULATE (object);
	GfsParticle * p = GFS_PARTICLE (drop);
	
	drop->vel = pars->drops[i].vel;
	p->pos = pars->drops[i].pos;
	drop->volume = pars->drops[i].volume;
	p->id = ++plist->idlast;
	drop->mass = sim->physical_params.alpha ? 1./
	  gfs_function_value (sim->physical_params.alpha, cell) : 1.;
	drop->mass *= drop->volume;
	for (c = 0; c < FTT_DIMENSION; c++)
	  (&drop->force.x)[c] = 0.;
      }       
    }   
  }  
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) reset_small_fraction, pars); 
  g_free (pars->drops);
  g_free (pars->sizes);
 
}

/* GfsDropletToParticle: object */

static void compute_v (FttCell * cell, GfsRemoveDroplets * d)
{
  GFS_VALUE (cell, d->v) = gfs_function_value (d->fc, cell);
}

static gboolean gfs_droplet_to_particle_event (GfsEvent * event, GfsSimulation * sim)
{ 
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_droplet_to_particle_class ())->parent_class)->event)
      (event, sim)) {  
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsParticleList * plist = GFS_PARTICLE_LIST (event);
    GfsDropletToParticle *d = DROPLET_TO_PARTICLE (event);
    d->v = d->fc ? gfs_function_get_variable (d->fc) : d->c;
    DropletsPar p ;   
  
    p.resetval = d->resetwith;
    p.tag = gfs_temporary_variable (domain);
    p.u = gfs_domain_velocity (domain);
    p.density = d->density;
    p.t = d->c;
  
    if (d->v){
      p.c = d->v;
      p.n = gfs_domain_tag_droplets (domain, p.c, p.tag);
      if (p.n > 0 && -d->min < (gint) p.n){
	p.c = d->c;
	convert_droplets (domain, &p, plist);
      }
    }
    else {      
      d->v = gfs_temporary_variable (domain);      
      gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
				(FttCellTraverseFunc) compute_v, d);
      p.c = d->v;      
      p.n = gfs_domain_tag_droplets (domain, p.c, p.tag);
      if (p.n > 0 && -d->min < (gint) p.n){
	p.c = d->c;	
	convert_droplets (domain, &p, plist);	
      }              
      gts_object_destroy (GTS_OBJECT (d->v));      
    } 

    gts_object_destroy (GTS_OBJECT (p.tag));
    return TRUE;
  }
  return FALSE;
}

static void gfs_droplet_to_particle_read (GtsObject ** o, GtsFile * fp)
{  
  if (GTS_OBJECT_CLASS (gfs_droplet_to_particle_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_droplet_to_particle_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsDropletToParticle * r = DROPLET_TO_PARTICLE (*o);  
  GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (r));

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (variable)");
    return;
  }

  if ((r->c = gfs_variable_from_name (domain->variables, fp->token->str)) == NULL) {
    gts_file_error (fp, "unknown variable `%s'", fp->token->str);
    return;
  }
  gts_file_next_token (fp);

  if (fp->type == '{') {
    GtsFileVariable var[] = {
      {GTS_INT, "min",  TRUE},
      {GTS_DOUBLE, "reset",    TRUE},
      {GTS_DOUBLE, "density",   TRUE},
      {GTS_NONE}
    };

    var[0].data = &r->min;
    var[1].data = &r->resetwith;
    var[2].data = &r->density;

    gts_file_assign_variables (fp, var);
  }
 
  if (fp->type != '\n') {
    r->fc = gfs_function_new (gfs_function_class (), 0.);
    gfs_function_read (r->fc, gfs_object_simulation (r), fp);
  }
}

static void gfs_droplet_to_particle_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_droplet_to_particle_class ())->parent_class->write) (o, fp);

  GfsDropletToParticle * r = DROPLET_TO_PARTICLE(o);

  fprintf (fp, " %s { min = %d reset = %g density = %g } ",
	   r->c->name, r->min, r->resetwith, r->density);
  if (r->fc)
    gfs_function_write (r->fc, fp);
}

static void gfs_droplet_to_particle_destroy (GtsObject * o)
{
  GfsDropletToParticle * drops = DROPLET_TO_PARTICLE (o);
  if (drops->fc)
    gts_object_destroy (GTS_OBJECT (drops->fc));
  
  (* GTS_OBJECT_CLASS (gfs_droplet_to_particle_class ())->parent_class->destroy) (o);
}

static void gfs_droplet_to_particle_init (GfsDropletToParticle * r)
{
  r->resetwith = 0.;
  r->min = 20;
  r->density = 1.;
}

static void gfs_droplet_to_particle_class_init (GfsEventClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = gfs_droplet_to_particle_event;
  GTS_OBJECT_CLASS (klass)->read = gfs_droplet_to_particle_read;
  GTS_OBJECT_CLASS (klass)->write = gfs_droplet_to_particle_write;  
  GTS_OBJECT_CLASS (klass)->destroy = gfs_droplet_to_particle_destroy;  
}

GfsEventClass * gfs_droplet_to_particle_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_droplet_to_particle_info = {
      "GfsDropletToParticle",
      sizeof (GfsDropletToParticle),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_droplet_to_particle_class_init,
      (GtsObjectInitFunc) gfs_droplet_to_particle_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particle_list_class ()),
				  &gfs_droplet_to_particle_info);
  }
  return klass;
}

/* GfsParticulateField: object */

static void voidfraction_from_particles (FttCell * cell, GfsVariable * v, GfsParticulate * part)
{
  GFS_VALUE (cell, v) += part->volume/ftt_cell_volume (cell);
}

static gboolean particulate_field_event (GfsEvent * event, 
					 GfsSimulation * sim) 
{
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_particulate_field_class ())->parent_class)->event)
      (event, sim)) {
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsVariable * v = GFS_VARIABLE (event);
    GfsParticulateField * pfield = GFS_PARTICULATE_FIELD (v);

    /* Reset variable */
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_reset, v);
    /* Loop over the list of particles in the selected object */
    GSList * i = GFS_EVENT_LIST (pfield->plist)->list->items;
    while (i) {
      FttCell * cellpart = gfs_domain_locate (domain, GFS_PARTICLE (i->data)->pos, -1, NULL);
      if (cellpart)
	pfield->voidfraction_func (cellpart, v, i->data);
      i = i->next;
    }
    return TRUE;
  }
  return FALSE;
}

static void particulate_field_read (GtsObject ** o, GtsFile * fp)
{
  (* GTS_OBJECT_CLASS (gfs_particulate_field_class ())->parent_class->read) (o, fp); 
  if (fp->type == GTS_ERROR)
    return;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (object name)");
    return;
  }

  GfsParticulateField * pfield = GFS_PARTICULATE_FIELD (*o);
  GtsObject * object = gfs_object_from_name (GFS_DOMAIN (gfs_object_simulation (*o)), 
					     fp->token->str);
  if (object == NULL) {
    gts_file_error (fp, "unknown object '%s'", fp->token->str);
    return;
  }
  if (!GFS_IS_PARTICLE_LIST (object)) {
    gts_file_error (fp, "object '%s' is not a GfsParticleList", fp->token->str);
    return;  
  }
  gts_file_next_token (fp);
  
  pfield->plist = GFS_PARTICLE_LIST (object);
}

static void particulate_field_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_particulate_field_class ())->parent_class->write) (o, fp); 
  fprintf (fp, " %s", GFS_EVENT (GFS_PARTICULATE_FIELD (o)->plist)->name);
}

static void particulate_field_class_init (GtsObjectClass * klass)
{
  GFS_EVENT_CLASS (klass)->event = particulate_field_event;
  klass->read =  particulate_field_read;
  klass->write = particulate_field_write;
}

static void particulate_field_init (GfsVariable * v)
{
  v->units = -FTT_DIMENSION;
  GFS_PARTICULATE_FIELD (v)->voidfraction_func = voidfraction_from_particles;
}

GfsVariableClass * gfs_particulate_field_class (void)
{
  static GfsVariableClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_particulate_field_info = {
      "GfsParticulateField",
      sizeof (GfsParticulateField),
      sizeof (GfsVariableClass),
      (GtsObjectClassInitFunc) particulate_field_class_init,
      (GtsObjectInitFunc) particulate_field_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_variable_class ()),
                                  &gfs_particulate_field_info);
  }

  return klass;
}

/** \beginobject{GfsFeedParticle} */

static void add_particulate (GfsDomain * domain, 
			     GfsFeedParticle * feedp, GfsParticleList * plist)
{
  GfsSimulation * sim = gfs_object_simulation (plist); 
  GfsEventList * l = GFS_EVENT_LIST (plist); 
  guint c;
  FttVector pos,vel;
  
  pos.x = gfs_function_value (feedp->posx, NULL); 
  pos.y = gfs_function_value (feedp->posy, NULL); 
  pos.z = gfs_function_value (feedp->posz, NULL); 
  FttCell * cell = gfs_domain_locate (domain, pos, -1, NULL);    
  if (cell) {
    /* Construct an Object */
    GtsObjectClass * klass = l->klass;
    g_assert (klass);
    GtsObject * object = gts_object_new (klass);
    gfs_object_simulation_set (object, sim);
    l->list->items = g_slist_reverse (l->list->items);	
    gts_container_add (GTS_CONTAINER (l->list), GTS_CONTAINEE (object));
    l->list->items = g_slist_reverse (l->list->items);
    GfsEvent * list = GFS_EVENT (l);	
    gfs_event_set (GFS_EVENT (object), 
                   list->start, list->end, list->step, list->istart, list->iend, list->istep);
    GfsParticulate * part = GFS_PARTICULATE (object);
    GfsParticle * p = GFS_PARTICLE (part);
    
    vel.x = gfs_function_value (feedp->velx, cell);
    vel.y = gfs_function_value (feedp->vely, cell);
    vel.z = gfs_function_value (feedp->velz, cell);
    part->vel = vel;
    p->pos = pos;
    part->volume = gfs_function_value (feedp->vol, cell);
    p->id = ++plist->idlast;
    part->mass = gfs_function_value (feedp->mass, cell);
    for (c = 0; c < FTT_DIMENSION; c++)
      (&part->force.x)[c] = 0.;
    assign_forces ( part , plist->forces);
  }       
}

static gboolean gfs_feed_particle_event (GfsEvent * event, GfsSimulation * sim)
{ 
  if ((* GFS_EVENT_CLASS (GTS_OBJECT_CLASS (gfs_feed_particle_class ())->parent_class)->event)
      (event, sim)) {  
    GfsDomain * domain = GFS_DOMAIN (sim);
    GfsParticleList * plist = GFS_PARTICLE_LIST (event);
    GfsFeedParticle * feedp = GFS_FEED_PARTICLE (event);
    gint i;
    guint np = gfs_function_value (feedp->np, NULL);

    for (i = 0; i < np; i++)
      add_particulate (domain, feedp, plist);
    return TRUE;
  }
  return FALSE;
}

static void gfs_feed_particle_destroy (GtsObject * o)
{
  gts_object_destroy (GTS_OBJECT (GFS_FEED_PARTICLE (o)->np));
  gts_object_destroy (GTS_OBJECT (GFS_FEED_PARTICLE (o)->posx));
  gts_object_destroy (GTS_OBJECT (GFS_FEED_PARTICLE (o)->posy));
  gts_object_destroy (GTS_OBJECT (GFS_FEED_PARTICLE (o)->posz));
  gts_object_destroy (GTS_OBJECT (GFS_FEED_PARTICLE (o)->velx));
  gts_object_destroy (GTS_OBJECT (GFS_FEED_PARTICLE (o)->vely));
  gts_object_destroy (GTS_OBJECT (GFS_FEED_PARTICLE (o)->velz));
  gts_object_destroy (GTS_OBJECT (GFS_FEED_PARTICLE (o)->mass));
  gts_object_destroy (GTS_OBJECT (GFS_FEED_PARTICLE (o)->vol));

  (* GTS_OBJECT_CLASS (gfs_feed_particle_class ())->parent_class->destroy) (o); 
}

static void gfs_feed_particle_read (GtsObject ** o, GtsFile * fp)
{  
  if (GTS_OBJECT_CLASS (gfs_feed_particle_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_feed_particle_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  GfsFeedParticle * feedp = GFS_FEED_PARTICLE(*o);

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
  else if (!strcmp (fp->token->str, "nparts")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }    
      gts_file_next_token (fp);
      gfs_function_read (feedp->np, gfs_object_simulation (*o), fp);    
    }
    else if (!strcmp (fp->token->str, "xfeed")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      gfs_function_read (feedp->posx, gfs_object_simulation (*o), fp);    
    }
    else if (!strcmp (fp->token->str, "yfeed")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      gfs_function_read (feedp->posy, gfs_object_simulation (*o), fp);    
    }
    else if (!strcmp (fp->token->str, "zfeed")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      gfs_function_read (feedp->posz, gfs_object_simulation (*o), fp);    
    }
    else if (!strcmp (fp->token->str, "velx")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      gfs_function_read (feedp->velx, gfs_object_simulation (*o), fp);    
    }
    else if (!strcmp (fp->token->str, "vely")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      gfs_function_read (feedp->vely, gfs_object_simulation (*o), fp);    
    }
    else if (!strcmp (fp->token->str, "velz")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      gfs_function_read (feedp->velz, gfs_object_simulation (*o), fp);    
    }
    else if (!strcmp (fp->token->str, "mass")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      gfs_function_read (feedp->mass, gfs_object_simulation (*o), fp);    
    }
    else if (!strcmp (fp->token->str, "volume")) {
      gts_file_next_token (fp);
      if (fp->type != '=') {
        gts_file_error (fp, "expecting '='");
        return;
      }
      gts_file_next_token (fp);
      gfs_function_read (feedp->vol, gfs_object_simulation (*o), fp);    
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

static void gfs_feed_particle_write (GtsObject * o, FILE * fp)
{
  (* GTS_OBJECT_CLASS (gfs_feed_particle_class ())->parent_class->write) (o, fp);

  GfsFeedParticle * feedp = GFS_FEED_PARTICLE(o);
  fputs (" {\n  nparts = ", fp);
  gfs_function_write (feedp->np, fp);
  fputs ("  xfeed =", fp);
  gfs_function_write (feedp->posx, fp);
  fputs (" yfeed =", fp);
  gfs_function_write (feedp->posy, fp);
  fputs (" zfeed =", fp);
  gfs_function_write (feedp->posz, fp);
  fputs ("\n  velx =", fp);
  gfs_function_write (feedp->velx, fp);
  fputs (" vely =", fp);
  gfs_function_write (feedp->vely, fp);
  fputs (" velz =", fp);
  gfs_function_write (feedp->velz, fp);
  fputs ("\n  mass =", fp);
  gfs_function_write (feedp->mass, fp);
  fputs ("\n  volume =", fp);
  gfs_function_write (feedp->vol, fp);
  fputs ("\n}", fp);
}

static void gfs_feed_particle_class_init (GfsEventClass * klass)
{
  GFS_EVENT_CLASS (klass)->event    = gfs_feed_particle_event;
  GTS_OBJECT_CLASS (klass)->destroy = gfs_feed_particle_destroy;
  GTS_OBJECT_CLASS (klass)->read    = gfs_feed_particle_read;
  GTS_OBJECT_CLASS (klass)->write   = gfs_feed_particle_write;  
}

static void gfs_feed_particle_init ( GfsFeedParticle * feedp)
{
  feedp->np   = gfs_function_new (gfs_function_class (), 0.);
  feedp->posx = gfs_function_new (gfs_function_class (), 0.);
  feedp->posy = gfs_function_new (gfs_function_class (), 0.);
  feedp->posz = gfs_function_new (gfs_function_class (), 0.);
  feedp->velx = gfs_function_new (gfs_function_class (), 0.);
  feedp->vely = gfs_function_new (gfs_function_class (), 0.);
  feedp->velz = gfs_function_new (gfs_function_class (), 0.);
  feedp->mass = gfs_function_new (gfs_function_class (), 0.);
  feedp->vol  = gfs_function_new (gfs_function_class (), 0.);
}

GfsEventClass * gfs_feed_particle_class (void)
{
  static GfsEventClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_feed_particle_info = {
      "GfsFeedParticle",
      sizeof (GfsFeedParticle),
      sizeof (GfsEventClass),
      (GtsObjectClassInitFunc) gfs_feed_particle_class_init,
      (GtsObjectInitFunc) gfs_feed_particle_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_particle_list_class ()),
				  &gfs_feed_particle_info);
  }
  return klass;
}

/** \endobject{GfsFeedParticle} */
