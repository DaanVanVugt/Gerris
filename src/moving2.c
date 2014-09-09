/* Gerris - The GNU Flow Solver
 * Copyright (C) 2005-2009 National Institute of Water and Atmospheric Research
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

/*
 * Code only used by second-order moving solid boundaries
 */

#define SOLD2(c, d)  (GFS_VALUE (c, sold2[d]))

static void sold2_fine_init (FttCell * parent, GfsVariable * v)
{
  FttCellChildren child;
  guint n;

  ftt_cell_children (parent, &child);
  for (n = 0; n < FTT_CELLS; n++)
    if (child.c[n])
      GFS_VALUE (child.c[n], v) = 1.;
}

static int cell_is_corner (FttCell * cell)
{
  FttDirection d, d1, d2;
  gdouble  norm;
  FttCellNeighbors neighbors;
  FttVector n1, n2;


  g_assert (cell);

  ftt_cell_neighbors (cell,&neighbors);

  d1 = d2 = -1;

  if (!GFS_IS_MIXED(cell))
    return 0;

  for (d = 0; d < FTT_NEIGHBORS; d ++)
    if (GFS_STATE(cell)->solid->s[d] != 1. && GFS_STATE(cell)->solid->s[d] != 0. && d1 == -1 && d2 == -1)
      d1 = d;
    else if (GFS_STATE(cell)->solid->s[d] != 1. && GFS_STATE(cell)->solid->s[d] != 0 && d2 == -1)
      d2 = d;
    else if (GFS_STATE(cell)->solid->s[d] != 1. && GFS_STATE(cell)->solid->s[d] != 0.)
      g_assert_not_reached ();

  if ( d1 == -1 || d2 == -1) {
    FttVector pos;
    ftt_cell_pos (cell,&pos);
    
    g_warning ("REA: %f, %f \n", pos.x, pos.y);
    g_warning ("d1: %i d2: %i  \n", d1,d2);
    
    g_assert_not_reached ();
  }

 


  gfs_solid_normal (neighbors.c[d1], &n1);
  norm = sqrt (n1.x*n1.x + n1.y*n1.y);
  if (norm != 0.) {
    n1.x /= norm;
    n1.y /= norm;
  }
  gfs_solid_normal (neighbors.c[d2], &n2);
  norm = sqrt (n2.x*n2.x + n2.y*n2.y);
  if (norm != 0.) {
    n2.x /= norm;
    n2.y /= norm;
  }

  if (d1/2 == d2/2)
    return 0;
  else {
    if (neighbors.c[d2])
      if ( neighbors.c[d1])
	if (GFS_IS_MIXED (neighbors.c[d2]) && GFS_IS_MIXED (neighbors.c[d1]))
	  if (fabs(n1.x*n2.x+n1.y*n2.y) < 0.70) {
	    if (GFS_STATE(neighbors.c[d1])->solid->s[d1] > 0 && GFS_STATE(neighbors.c[d1])->solid->s[d1] < 1)
	      return 1;
	    if (GFS_STATE(neighbors.c[d2])->solid->s[d2] > 0 && GFS_STATE(neighbors.c[d2])->solid->s[d2] < 1)
	      return 1;
	  }
    return 0;
  }
}

static int cell_was_corner (FttCell * cell, GfsVariable * old_solid_v, GfsVariable ** sold2)
{
  FttDirection d, d1, d2;
  
  g_assert (cell);

  d1 = d2 = -1;

  if (!OLD_SOLID (cell))
    return 0;

  for (d = 0; d < FTT_NEIGHBORS; d ++)
    if (SOLD2 (cell, d) != 1. && SOLD2 (cell, d) != 0. && d1 == -1 && d2 == -1)
      d1 = d;
    else if (SOLD2 (cell, d) != 1. && SOLD2 (cell, d) != 0 && d2 == -1)
      d2 = d;
    else if (SOLD2 (cell, d) != 1. && SOLD2 (cell, d) != 0.)
      g_assert_not_reached (); 

  if (d1/2 == d2/2)
    return 0;
  else {
    FttCellNeighbors neighbors;

    ftt_cell_neighbors (cell, &neighbors);

    if (neighbors.c[d1] && neighbors.c[d2]) {
      FttVector n1, n2;
      FttComponent c;
      gdouble norm;

      for (c = 0; c < FTT_DIMENSION; c++) {
	(&n1.x)[c] = (SOLD2 (neighbors.c[d1], 2*c + 1) - SOLD2 (neighbors.c[d1], 2*c));
	(&n2.x)[c] = (SOLD2 (neighbors.c[d2], 2*c + 1) - SOLD2 (neighbors.c[d2], 2*c));
      }
      norm = sqrt (n1.x*n1.x + n1.y*n1.y);
      if (norm != 0.) {
	n1.x /= norm;
	n1.y /= norm;
      }
      norm = sqrt (n2.x*n2.x + n2.y*n2.y);
      if (norm != 0.) {
	n2.x /= norm;
	n2.y /= norm;
      }    
      if (fabs(n1.x*n2.x+n1.y*n2.y) < 0.70) {
	if (SOLD2 (neighbors.c[d1], d1) > 0 && SOLD2 (neighbors.c[d1], d1) < 1)
	  return 1.;
	else if (SOLD2 (neighbors.c[d2], d2) > 0 && SOLD2 (neighbors.c[d2], d2) < 1)
	  return 1;
      }
    }
    return 0;
  }
}

static double new_fluid_old_solid (FttCell * cell, FttDirection d1, 
				   GfsVariable * old_solid,
				   GfsVariable ** sold2) 
{  
  FttDirection d;
  FttCellNeighbors neighbors;
  double s1, s2;

  g_assert(cell);

  s1 = 1.-SOLD2 (cell, d1);
  ftt_cell_neighbors (cell,&neighbors);
  for  (d = 0; d < 2*FTT_DIMENSION;d++)
    if (d != 2*(d1/2) && d != 2*(d1/2)+1)
      if (neighbors.c[d])
	if(GFS_IS_MIXED(neighbors.c[d]))
	  if (!cell_is_corner (neighbors.c[d]) && 
	      !cell_was_corner (neighbors.c[d], old_solid, sold2)) {
	    if (GFS_STATE(neighbors.c[d])->solid->s[d1] != 1.) {
	      if (SOLD2 (neighbors.c[d], d1) == 0.)
		{
		  s2 = GFS_STATE(neighbors.c[d])->solid->s[d1];
		  return s2/(s1+s2);
		}		  
	    }
	  } 
  return -1.;
}

static double new_solid_old_fluid (FttCell * cell, FttDirection d1, 
				   GfsVariable * old_solid,
				   GfsVariable ** sold2) 
{  
  FttDirection d;
  FttCellNeighbors neighbors;
  double s1, s2;

  g_assert(cell);

  s1 = 1.-GFS_STATE (cell)->solid->s[d1];
		    
  ftt_cell_neighbors (cell,&neighbors);
  for  (d = 0; d < 2*FTT_DIMENSION;d++)
    if (d != 2*(d1/2) && d != 2*(d1/2)+1)
      if (neighbors.c[d])
	if (!cell_is_corner(neighbors.c[d]) && 
	    !cell_was_corner(neighbors.c[d], old_solid, sold2))
	  if (GFS_STATE(neighbors.c[d])->solid)	 
	    if (GFS_STATE(neighbors.c[d])->solid->s[d1] == 0. && SOLD2 (neighbors.c[d], d1) != 1.) {
	      
	      s2 = SOLD2 (neighbors.c[d], d1);
	      return s1/(s1+s2);
	    }
  return -1.;
}

static double new_solid_old_solid (FttCell * cell, FttDirection d1,
				   GfsVariable * old_solid,
				   GfsVariable ** sold2)
{
  FttDirection d;
  FttCellNeighbors neighbors;
  double s1, s2;

  g_assert(cell);

  s1 = GFS_STATE (cell)->solid->s[d1];
  ftt_cell_neighbors (cell,&neighbors);
  for  (d = 0; d < 2*FTT_DIMENSION;d++)
    if (d != 2*(d1/2) && d != 2*(d1/2)+1)
      if (neighbors.c[d] &&
	  !cell_is_corner(neighbors.c[d]) && 
	  !cell_was_corner(neighbors.c[d], old_solid, sold2)) {
	if ((GFS_IS_MIXED(neighbors.c[d]) && GFS_STATE(neighbors.c[d])->solid->s[d1] == 1.) ||
	    !GFS_IS_MIXED(neighbors.c[d])) {
	  if (SOLD2 (neighbors.c[d], d1) != 1.){
	    s2 = 1.-SOLD2 (neighbors.c[d], d1);
	    return s1/(s1+s2);
	  }
	}
	else if ((GFS_STATE(cell)->solid->s[d1] == 0. && GFS_IS_MIXED(neighbors.c[d])) ) {
	  s1 = SOLD2 (cell, d1);
	  s2 = 1.-GFS_STATE(neighbors.c[d])->solid->s[d1];
	  return s2/(s1+s2);
	}
      }
  return -1.;
}

static void second_order_face_fractions (FttCell * cell, GfsSimulationMoving * sim)
{
#ifndef FTT_2D /* 3D */
  g_assert_not_implemented ();
#endif

  GfsVariable * old_solid_v = sim->old_solid;
  GfsVariable ** sold2 = sim->sold2;
  gdouble dt1, dt2, dto1, dto2, s1, s2;
  gint d1, d2, d, do1, do2;
  FttCellNeighbors neighbors;

  dt1 = dt2 = dto1 = dto2 = -2;
  d1 = d2 = do1 = do2 = -1;
  s1 = s2 = -1;

  g_assert(cell);
      
  ftt_cell_neighbors (cell,&neighbors);

  if (!OLD_SOLID (cell) && !GFS_IS_MIXED(cell))
    return;

  if (!OLD_SOLID (cell)) {
    FttDirection c;
    OLD_SOLID (cell) = g_malloc0 (sizeof (GfsSolidVector));

    OLD_SOLID (cell)->a = 1.;
    for (c = 0; c < FTT_NEIGHBORS; c++)
      OLD_SOLID (cell)->s[c] = 1.;
  }

  /* Find directions of intersection */
  if (GFS_IS_MIXED(cell))
    for (d = 0; d < FTT_NEIGHBORS; d ++) {
      if (GFS_STATE(cell)->solid->s[d] != 1. && GFS_STATE(cell)->solid->s[d] != 0. &&
	  d1 == -1 && d2 == -1)
	d1 = d;
      else if (GFS_STATE(cell)->solid->s[d] != 1. && GFS_STATE(cell)->solid->s[d] != 0 && d2 == -1)
	d2 = d;
      else if (GFS_STATE(cell)->solid->s[d] != 1. && GFS_STATE(cell)->solid->s[d] != 0.)
	g_assert_not_reached (); 
    }
  
  for (d = 0; d < FTT_NEIGHBORS; d ++) {
    if (SOLD2 (cell, d) != 1. && SOLD2 (cell, d) != 0. && do1 == -1 && do2 == -1)
      do1 = d;
    else if (SOLD2 (cell, d) != 1. && SOLD2 (cell, d) != 0 && do2 == -1)
      do2 = d;
    else if (SOLD2 (cell, d) != 1. && SOLD2 (cell, d) != 0.)
      g_assert_not_reached (); 
  }

  /* Treats easy cases */
  if (d1 != -1 && d1 == do1)
    OLD_SOLID (cell)->s[d1] = SOLD2 (cell, d1);
  if (d2 != -1 && d2 == do2)
    OLD_SOLID (cell)->s[d2] = SOLD2 (cell, d2);

  if (d1 == do1 && d2 == do2)
    return;
    
  /* Finds timescale for d1/do1 */
  if (d1 != -1) {
    if (SOLD2 (cell, d1) == 1.) {
      dt1 = new_solid_old_fluid (cell, d1, old_solid_v, sold2);
      if (dt1 == -1)
	if (neighbors.c[d1]){
	  FttDirection dop = ftt_opposite_direction[d1];
	  dt1 = new_solid_old_fluid (neighbors.c[d1], dop, old_solid_v, sold2);
	}
    }   
    else if (SOLD2 (cell, d1) == 0.){
      dt1 = new_solid_old_solid (cell, d1, old_solid_v, sold2);
    }  
  }
  
  if (do1 != -1 && do1 != d1 && do1 != d2) {
    if (GFS_IS_MIXED(cell) && GFS_STATE(cell)->solid->s[do1] == 0.)
      dto1 = new_solid_old_solid (cell, do1, old_solid_v, sold2);
    else
      dto1 = new_fluid_old_solid (cell, do1, old_solid_v, sold2);
  }

  /* Finds timescale for d2/do2 */

  if (d2 != -1) {
    if (SOLD2 (cell, d2) == 1.) {
      dt2 = new_solid_old_fluid (cell, d2, old_solid_v, sold2);
      if (dt2 == -1 && neighbors.c[d2]) {
	FttDirection dop = ftt_opposite_direction[d2];
	dt2 = new_solid_old_fluid (neighbors.c[d2], dop, old_solid_v, sold2);
      }
    }
    else if (SOLD2 (cell, d2) == 0.)
      dt2 = new_solid_old_solid (cell, d2, old_solid_v, sold2);
  }
 
  if (do2 != -1 && do2 != d1 && do2 != d2) {
    if (GFS_IS_MIXED(cell) && GFS_STATE(cell)->solid->s[do2] == 0.)
      dto2 = new_solid_old_solid (cell, do2, old_solid_v, sold2);
    else
      dto2 = new_fluid_old_solid (cell, do2, old_solid_v, sold2);
  }
 
  /* Uses time-scale from other faces if one is missing */
  if (dt1 == -1) {
    if (dto1 != -2)
      dt1 = dto1;
    else if (dt2 != -2)
      dt1 = dt2;
    else if (dto2 != -2)
      dt1 = dto2;
  }

  if (dt2 == -1) {
    if (dt1 != -2)
      dt2 = dt1;
    else if (dto2 != -2)
      dt2 = dto2;
    else if (dto1 != -2)
      dt2 = dto1;
  }

  if (dto1 == -1) {
    if (dt1 != -2)
      dto1 = dt1;
    else if (dt2 != -2)
      dto1 = dt2;
    else if (dto2 != -2)
      dto1 = dto2;
  }

  if (dto2 == -1) {
    if (dt1 != -2)
      dto2 = dt1;
    else if (dt2 != -2)
      dto2 = dt2;
    else if (dto1 != -2)
      dto2 = dto1;
  }

  /* Treats cell is corner */
  if (dt1 != -2 && dt2 != -2) {
    if (dt1 != dt2 && d1/2 != d2/2) {
      if (cell_is_corner (cell)) {
	if (dt1 < dt2)
	  dt2 = dt1;
	else
	  dt1 = dt2;
      }}}

  /* Treats cell was corner */
  if (dto1 != -2 && dto2 != -2 && 
      dto1 != dto2 && do1/2 != do2/2 &&
      cell_was_corner (cell, old_solid_v, sold2)) {
    if (dto1 < dto2)
      dto2 = dto1;
    else
      dto1 = dto2;
  }
  
  /* Compute the t^n+1/2 contribution of the face */
  if (do1 > -1)
    if (do1 != d1 && do1 != d2) {
      OLD_SOLID (cell)->s[do1]=SOLD2 (cell, do1)*(1-dto1)+dto1;
      if (neighbors.c[do1])
	if (!OLD_SOLID (neighbors.c[do1]) || !GFS_IS_MIXED(neighbors.c[do1])) {
	  if (!OLD_SOLID (neighbors.c[do1])) {
	    FttDirection c;
	    OLD_SOLID (neighbors.c[do1]) = g_malloc0 (sizeof (GfsSolidVector));
	    OLD_SOLID (neighbors.c[do1])->a = 1.;
	    for (c = 0; c < FTT_NEIGHBORS; c++)
	      OLD_SOLID (neighbors.c[do1])->s[c] = 1.;
	  }	  
	  OLD_SOLID (neighbors.c[do1])->s[ftt_opposite_direction[do1]] = 
	    SOLD2 (cell, do1)*(1-dto1)+dto1;
	}
    }

  if (do2 > -1)
    if (do2 != d1 && do2 != d2) {
      OLD_SOLID (cell)->s[do2]=SOLD2 (cell, do2)*(1-dto2)+dto2;
      if (neighbors.c[do2])
	if (!OLD_SOLID (neighbors.c[do2]) || !GFS_IS_MIXED(neighbors.c[do2])) {
	  if (!OLD_SOLID (neighbors.c[do2])) {
	    FttDirection c;
	    OLD_SOLID (neighbors.c[do2]) = g_malloc0 (sizeof (GfsSolidVector));
	    OLD_SOLID (neighbors.c[do2])->a = 1.;
	    for (c = 0; c < FTT_NEIGHBORS; c++)
	      OLD_SOLID (neighbors.c[do2])->s[c] = 1.;
	  }	  
	  OLD_SOLID (neighbors.c[do2])->s[ftt_opposite_direction[do2]] = 
	    SOLD2 (cell, do2)*(1-dto2)+dto2;
	}
    }


  if (d1 > -1) {
    if (SOLD2 (cell, d1) == 0.)
      OLD_SOLID (cell)->s[d1] = GFS_STATE(cell)->solid->s[d1]*(dt1-1.); 
    else if (SOLD2 (cell, d1) == 1.)
      OLD_SOLID (cell)->s[d1] = (dt1-1.)*GFS_STATE(cell)->solid->s[d1]+2.-dt1;
  }

  if (d2 > -1) {
    if (SOLD2 (cell, d2) == 0.)
      OLD_SOLID (cell)->s[d2] = GFS_STATE(cell)->solid->s[d2]*(dt2-1.); 
    else if (SOLD2 (cell, d2) == 1.)
      OLD_SOLID (cell)->s[d2] = (dt2-1.)*GFS_STATE(cell)->solid->s[d2]+2.-dt2;
  }

  if (d1/2 == d2/2 && do1 == -1 && do2 == -1)  /* third face has to be treated for 
						  the timescale determined on the other faces */  
    for (d = 0; d < FTT_NEIGHBORS; d ++)
      if (d/2 != d1/2 && SOLD2 (cell, d) == 0.)
	OLD_SOLID (cell)->s[d] = -1.+dt1+dt2;
    

  if (do1/2 == do2/2 && d1 == -1 && d2 == -1)
    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (d/2 != do1/2 && SOLD2 (cell, d) == 0.)
	OLD_SOLID (cell)->s[d] = -1.+dto1+dto2;
}

static void set_sold2 (FttCell * cell, GfsSimulationMoving * sim)
{
  GfsVariable * old_solid_v = sim->old_solid;
  GfsVariable ** sold2 = sim->sold2;
  FttDirection d;

  if (OLD_SOLID (cell))
    for (d = 0; d < FTT_NEIGHBORS; d++)
      SOLD2 (cell, d) = OLD_SOLID (cell)->s[d];
  else
    for (d = 0; d < FTT_NEIGHBORS; d++)
      SOLD2 (cell, d) = 1.;
}

static void redistribute_old_face_in_merged (FttCell * cell, 
					     FttCell * merged, FttDirection d, 
					     GfsVariable * old_solid_v)
{  
  gint i;
  gdouble sc, sm;
  
  g_assert (cell != NULL);
  g_assert (merged != NULL);

  sc = ftt_cell_volume(cell);
  sm = ftt_cell_volume(merged);
  
  if (sc != sm)
    printf("Face redistribution not implemented yet for adaptive grid \n");
  g_assert (sc == sm);
  
  for (i = 0; i < FTT_DIMENSION;i++)
    if (i != d/2) {
      FttCellNeighbors neighbors;
      FttVector pos;
      
      ftt_cell_pos(cell,&pos);
      ftt_cell_neighbors (merged,&neighbors);

      GfsSolidVector * old_solid_merged = OLD_SOLID (merged);
      if (!old_solid_merged) {
	FttDirection c;
	OLD_SOLID (merged) = old_solid_merged = g_malloc0 (sizeof (GfsSolidVector));
	old_solid_merged->a = 1.;
	for (c = 0; c < FTT_NEIGHBORS; c++)
	  old_solid_merged->s[c] = 1.;
      }
          
      old_solid_merged->s[2*i] += OLD_SOLID (cell)->s[2*i];

      if (neighbors.c[2*i]) {
	GfsSolidVector * old_solid = OLD_SOLID (neighbors.c[2*i]);
	if (!old_solid) {
	  FttDirection c;
	  OLD_SOLID (neighbors.c[2*i]) = old_solid = g_malloc0 (sizeof (GfsSolidVector));
	  old_solid->a = 1.;
	  for (c = 0; c < FTT_NEIGHBORS; c++)
	    old_solid->s[c] = 1.;
	}
	old_solid->s[ftt_opposite_direction[2*i]] += OLD_SOLID (cell)->s[2*i];
      }
      
      old_solid_merged->s[2*i+1] += OLD_SOLID (cell)->s[2*i+1];
      
      if (neighbors.c[2*i+1]) {
	GfsSolidVector * old_solid = OLD_SOLID (neighbors.c[2*i+1]);
	if (!old_solid) {
	  FttDirection c;
	  OLD_SOLID (neighbors.c[2*i+1]) = old_solid = g_malloc0 (sizeof (GfsSolidVector));
	  old_solid->a = 1.;
	  for (c = 0; c < FTT_NEIGHBORS; c++)
	    old_solid->s[c] = 1.;
	}
	old_solid->s[ftt_opposite_direction[2*i+1]] += OLD_SOLID (cell)->s[2*i+1];	
      }
    }
}

static void redistribute_old_face (FttCell * cell, FttCell * merged, GfsVariable * old_solid) 
{
  FttCellNeighbors neighbors;
  FttDirection d;

  ftt_cell_neighbors (cell,&neighbors);
  for (d = 0; d< FTT_NEIGHBORS; d++)
    if (neighbors.c[d])
      redistribute_old_face_in_merged (cell, neighbors.c[d], d, old_solid);
}

static double face_fraction_half (const FttCellFace * face, const GfsAdvectionParams * par)
{
  GfsVariable * old_solid_v = GFS_SIMULATION_MOVING (par->v->domain)->old_solid;
  if (face->cell && OLD_SOLID (face->cell))
    return OLD_SOLID (face->cell)->s[face->d];
  return 1.;
}

/* see gfs_face_advection_flux() for the initial implementation with static boundaries */
static void moving_face_advection_flux (const FttCellFace * face,
					const GfsAdvectionParams * par)
{
  gdouble flux;
  
  /* fixme: what's up with face mapping? */
  flux = face_fraction_half (face, par)*GFS_FACE_NORMAL_VELOCITY (face)*par->dt*
    gfs_face_upwinded_value (face, GFS_FACE_UPWINDING, NULL)/ftt_cell_size (face->cell);
  if (!FTT_FACE_DIRECT (face))
    flux = - flux;
  GFS_VALUE (face->cell, par->fv) -= flux;

  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE:
    GFS_VALUE (face->neighbor, par->fv) += flux;
    break;
  case FTT_FINE_COARSE:
    GFS_VALUE (face->neighbor, par->fv) += flux/FTT_CELLS;
    break;
  default:
    g_assert_not_reached ();
  }
}

/* see gfs_face_velocity_advection_flux() for the initial implementation with static boundaries */
static void moving_face_velocity_advection_flux (const FttCellFace * face,
						 const GfsAdvectionParams * par)
{
  gdouble flux;
  FttComponent c = par->v->component;

  g_return_if_fail (c >= 0 && c < FTT_DIMENSION);

  /* fixme: what's up with face mapping? */
  flux = face_fraction_half (face, par)*GFS_FACE_NORMAL_VELOCITY (face)*
    par->dt/ftt_cell_size (face->cell);
#if 0
  if (c == face->d/2) /* normal component */
    flux *= GFS_FACE_NORMAL_VELOCITY (face);
  else /* tangential component */
#else
    flux *= gfs_face_upwinded_value (face, par->upwinding, par->u)
      /* pressure correction */
      - gfs_face_interpolated_value (face, par->g[c]->i)*par->dt/2.;
#endif
  if (!FTT_FACE_DIRECT (face))
    flux = - flux;
  GFS_VALUE (face->cell, par->fv) -= flux;

  switch (ftt_face_type (face)) {
  case FTT_FINE_FINE:
    GFS_VALUE (face->neighbor, par->fv) += flux;
    break;
  case FTT_FINE_COARSE:
    GFS_VALUE (face->neighbor, par->fv) += flux/FTT_CELLS;
    break;
  default:
    g_assert_not_reached ();
  }
}

static void swap_fractions (FttCell * cell, GfsVariable * old_solid_v) {
  FttDirection c;
  
  g_assert (cell);
  
  if (FTT_CELL_IS_LEAF(cell)) {
    if (OLD_SOLID (cell)) {
      GfsSolidVector * solid_old = OLD_SOLID (cell);
      
      if (GFS_STATE (cell)->solid) {
	GfsSolidVector * solid = GFS_STATE (cell)->solid;
	
	OLD_SOLID (cell)->merged = GFS_STATE (cell)->solid->merged;
	
	for (c = 0; c < 2*FTT_DIMENSION; c++) 
	  if (solid->s[c] == 0.)
	    solid_old->s[c] = 0;
	  else
	    solid_old->s[c] = (solid_old->s[c]+solid->s[c])/2. ;	
      }
      else {
	OLD_SOLID (cell)->merged = NULL;
	
	for (c = 0; c < 2*FTT_DIMENSION; c++)
	  solid_old->s[c] = (solid_old->s[c]+1.)/2. ;
	
      }
    }
    else if (GFS_STATE (cell)->solid) {
      GfsSolidVector * solid = GFS_STATE (cell)->solid;
      GfsSolidVector * solid_old = OLD_SOLID (cell) = g_malloc0 (sizeof (GfsSolidVector));
      OLD_SOLID (cell)->a= 1.;
      OLD_SOLID (cell)->merged = GFS_STATE (cell)->solid->merged;
      
      for (c = 0; c < 2*FTT_DIMENSION; c++) 
	solid_old->s[c] = 1.;	
      
      for (c = 0; c < 2*FTT_DIMENSION; c++) 
	if (solid->s[c] == 0.)
	  solid_old->s[c] = 0;
	else
	  solid_old->s[c] = (solid_old->s[c]+solid->s[c])/2. ;
    }
  }
  
  if (OLD_SOLID (cell)) {
    if (GFS_STATE(cell)->solid) {
      GfsSolidVector * tmp = OLD_SOLID (cell);
      OLD_SOLID (cell)->merged = GFS_STATE (cell)->solid->merged;
      OLD_SOLID (cell) = GFS_STATE(cell)->solid;
      GFS_STATE(cell)->solid = tmp;
      tmp = NULL;
    }
    else {
      OLD_SOLID (cell)->merged = NULL;
      GFS_STATE(cell)->solid = OLD_SOLID (cell);
      OLD_SOLID (cell) = NULL;
    }
  }
  else if (GFS_STATE(cell)->solid) {
    OLD_SOLID (cell) = GFS_STATE(cell)->solid;
    GFS_STATE(cell)->solid = NULL;
  }
  
  
  /* Check for negative fractions and fix */
  if (GFS_STATE(cell)->solid)
    for (c = 0; c < 2*FTT_DIMENSION; c++)
      if (GFS_STATE(cell)->solid->s[c] < 0.) {
	if (OLD_SOLID (cell)) 
	  if (OLD_SOLID (cell)->s[c] >= 0.)
	    GFS_STATE(cell)->solid->s[c] = OLD_SOLID (cell)->s[c];
	  else
	    GFS_STATE(cell)->solid->s[c] = 1.;
	else
	  GFS_STATE(cell)->solid->s[c] = 0.;
      }
  
  if (OLD_SOLID (cell)) 
    for (c = 0; c < 2*FTT_DIMENSION; c++)
      if (OLD_SOLID (cell)->s[c] < 0.){
	if (GFS_STATE(cell)->solid)
	  if (GFS_STATE(cell)->solid->s[c] >= 0.)
	    OLD_SOLID (cell)->s[c] = GFS_STATE(cell)->solid->s[c];
	  else
	    OLD_SOLID (cell)->s[c] = 1.;
	else
	  OLD_SOLID (cell)->s[c] = 0.;
      }
}

static void old_solid_fractions_from_children (FttCell * cell)
{
  if (!FTT_CELL_IS_LEAF (cell)) {
    FttCellChildren child;
    guint i;
    
    ftt_cell_children (cell, &child);
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i])
	old_solid_fractions_from_children (child.c[i]);
    
      gfs_cell_init_solid_fractions_from_children (cell);
  }
}

static void foreach_box (GfsBox * box, gpointer data)
{
  old_solid_fractions_from_children (box->root);
}

static void swap_face_fractions (GfsSimulation * sim)
{
  GfsDomain * domain = GFS_DOMAIN (sim);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) swap_fractions, 
			    GFS_SIMULATION_MOVING (sim)->old_solid);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) foreach_box, NULL);
}

static void swap_fractions_back (FttCell * cell, GfsVariable * old_solid_v) 
{
  if (OLD_SOLID (cell))
    if (GFS_STATE(cell)->solid) {
      GfsSolidVector * tmp = OLD_SOLID (cell);
      OLD_SOLID (cell) = GFS_STATE(cell)->solid;
      GFS_STATE(cell)->solid = tmp;
      tmp = NULL;
    }
    else {
      GFS_STATE(cell)->solid = OLD_SOLID (cell);
      OLD_SOLID (cell) = NULL;
    }
  else
    if (GFS_STATE(cell)->solid) {
      OLD_SOLID (cell) = GFS_STATE(cell)->solid;
      GFS_STATE(cell)->solid = NULL;
    }
}

static void swap_face_fractions_back (GfsSimulation * sim) 
{
  gfs_domain_cell_traverse (GFS_DOMAIN (sim), FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) swap_fractions_back,
			    GFS_SIMULATION_MOVING (sim)->old_solid);
}

static void moving_divergence_distribution_second_order (GSList * merged, DivergenceData * p)
{
  if (merged->next != NULL && merged->next->data != merged->data) {
    gdouble total_volume = 0., total_div = 0.;
    GfsVariable * old_solid_v = GFS_SIMULATION_MOVING (p->domain)->old_solid;
    GSList * i = merged;

    while (i) {
      FttCell * cell = i->data;
      g_assert (FTT_CELL_IS_LEAF (cell));
      gdouble a = OLD_SOLID (cell) ? OLD_SOLID (cell)->a : 1.;
      total_volume += a*ftt_cell_volume (cell);
      total_div += GFS_VALUE (cell, p->div);
      i = i->next;
    }
    total_div /= total_volume;

    i = merged;
    while (i) {
      FttCell * cell = i->data;
      gdouble a = OLD_SOLID (cell) ? OLD_SOLID (cell)->a : 1.;
      GFS_VALUE (cell, p->div) = total_div*a*ftt_cell_volume (cell);
      i = i->next;
    }
  }
}
