/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001 National Institute of Water and Atmospheric Research
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
/*! \file
 * \brief Quad/Octrees.
 */

#include <stdlib.h>
#include "ftt.h"

#define  FTT_CELL_IS_DESTROYED(c) (((c)->flags & FTT_FLAG_DESTROYED) != 0)

gchar * ftt_direction_name[FTT_NEIGHBORS] = {
  "right", "left", "top", "bottom"
#if (!FTT_2D)
  , "front", "back"
#endif /* FTT_3D */
};

gint ftt_opposite_direction[FTT_NEIGHBORS] =
#if      FTT_2D
  {1, 0, 3, 2};
#else  /* FTT_3D */
  {1, 0, 3, 2, 5, 4};
#endif /* FTT_3D */

typedef struct _FttOct      FttOct;
typedef struct _FttRootCell FttRootCell;

static void oct_new (FttCell * parent,
		     gboolean check_neighbors,
		     FttCellInitFunc init,
		     gpointer data)
{
  FttOct * oct;
  guint n;

  g_assert (parent != NULL);
  g_assert (parent->children == NULL);

  oct = g_malloc0 (sizeof (FttOct));
  oct->level = ftt_cell_level (parent);
  oct->parent = parent;

  ftt_cell_pos (parent, &(oct->pos));
  ftt_cell_neighbors (parent, &(oct->neighbors));

  for (n = 0; n < FTT_CELLS; n++) {
    oct->cell[n].parent = oct;
    oct->cell[n].flags = n;
  }

  if (check_neighbors)
    for (n = 0; n < FTT_NEIGHBORS; n++) {
      FttCell * neighbor = oct->neighbors.c[n];
      
      if (neighbor && ftt_cell_level (neighbor) < oct->level) {
	oct_new (neighbor, check_neighbors, init, data);
	oct->neighbors.c[n] = ftt_cell_neighbor (parent, n);
      }
    }

  g_assert (parent->children == NULL);
  parent->children = oct;

  if (init)
    (* init) (parent, data);
}

/**
 * ftt_cell_new:
 * @init: a #FttCellInitFunc or %NULL.
 * @data: user data to pass to @init.
 *
 * Returns: a new root #FttCell, initialized by calling @init (if not %NULL).
 */
FttCell * ftt_cell_new (FttCellInitFunc init,
			gpointer data)
{
  FttCell * cell;

  cell = g_malloc0 (sizeof (FttRootCell));
  if (init)
    (* init) (cell, data);

  return cell;
}

/**
 * ftt_cell_check:
 * @cell: a #FttCell.
 *
 * Returns: %TRUE if cell is consistent, %FALSE otherwise.
 */
gboolean ftt_cell_check (const FttCell * cell)
{
  FttCellNeighbors neighbor;
  guint i, level;

  g_return_val_if_fail (cell != NULL, FALSE);

  ftt_cell_neighbors (cell, &neighbor);
  level = ftt_cell_level (cell);
  for (i = 0; i < FTT_NEIGHBORS; i++)
    if (neighbor.c[i] && 
	!FTT_CELL_IS_LEAF (neighbor.c[i]) &&
	ftt_cell_level (neighbor.c[i]) == level &&
	neighbor.c[i]->children->neighbors.c[FTT_OPPOSITE_DIRECTION (i)] != cell) {
      g_warning ("ftt_cell_check (%p): neighbor %d = %p: %d/%d",
		 cell, 
		 i,		 
	  neighbor.c[i]->children->neighbors.c[FTT_OPPOSITE_DIRECTION (i)],
		 ftt_cell_level (neighbor.c[i]),
		 ftt_cell_level (neighbor.c[i]->children->neighbors.c[FTT_OPPOSITE_DIRECTION (i)]));
      return FALSE;
    }

  return TRUE;
}

/**
 * ftt_cell_refine_single:
 * @cell: a #FttCell.
 * @init: a #FttCellInitFunc or %NULL.
 * @init_data: user data to pass to @init.
 *
 * Refines @cell and eventually its neighbors to ensure that the
 * neighborhood properties are preserved. The new refined cells
 * created are initialized using @init (if not %NULL).  
 */
void ftt_cell_refine_single (FttCell * cell,
			     FttCellInitFunc init,
			     gpointer init_data)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (FTT_CELL_IS_LEAF (cell));

  oct_new (cell, TRUE, init, init_data);
}

/**
 * ftt_cell_refine:
 * @root: a #FttCell.
 * @refine: a #FttCellRefineFunc.
 * @refine_data: user data to pass to @refine.
 * @init: a #FttCellInitFunc or %NULL.
 * @init_data: user data to pass to @init.
 *
 * Recursively refines the tree starting from @root. Each leaf of the
 * tree is tested for refinement using the @refine function. The new
 * refined cells created are initialized using @init (if not %NULL)
 * and are themselves recursively refined.  
 */
void ftt_cell_refine (FttCell * root,
		      FttCellRefineFunc refine,
		      gpointer refine_data,
		      FttCellInitFunc init,
		      gpointer init_data)
{
  guint n;
  FttOct * oct;

  g_return_if_fail (root != NULL);
  g_return_if_fail (refine != NULL);

  if (FTT_CELL_IS_LEAF (root) && !(* refine) (root, refine_data))
    return;

  if (FTT_CELL_IS_LEAF (root))
    oct_new (root, TRUE, init, init_data);

  g_assert (!FTT_CELL_IS_DESTROYED (root));
  oct = root->children;
  for (n = 0; n < FTT_CELLS; n++)
    if (!FTT_CELL_IS_DESTROYED (&(oct->cell[n])))
      ftt_cell_refine (&(oct->cell[n]), refine, refine_data, init, init_data);
}

/**
 * ftt_cell_draw:
 * @cell: a #FttCell.
 * @fp: a file pointer.
 *
 * Outputs in @fp an OOGL (geomview) representation of @cell.  
 */
void ftt_cell_draw (const FttCell * cell, FILE * fp)
{
  gdouble size;
  FttVector p;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (fp != NULL);

  size = ftt_cell_size (cell)/2.;
  ftt_cell_pos (cell, &p);
  fprintf (fp, 
	   "OFF 8 6 12\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n"
	   "%g %g %g\n",
	   p.x - size, p.y - size, p.z - size,
	   p.x + size, p.y - size, p.z - size,
	   p.x + size, p.y + size, p.z - size,
	   p.x - size, p.y + size, p.z - size,
	   p.x - size, p.y - size, p.z + size,
	   p.x + size, p.y - size, p.z + size,
	   p.x + size, p.y + size, p.z + size,
	   p.x - size, p.y + size, p.z + size);
  fputs ("4 3 2 1 0\n"
	 "4 4 5 6 7\n"
	 "4 2 3 7 6\n"
	 "4 0 1 5 4\n"
	 "4 0 4 7 3\n"
	 "4 1 2 6 5\n",
	 fp);
}

/**
 * ftt_face_draw:
 * @face: a #FttCellFace.
 * @fp: a file pointer.
 *
 * Outputs in @fp an OOGL (geomview) representation of @face.  
 */
void ftt_face_draw (const FttCellFace * face, FILE * fp)
{
  gdouble size;
  FttVector p;
#if FTT_2D
  static FttVector dp[FTT_NEIGHBORS][2] = {
    {{1.,-1.,0.},{1.,1.,0.}},
    {{-1.,1.,0.},{-1.,-1,0.}},
    {{1.,1.,0.},{-1.,1.,0.}},
    {{-1.,-1.,0.},{1.,-1.,0.}}
  };
#else  /* FTT_3D */
  static FttVector dp[FTT_NEIGHBORS][4] = {
    {{1.,-1.,1.},{1.,-1.,-1.},{1.,1.,-1.},{1.,1.,1.}},
    {{-1.,-1.,1.},{-1.,-1.,-1.},{-1.,1.,-1.},{-1.,1.,1.}},
    {{1.,1.,1.},{1.,1.,-1.},{-1,1.,-1.},{-1.,1.,1.}},
    {{1.,-1.,1.},{1.,-1.,-1.},{-1,-1.,-1.},{-1.,-1.,1.}},
    {{1.,-1.,1.},{1.,1.,1.},{-1.,1.,1.},{-1.,-1.,1.}},
    {{1.,-1.,-1.},{1.,1.,-1.},{-1.,1.,-1.},{-1.,-1.,-1.}},
  };
#endif /* FTT_3D */

  g_return_if_fail (face != NULL);
  g_return_if_fail (fp != NULL);

  size = ftt_cell_size (face->cell)/2.;
  ftt_cell_pos (face->cell, &p);
#if FTT_2D
  fprintf (fp, "VECT 1 2 0 2 0 %g %g 0 %g %g 0\n",
	   p.x + dp[face->d][0].x*size, 
	   p.y + dp[face->d][0].y*size,
	   p.x + dp[face->d][1].x*size, 
	   p.y + dp[face->d][1].y*size);
#else /* FTT_3D */
  fprintf (fp, 
	   "OFF 4 1 4 "
	   "%g %g %g "
	   "%g %g %g "
	   "%g %g %g "
	   "%g %g %g "
	   "4 0 1 2 3\n",
	   p.x + dp[face->d][0].x*size,
	   p.y + dp[face->d][0].y*size,
	   p.z + dp[face->d][0].z*size,
	   p.x + dp[face->d][1].x*size,
	   p.y + dp[face->d][1].y*size,
	   p.z + dp[face->d][1].z*size,
	   p.x + dp[face->d][2].x*size,
	   p.y + dp[face->d][2].y*size,
	   p.z + dp[face->d][2].z*size,
	   p.x + dp[face->d][3].x*size,
	   p.y + dp[face->d][3].y*size,
	   p.z + dp[face->d][3].z*size);
#endif /* FTT_3D */
}

static gdouble coords[FTT_CELLS][3] =
#if FTT_2D
 {{-1., 1.,0.},
  { 1., 1.,0.},
  {-1.,-1.,0.},
  { 1.,-1.,0.}};
#else  /* FTT_3D */
 {{-1., 1., 1.},
  { 1., 1., 1.},
  {-1.,-1., 1.},
  { 1.,-1., 1.},
  {-1., 1.,-1.},
  { 1., 1.,-1.},
  {-1.,-1.,-1.},
  { 1.,-1.,-1.}};
#endif /* FTT_3D */

/**
 * ftt_cell_relative_pos:
 * @cell: a #FttCell (not a root cell).
 * @pos: a #FttVector.
 *
 * Fills @pos with the coordinates of the center of @cell relative to
 * the center of its parent cell. The length unit is the size of the
 * parent cell.
 */
void ftt_cell_relative_pos (const FttCell * cell,
			    FttVector * pos)
{
  guint n;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (pos != NULL);
  g_return_if_fail (!FTT_CELL_IS_ROOT (cell));

  n = FTT_CELL_ID (cell);
  pos->x = coords[n][0]/4.;
  pos->y = coords[n][1]/4.;
  pos->z = coords[n][2]/4.;
}

/**
 * ftt_cell_pos:
 * @cell: a #FttCell.
 * @pos: a #FttVector.
 *
 * Fills @pos with the coordinates of the center of @cell.  
 */
void ftt_cell_pos (const FttCell * cell, 
		   FttVector * pos)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (pos != NULL);

  if (FTT_CELL_IS_ROOT (cell))
    *pos = FTT_ROOT_CELL (cell)->pos;
  else {
    gdouble size;
    guint n;

    size = ftt_cell_size (cell)/2.;
    n = FTT_CELL_ID (cell);
    pos->x = cell->parent->pos.x + coords[n][0]*size;
    pos->y = cell->parent->pos.y + coords[n][1]*size;
    pos->z = cell->parent->pos.z + coords[n][2]*size;
  }
}

/**
 * ftt_corner_relative_pos:
 * @cell: a #FttCell.
 * @d: a set of perpendicular directions.
 * @pos: a #FttVector.
 *
 * Fills @pos with the coordinates (normalised by the size of @cell)
 * of the corner of @cell defined by @d relative to the position of
 * the center of @cell.
 */
void ftt_corner_relative_pos (const FttCell * cell,
			      FttDirection d[FTT_DIMENSION],
			      FttVector * pos)
{
  static gdouble coords[FTT_NEIGHBORS][3] =
#if FTT_2D
    {{0.5,0.,0.},{-0.5,0.,0.},{0.,0.5,0.},{0.,-0.5,0.}};
#else  /* FTT_3D */
    {{0.5,0.,0.},{-0.5,0.,0.},{0.,0.5,0.},{0.,-0.5,0.},{0.,0.,0.5},{0.,0.,-0.5}};
#endif /* FTT_3D */

  g_return_if_fail (cell != NULL);
  g_return_if_fail (pos != NULL);

#if FTT_2D
  pos->x = coords[d[0]][0] + coords[d[1]][0];
  pos->y = coords[d[0]][1] + coords[d[1]][1];
  pos->z = 0.;
#else  /* FTT_3D */
  pos->x = coords[d[0]][0] + coords[d[1]][0] + coords[d[2]][0];
  pos->y = coords[d[0]][1] + coords[d[1]][1] + coords[d[2]][1];
  pos->z = coords[d[0]][2] + coords[d[1]][2] + coords[d[2]][2];
#endif /* FTT_3D */
}

/**
 * ftt_corner_pos:
 * @cell: a #FttCell.
 * @d: a set of perpendicular directions.
 * @pos: a #FttVector.
 *
 * Fills @pos with the coordinates of the corner of @cell defined by
 * @d.
 */
void ftt_corner_pos (const FttCell * cell,
		     FttDirection d[FTT_DIMENSION],
		     FttVector * pos)
{
  gdouble size;
  FttVector p;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (pos != NULL);

  ftt_corner_relative_pos (cell, d, pos);
  ftt_cell_pos (cell, &p);
  size = ftt_cell_size (cell);
  pos->x = p.x + size*pos->x;
  pos->y = p.y + size*pos->y;
  pos->z = p.z + size*pos->z;
}

/**
 * ftt_face_pos:
 * @face: a #FttCellFace.
 * @pos: a #FttVector.
 *
 * Fills @pos with the coordinates of the center of @face.
 */
void ftt_face_pos (const FttCellFace * face, FttVector * pos)
{
  gdouble size;
  static gdouble coords[FTT_NEIGHBORS][3] =
#if FTT_2D
  {{1.,0.,0.},{-1.,0.,0.},{0.,1.,0.},{0.,-1.,0.}};
#else  /* FTT_3D */
  {{1.,0.,0.},{-1.,0.,0.},{0.,1.,0.},{0.,-1.,0.},{0.,0.,1.},{0.,0.,-1.}};
#endif /* FTT_3D */

  g_return_if_fail (face != NULL);
  g_return_if_fail (pos != NULL);

  ftt_cell_pos (face->cell, pos);
  size = ftt_cell_size (face->cell)/2.;
  pos->x += size*coords[face->d][0];
  pos->y += size*coords[face->d][1];
  pos->z += size*coords[face->d][2];
}

static void update_children_pos (FttCell * parent)
{
  if (!FTT_CELL_IS_LEAF (parent)) {
    FttOct * oct = parent->children;
    guint n;

    ftt_cell_pos (parent, &(oct->pos));
    for (n = 0; n < FTT_CELLS; n++)
      if (!FTT_CELL_IS_DESTROYED (&(oct->cell[n])))
	update_children_pos (&(oct->cell[n]));
  }
}

/**
 * ftt_cell_set_pos:
 * @root: a #FttCell, root of a cell tree.
 * @pos: a #FttVector.
 *
 * Sets the position of the center of the @root cell of a cell tree to
 * @pos. Updates the positions of its children recursively.  
 */
void ftt_cell_set_pos (FttCell * root,
		       const FttVector * pos)
{
  g_return_if_fail (root != NULL);
  g_return_if_fail (FTT_CELL_IS_ROOT (root));
  g_return_if_fail (pos != NULL);

  FTT_ROOT_CELL (root)->pos = *pos;
  update_children_pos (root);
}

static void update_children_level (FttCell * parent)
{
  if (!FTT_CELL_IS_LEAF (parent)) {
    FttOct * oct = parent->children;
    guint n;

    oct->level = ftt_cell_level (parent);
    for (n = 0; n < FTT_CELLS; n++)
      if (!FTT_CELL_IS_DESTROYED (&(oct->cell[n])))
	update_children_level (&(oct->cell[n]));
  }
}

/**
 * ftt_cell_set_level:
 * @root: a #FttCell, root of a cell tree.
 * @level: the new level.
 *
 * Sets the level of the @root cell of a cell tree to @level.
 * Updates the levels of its children recursively.  
 */
void ftt_cell_set_level (FttCell * root, guint level)
{
  g_return_if_fail (root != NULL);
  g_return_if_fail (FTT_CELL_IS_ROOT (root));

  FTT_ROOT_CELL (root)->level = level;
  update_children_level (root);
}

static void update_neighbor (FttCell * cell,
			     FttDirection d,
			     FttCellInitFunc init,
			     gpointer init_data)
{
  if (!FTT_CELL_IS_LEAF (cell)) {
    FttCell * neighbor = ftt_cell_neighbor_not_cached (cell, d);
    
    if (neighbor) {
      FttOct * oct = cell->children;
      FttCellChildren children;
      guint i, n;

      g_assert (oct->neighbors.c[d] == NULL ||
		oct->neighbors.c[d] == neighbor);
      oct->neighbors.c[d] = neighbor;
      
      if (ftt_cell_level (neighbor) < oct->level) {
	oct_new (neighbor, TRUE, init, init_data);
	oct->neighbors.c[d] = ftt_cell_neighbor (cell, d);
      }
      
      g_assert (ftt_cell_level (oct->neighbors.c[d]) == oct->level);
      n = ftt_cell_children_direction (cell, d, &children);
      for (i = 0; i < n; i++)
	if (children.c[i])
	  update_neighbor (children.c[i], d, init, init_data);
    }
  }
}

/**
 * ftt_cell_set_neighbor:
 * @root: a #FttCell, root of a cell tree.
 * @neighbor: a #FttCell, root of a cell tree.
 * @d: a direction.
 * @init: a #FttCellInitFunc or %NULL.
 * @init_data: user data to pass to @init.
 * 
 * Sets the cell tree defined by @neighbor as the neighbor in
 * direction @d of the cell tree defined by @root.
 *
 * Any new cell created during the process is initialized using the
 * user-defined function @init.
 *
 * Both @root and @neighbor must be the roots of their respective cell
 * trees.  
 */
void ftt_cell_set_neighbor (FttCell * root,
			    FttCell * neighbor,
			    FttDirection d,
			    FttCellInitFunc init,
			    gpointer init_data)
{
  FttDirection od;

  g_return_if_fail (d < FTT_NEIGHBORS);

  g_return_if_fail (root != NULL);
  g_return_if_fail (FTT_CELL_IS_ROOT (root));

  g_return_if_fail (neighbor != NULL);
  g_return_if_fail (FTT_CELL_IS_ROOT (neighbor));

  g_return_if_fail (ftt_cell_level (root) == ftt_cell_level (neighbor));

  g_return_if_fail (FTT_ROOT_CELL (root)->neighbors.c[d] == NULL);
  FTT_ROOT_CELL (root)->neighbors.c[d] = neighbor;
  update_neighbor (root, d, init, init_data);

  od = FTT_OPPOSITE_DIRECTION (d);
  g_return_if_fail (FTT_ROOT_CELL (neighbor)->neighbors.c[od] == NULL);
  FTT_ROOT_CELL (neighbor)->neighbors.c[od] = root;
  update_neighbor (neighbor, od, init, init_data);
}

static void update_neighbor_match (FttCell * cell,
				   FttDirection d,
				   FttCellInitFunc init,
				   gpointer init_data)
{
  if (!FTT_CELL_IS_LEAF (cell)) {
    FttCell * neighbor = ftt_cell_neighbor_not_cached (cell, d);
    
    if (neighbor) {
      FttOct * oct = cell->children;
      FttCellChildren children;
      guint i, n;

      oct->neighbors.c[d] = neighbor;
      
      if (ftt_cell_level (neighbor) < oct->level) {
	oct_new (neighbor, TRUE, init, init_data);
	oct->neighbors.c[d] = ftt_cell_neighbor (cell, d);
      }
      else if (FTT_CELL_IS_LEAF (neighbor))
	oct_new (neighbor, TRUE, init, init_data);
      
      g_assert (ftt_cell_level (oct->neighbors.c[d]) == oct->level);
      n = ftt_cell_children_direction (cell, d, &children);
      for (i = 0; i < n; i++)
	if (children.c[i])
	  update_neighbor_match (children.c[i], d, init, init_data);
    }
  }
  else { /* leaf cell */
    FttCell * neighbor = ftt_cell_neighbor_not_cached (cell, d);
    
    if (neighbor) {
      g_assert (ftt_cell_level (cell) == ftt_cell_level (neighbor));
      if (!FTT_CELL_IS_LEAF (neighbor)) {
	FttCellChildren children;
	guint i, n;

	oct_new (cell, TRUE, init, init_data);
	n = ftt_cell_children_direction (cell, d, &children);
	for (i = 0; i < n; i++)
	  if (children.c[i])
	    update_neighbor_match (children.c[i], d, init, init_data);
      }
    }
  }
}

/**
 * ftt_cell_set_neighbor_match:
 * @root: a #FttCell, root of a cell tree.
 * @neighbor: a #FttCell, root of a cell tree.
 * @d: a direction.
 * @init: a #FttCellInitFunc or %NULL.
 * @init_data: user data to pass to @init.
 * 
 * Sets the cell tree defined by @neighbor as the neighbor in
 * direction @d of the cell tree defined by @root.
 *
 * The boundary between both trees is matched i.e. the type of the
 * face between any pair of cells belonging to each tree is always
 * %FTT_FINE_FINE. Any new cell created during the process is
 * initialized using the user-defined function @init.
 *
 * Both @root and @neighbor must be the roots of their respective cell
 * trees.
 */
void ftt_cell_set_neighbor_match (FttCell * root,
				  FttCell * neighbor,
				  FttDirection d,
				  FttCellInitFunc init,
				  gpointer init_data)
{
  FttDirection od;

  g_return_if_fail (d < FTT_NEIGHBORS);

  g_return_if_fail (root != NULL);
  g_return_if_fail (FTT_CELL_IS_ROOT (root));

  g_return_if_fail (neighbor != NULL);
  g_return_if_fail (FTT_CELL_IS_ROOT (neighbor));

  g_return_if_fail (ftt_cell_level (root) == ftt_cell_level (neighbor));

  FTT_ROOT_CELL (root)->neighbors.c[d] = neighbor;
  update_neighbor_match (root, d, init, init_data);

  od = FTT_OPPOSITE_DIRECTION (d);
  FTT_ROOT_CELL (neighbor)->neighbors.c[od] = root;
  update_neighbor_match (neighbor, od, init, init_data);
}

static void cell_traverse_pre_order_all (FttCell * cell,
					 gint max_depth,
					 FttCellTraverseFunc func,
					 gpointer data)
{
  FttCell * parent;

  if (max_depth >= 0 && ftt_cell_level (cell) > max_depth)
    return;

  parent = ftt_cell_parent (cell);
  (* func) (cell, data);
  /* check that cell has not been deallocated by @func */
  g_assert (parent == NULL || parent->children != NULL);

  if (!FTT_CELL_IS_LEAF (cell)) {
    FttOct * children = cell->children;
    guint n;

    for (n = 0; n < FTT_CELLS; n++) {
      FttCell * c = &(children->cell[n]);

      if (!FTT_CELL_IS_DESTROYED (c))
	cell_traverse_pre_order_all (c, max_depth, func, data);
    }
  }
}

static void cell_traverse_post_order_all (FttCell * cell,
					  gint max_depth,
					  FttCellTraverseFunc func,
					  gpointer data)
{
  if (max_depth >= 0 && ftt_cell_level (cell) > max_depth)
    return;

  if (!FTT_CELL_IS_LEAF (cell)) {
    FttOct * children = cell->children;
    guint n;

    for (n = 0; n < FTT_CELLS; n++) {
      FttCell * c = &(children->cell[n]);

      if (!FTT_CELL_IS_DESTROYED (c))
	cell_traverse_post_order_all (c, max_depth, func, data);
    }
  }

  (* func) (cell, data);
}

static void cell_traverse_leafs (FttCell * cell,
				 gint max_depth,
				 FttCellTraverseFunc func,
				 gpointer data)
{
  if (max_depth >= 0 && ftt_cell_level (cell) > max_depth)
    return;

  if (FTT_CELL_IS_LEAF (cell))
    (* func) (cell, data);
  if (!FTT_CELL_IS_LEAF (cell)) {
    FttOct * children = cell->children;
    guint n;

    for (n = 0; n < FTT_CELLS; n++) {
      FttCell * c = &(children->cell[n]);

      if (!FTT_CELL_IS_DESTROYED (c))
	cell_traverse_leafs (c, max_depth, func, data);
    }
  }
}

static void cell_traverse_pre_order_nonleafs (FttCell * cell,
					      gint max_depth,
					      FttCellTraverseFunc func,
					      gpointer data)
{
  if (max_depth >= 0 && ftt_cell_level (cell) > max_depth)
    return;

  if (!FTT_CELL_IS_LEAF (cell)) {
    FttCell * parent = ftt_cell_parent (cell);

    (* func) (cell, data);
    /* check that cell has not been deallocated by @func */
    g_assert (parent == NULL || parent->children != NULL);
    if (!FTT_CELL_IS_LEAF (cell)) {
      FttOct * children = cell->children;
      guint n;

      for (n = 0; n < FTT_CELLS; n++) {
	FttCell * c = &(children->cell[n]);
	
	if (!FTT_CELL_IS_DESTROYED (c))
	  cell_traverse_pre_order_nonleafs (c, max_depth, func, data);
      }
    }
  }
}

static void cell_traverse_post_order_nonleafs (FttCell * cell,
					       gint max_depth,
					       FttCellTraverseFunc func,
					       gpointer data)
{
  if (max_depth >= 0 && ftt_cell_level (cell) > max_depth)
    return;

  if (!FTT_CELL_IS_LEAF (cell)) {
    FttOct * children = cell->children;
    guint n;

    for (n = 0; n < FTT_CELLS; n++) {
      FttCell * c = &(children->cell[n]);

      if (!FTT_CELL_IS_DESTROYED (c))
	cell_traverse_post_order_nonleafs (c, max_depth, func, data);
    }

    (* func) (cell, data);
  }
}

static void cell_traverse_level (FttCell * cell,
				 gint max_depth,
				 FttCellTraverseFunc func,
				 gpointer data)
{
  if (ftt_cell_level (cell) == max_depth)
    (* func) (cell, data);
  else if (!FTT_CELL_IS_LEAF (cell)) {
    FttOct * children = cell->children;
    guint n;

    for (n = 0; n < FTT_CELLS; n++) {
      FttCell * c = &(children->cell[n]);

      if (!FTT_CELL_IS_DESTROYED (c))
	cell_traverse_level (c, max_depth, func, data);
    }
  }
}

static void cell_traverse_level_leafs (FttCell * cell,
				       gint max_depth,
				       FttCellTraverseFunc func,
				       gpointer data)
{
  if (ftt_cell_level (cell) == max_depth || FTT_CELL_IS_LEAF (cell))
    (* func) (cell, data);
  else if (!FTT_CELL_IS_LEAF (cell)) {
    FttOct * children = cell->children;
    guint n;

    for (n = 0; n < FTT_CELLS; n++) {
      FttCell * c = &(children->cell[n]);

      if (!FTT_CELL_IS_DESTROYED (c))
	cell_traverse_level_leafs (c, max_depth, func, data);
    }
  }
}

static void cell_traverse_level_non_leafs (FttCell * cell,
					   gint max_depth,
					   FttCellTraverseFunc func,
					   gpointer data)
{
  if (ftt_cell_level (cell) == max_depth && !FTT_CELL_IS_LEAF (cell))
    (* func) (cell, data);
  else if (!FTT_CELL_IS_LEAF (cell)) {
    FttOct * children = cell->children;
    guint n;

    for (n = 0; n < FTT_CELLS; n++) {
      FttCell * c = &(children->cell[n]);

      if (!FTT_CELL_IS_DESTROYED (c))
	cell_traverse_level_non_leafs (c, max_depth, func, data);
    }
  }
}

/**
 * ftt_cell_traverse:
 * @root: the root #FttCell of the tree to traverse.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @max_depth: the maximum depth of the traversal. Cells below this
 * depth will not be traversed. If @max_depth is -1 all cells in the
 * tree are visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 *
 * Traverses a cell tree starting at the given root #FttCell. Calls
 * the given function for each cell visited.  
 */
void ftt_cell_traverse (FttCell * root,
			FttTraverseType order,
			FttTraverseFlags flags,
			gint max_depth,
			FttCellTraverseFunc func,
			gpointer data)
{
  g_return_if_fail (root != NULL);
  g_return_if_fail (func != NULL);

  if (max_depth >= 0 && ftt_cell_level (root) > max_depth)
    return;

  if (flags == FTT_TRAVERSE_ALL) {
    if (order == FTT_PRE_ORDER)
      cell_traverse_pre_order_all (root, max_depth, func, data);
    else
      cell_traverse_post_order_all (root, max_depth, func, data);
  }
  else if ((flags & FTT_TRAVERSE_LEVEL) != 0) {
    if ((flags & FTT_TRAVERSE_LEAFS) != 0)
      cell_traverse_level_leafs (root, max_depth, func, data);
    else if ((flags & FTT_TRAVERSE_NON_LEAFS) != 0)
      cell_traverse_level_non_leafs (root, max_depth, func, data);
    else
      cell_traverse_level (root, max_depth, func, data);
  }
  else if ((flags & FTT_TRAVERSE_LEAFS) != 0)
    cell_traverse_leafs (root, max_depth, func, data);
  else {
    g_return_if_fail ((flags & FTT_TRAVERSE_NON_LEAFS) != 0);

    if (order == FTT_PRE_ORDER)
      cell_traverse_pre_order_nonleafs (root, max_depth, func, data);
    else
      cell_traverse_post_order_nonleafs (root, max_depth, func, data);
  }
}

/**
 * ftt_cell_traverse_condition:
 * @root: the root #FttCell of the tree to traverse.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @max_depth: the maximum depth of the traversal. Cells below this
 * depth will not be traversed. If @max_depth is -1 all cells in the
 * tree are visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 * @condition: the condition.
 * @cdata: user data to pass to @condition.
 *
 * Traverses a cell tree starting at the given root #FttCell. Calls
 * the given function for each cell visited.
 *
 * Traversal of any branch of the tree is stopped whenever @condition
 * is not verified.
 */
void ftt_cell_traverse_condition (FttCell * root,
				  FttTraverseType order,
				  FttTraverseFlags flags,
				  gint max_depth,
				  FttCellTraverseFunc func,
				  gpointer data,
				  gboolean (* condition) (FttCell *, gpointer),
				  gpointer cdata)
{
  g_return_if_fail (root != NULL);
  g_return_if_fail (func != NULL);
  g_return_if_fail (condition != NULL);

  if ((max_depth >= 0 && ftt_cell_level (root) > max_depth) ||
      !(* condition) (root, cdata))
    return;

  if (order == FTT_PRE_ORDER &&
      (flags == FTT_TRAVERSE_ALL ||
       ((flags & FTT_TRAVERSE_LEAFS) != 0 && FTT_CELL_IS_LEAF (root)) ||
       ((flags & FTT_TRAVERSE_NON_LEAFS) != 0 && !FTT_CELL_IS_LEAF (root))))
    (* func) (root, data);
  if (!FTT_CELL_IS_LEAF (root)) {
    struct _FttOct * children = root->children;
    guint n;

    for (n = 0; n < FTT_CELLS; n++) {
      FttCell * c = &(children->cell[n]);

      if (!FTT_CELL_IS_DESTROYED (c))
	ftt_cell_traverse_condition (c, order, flags, max_depth, func, data, condition, cdata);
    }
  }
  if (order == FTT_POST_ORDER &&
      (flags == FTT_TRAVERSE_ALL ||
       ((flags & FTT_TRAVERSE_LEAFS) != 0 && FTT_CELL_IS_LEAF (root)) ||
       ((flags & FTT_TRAVERSE_NON_LEAFS) != 0 && !FTT_CELL_IS_LEAF (root))))
    (* func) (root, data);
}

/**
 * ftt_cell_bbox:
 * @cell: a #FttCell.
 * @bb: a #GtsBBox.
 *
 * Fills @bb with the bounding box of @cell.
 */
void ftt_cell_bbox (const FttCell * cell, GtsBBox * bb)
{
  FttVector p;
  gdouble h;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (bb != NULL);
  
  h = ftt_cell_size (cell)/1.99999;
  ftt_cell_pos (cell, &p);
  bb->x1 = p.x - h; bb->y1 = p.y - h;
  bb->x2 = p.x + h; bb->y2 = p.y + h; 
#if FTT_2D
  bb->z1 = bb->z2 = 0.;
#else  /* 3D */
  bb->z1 = p.z - h; bb->z2 = p.z + h;
#endif /* 3D */
}

static gboolean cell_is_in_box (FttCell * cell, gpointer data)
{
  GtsBBox * box = data;
  GtsBBox bb;

  ftt_cell_bbox (cell, &bb);
  return gts_bboxes_are_overlapping (&bb, box);
}

/**
 * ftt_cell_traverse_box:
 * @root: the root #FttCell of the tree to traverse.
 * @box: a #GtsBBox.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @max_depth: the maximum depth of the traversal. Cells below this
 * depth will not be traversed. If @max_depth is -1 all cells in the
 * tree are visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 *
 * Traverses a cell tree starting at the given root #FttCell. Calls
 * the given function for each cell visited. Only the cells partly or
 * totally contained within @box are visited.  
 */
void ftt_cell_traverse_box (FttCell * root,
			    GtsBBox * box,
			    FttTraverseType order,
			    FttTraverseFlags flags,
			    gint max_depth,
			    FttCellTraverseFunc func,
			    gpointer data)
{
  g_return_if_fail (root != NULL);
  g_return_if_fail (box != NULL);
  g_return_if_fail (func != NULL);

  ftt_cell_traverse_condition (root, order, flags, max_depth, func, data, cell_is_in_box, box);
}

static void cell_traverse_boundary_pre_order_all (FttCell * cell,
						  FttDirection d,
						  gint max_depth,
						  FttCellTraverseFunc func,
						  gpointer data)
{
  FttCell * parent;

  if (max_depth >= 0 && ftt_cell_level (cell) > max_depth)
    return;

  parent = ftt_cell_parent (cell);
  (* func) (cell, data);
  /* check that cell has not been deallocated by @func */
  g_assert (parent == NULL || parent->children != NULL);

  if (!FTT_CELL_IS_LEAF (cell)) {
    FttCellChildren children;
    guint i, n;

    n = ftt_cell_children_direction (cell, d, &children);
    for (i = 0; i < n; i++)
      if (children.c[i])
	cell_traverse_boundary_pre_order_all (children.c[i], d, 
					      max_depth, func, data);
  }
}

static void cell_traverse_boundary_post_order_all (FttCell * cell,
						   FttDirection d,
						   gint max_depth,
						   FttCellTraverseFunc func,
						   gpointer data)
{
  if (max_depth >= 0 && ftt_cell_level (cell) > max_depth)
    return;

  if (!FTT_CELL_IS_LEAF (cell)) {
    FttCellChildren children;
    guint i, n;

    n = ftt_cell_children_direction (cell, d, &children);
    for (i = 0; i < n; i++)
      if (children.c[i])
	cell_traverse_boundary_post_order_all (children.c[i], d, 
					       max_depth, func, data);
  }

  (* func) (cell, data);
}

static void cell_traverse_boundary_leafs (FttCell * cell,
					  FttDirection d,
					  gint max_depth,
					  FttCellTraverseFunc func,
					  gpointer data)
{
  if (max_depth >= 0 && ftt_cell_level (cell) > max_depth)
    return;

  if (FTT_CELL_IS_LEAF (cell))
    (* func) (cell, data);
  else {
    FttCellChildren children;
    guint i, n;

    n = ftt_cell_children_direction (cell, d, &children);
    for (i = 0; i < n; i++)
      if (children.c[i])
	cell_traverse_boundary_leafs (children.c[i], d, 
				      max_depth, func, data);
  }
}

static void cell_traverse_boundary_pre_order_nonleafs (FttCell * cell,
						       FttDirection d,
						       gint max_depth,
				   FttCellTraverseFunc func,
						       gpointer data)
{
  if (max_depth >= 0 && ftt_cell_level (cell) > max_depth)
    return;

  if (!FTT_CELL_IS_LEAF (cell)) {
    FttCell * parent = ftt_cell_parent (cell);

    (* func) (cell, data);
    /* check that cell has not been deallocated by @func */
    g_assert (parent == NULL || parent->children != NULL);
    if (!FTT_CELL_IS_LEAF (cell)) {
      FttCellChildren children;
      guint i, n;

      n = ftt_cell_children_direction (cell, d, &children);
      for (i = 0; i < n; i++)
	if (children.c[i])
	  cell_traverse_boundary_pre_order_nonleafs (children.c[i], d, 
						     max_depth, func, data);
    }
  }
}

static void cell_traverse_boundary_post_order_nonleafs (FttCell * cell,
							FttDirection d,
							gint max_depth,
				    FttCellTraverseFunc func,
							gpointer data)
{
  if (max_depth >= 0 && ftt_cell_level (cell) > max_depth)
    return;

  if (!FTT_CELL_IS_LEAF (cell)) {
    FttCellChildren children;
    guint i, n;

    n = ftt_cell_children_direction (cell, d, &children);
    for (i = 0; i < n; i++)
      if (children.c[i])
	cell_traverse_boundary_post_order_nonleafs (children.c[i], d, 
						    max_depth, func, data);
    (* func) (cell, data);
  }
}

static void cell_traverse_boundary_level (FttCell * cell,
					  FttDirection d,
					  gint max_depth,
					  FttCellTraverseFunc func,
					  gpointer data)
{
  if (ftt_cell_level (cell) == max_depth)
    (* func) (cell, data);
  else if (!FTT_CELL_IS_LEAF (cell)) {
    FttCellChildren children;
    guint i, n;

    n = ftt_cell_children_direction (cell, d, &children);
    for (i = 0; i < n; i++)
      if (children.c[i])
	cell_traverse_boundary_level (children.c[i], d, 
				      max_depth, func, data);
  }
}

static void cell_traverse_boundary_level_leafs (FttCell * cell,
						FttDirection d,
						gint max_depth,
						FttCellTraverseFunc func,
						gpointer data)
{
  if (ftt_cell_level (cell) == max_depth || FTT_CELL_IS_LEAF (cell))
    (* func) (cell, data);
  else if (!FTT_CELL_IS_LEAF (cell)) {
    FttCellChildren children;
    guint i, n;

    n = ftt_cell_children_direction (cell, d, &children);
    for (i = 0; i < n; i++)
      if (children.c[i])
	cell_traverse_boundary_level_leafs (children.c[i], d, 
					    max_depth, func, data);
  }
}

static void cell_traverse_boundary_level_non_leafs (FttCell * cell,
						    FttDirection d,
						    gint max_depth,
						    FttCellTraverseFunc func,
						    gpointer data)
{
  if (ftt_cell_level (cell) == max_depth && !FTT_CELL_IS_LEAF (cell))
    (* func) (cell, data);
  else if (!FTT_CELL_IS_LEAF (cell)) {
    FttCellChildren children;
    guint i, n;

    n = ftt_cell_children_direction (cell, d, &children);
    for (i = 0; i < n; i++)
      if (children.c[i])
	cell_traverse_boundary_level_non_leafs (children.c[i], d, 
						max_depth, func, data);
  }
}

/**
 * ftt_cell_traverse_boundary:
 * @root: the root #FttCell of the tree to traverse.
 * @d: the direction of the boundary to traverse.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @max_depth: the maximum depth of the traversal. Cells below this
 * depth will not be traversed. If @max_depth is -1 all cells in the
 * tree are visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 *
 * Traverses the boundary of a cell tree in direction @d starting at
 * the given root #FttCell. Calls the given function for each node
 * visited.  
 */
void ftt_cell_traverse_boundary (FttCell * root,
				 FttDirection d,
				 FttTraverseType order,
				 FttTraverseFlags flags,
				 gint max_depth,
				 FttCellTraverseFunc func,
				 gpointer data)
{
  g_return_if_fail (root != NULL);
  g_return_if_fail (d < FTT_NEIGHBORS);
  g_return_if_fail (func != NULL);

  if (max_depth >= 0 && ftt_cell_level (root) > max_depth)
    return;

  if (flags == FTT_TRAVERSE_ALL) {
    if (order == FTT_PRE_ORDER)
      cell_traverse_boundary_pre_order_all (root, d, max_depth, func, data);
    else
      cell_traverse_boundary_post_order_all (root, d, max_depth, func, data);
  }
  else if ((flags & FTT_TRAVERSE_LEVEL) != 0) {
    if ((flags & FTT_TRAVERSE_LEAFS) != 0)
      cell_traverse_boundary_level_leafs (root, d, max_depth, func, data);
    else if ((flags & FTT_TRAVERSE_NON_LEAFS) != 0)
      cell_traverse_boundary_level_non_leafs (root, d, max_depth, func, data);
    else
      cell_traverse_boundary_level (root, d, max_depth, func, data);
  }
  else if ((flags & FTT_TRAVERSE_LEAFS) != 0)
    cell_traverse_boundary_leafs (root, d, max_depth, func, data);
  else {
    g_return_if_fail ((flags & FTT_TRAVERSE_NON_LEAFS) != 0);

    if (order == FTT_PRE_ORDER)
      cell_traverse_boundary_pre_order_nonleafs (root, d, 
						 max_depth, func, data);
    else
      cell_traverse_boundary_post_order_nonleafs (root, d, 
						  max_depth, func, data);
  }
}

static void oct_destroy (FttOct * oct,
			 FttCellCleanupFunc cleanup,
			 gpointer data)
{
  guint n;

  g_return_if_fail (oct != NULL);
  g_return_if_fail (oct->parent->children == oct);

  oct->parent->children = NULL;
  for (n = 0; n < FTT_CELLS; n++)
    ftt_cell_destroy (&(oct->cell[n]), cleanup, data);
  g_free (oct);
}

/**
 * ftt_cell_destroy:
 * @cell: a #FttCell.
 * @cleanup: a #FttCellCleanupFunc to call before destroying @cell or %NULL.
 * @data: user data to pass to @cleanup.
 *
 * Frees all memory allocated for @cell and its descendants.
 *
 * The user-defined function @cleanup is called prior to freeing memory.
 */
void ftt_cell_destroy (FttCell * cell,
		       FttCellCleanupFunc cleanup,
		       gpointer data)
{
  FttCellNeighbors neighbor;
  guint i, level;

  g_return_if_fail (cell != NULL);

  if (FTT_CELL_IS_DESTROYED (cell))
    return;

  ftt_cell_neighbors (cell, &neighbor);
  level = ftt_cell_level (cell);

  if (cleanup)
    (* cleanup) (cell, data);
  cell->flags |= FTT_FLAG_DESTROYED;

  /* destroy children */
  if (!FTT_CELL_IS_LEAF (cell)) {
    oct_destroy (cell->children, cleanup, data);
    cell->children = NULL;
  }

  /* update relationships for neighbors */
  for (i = 0; i < FTT_NEIGHBORS; i++)
    if (neighbor.c[i] && ftt_cell_level (neighbor.c[i]) == level) {
      FttDirection od = FTT_OPPOSITE_DIRECTION (i);

      if (FTT_CELL_IS_ROOT (neighbor.c[i])) {
	FttCell * opneighbor = FTT_ROOT_CELL (neighbor.c[i])->neighbors.c[od];

	g_assert (opneighbor == cell);
	FTT_ROOT_CELL (neighbor.c[i])->neighbors.c[od] = NULL;
      }
      if (!FTT_CELL_IS_LEAF (neighbor.c[i]))
	neighbor.c[i]->children->neighbors.c[od] = NULL;
    }
  
  if (FTT_CELL_IS_ROOT (cell))
    g_free (cell);
  else if (!FTT_CELL_IS_LEAF (cell->parent->parent)) {
    /* if parent Oct is not already destroyed and empty destroy it */
    FttOct * parent = cell->parent;
    gboolean empty = TRUE;

    for (i = 0; i < FTT_CELLS && empty; i++)
      if (!FTT_CELL_IS_DESTROYED (&(parent->cell[i])))
	empty = FALSE;
    if (empty)
      oct_destroy (parent, NULL, NULL);
  }
}

/**
 * ftt_cell_destroy_root:
 * @root: the root cell of a cell tree.
 * @children: a #FttCellChildren.
 * @cleanup: a #FttCellCleanupFunc to call before destroying a cell.
 * @data: user data to pass to @cleanup.
 *
 * Destroys the root cell of a cell tree but not its children. Each
 * child becomes the root cell of a new cell tree. The new (orphaned)
 * children are returned in @children.
 *
 * Note that the function will fail if @root is also a leaf cell.
 */
void ftt_cell_destroy_root (FttCell * root,
			    FttCellChildren * children,
			    FttCellCleanupFunc cleanup,
			    gpointer data)
{
  guint i;
  FttCellNeighbors neighbor;
  FttCellChildren child;

  g_return_if_fail (root != NULL);
  g_return_if_fail (FTT_CELL_IS_ROOT (root));
  g_return_if_fail (!FTT_CELL_IS_LEAF (root));
  g_return_if_fail (!FTT_CELL_IS_DESTROYED (root));
  g_return_if_fail (children != NULL);

  if (cleanup)
    (* cleanup) (root, data);
  root->flags |= FTT_FLAG_DESTROYED;

  ftt_cell_neighbors (root, &neighbor);
  for (i = 0; i < FTT_NEIGHBORS; i++)
    if (neighbor.c[i]) {
      FttDirection od = FTT_OPPOSITE_DIRECTION (i);
      
      g_assert (FTT_CELL_IS_ROOT (neighbor.c[i]));
      g_assert (FTT_ROOT_CELL (neighbor.c[i])->neighbors.c[od] == root);
      FTT_ROOT_CELL (neighbor.c[i])->neighbors.c[od] = NULL;

      if (!FTT_CELL_IS_LEAF (neighbor.c[i]))
	neighbor.c[i]->children->neighbors.c[od] = NULL;
    }

  ftt_cell_children (root, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i]) {
      FttCell * newc;
      FttDirection d;

      newc = g_malloc0 (sizeof (FttRootCell));
      newc->data = child.c[i]->data;
      newc->children = child.c[i]->children;
      ftt_cell_pos (child.c[i], &FTT_ROOT_CELL (newc)->pos);
      FTT_ROOT_CELL (newc)->level = ftt_cell_level (child.c[i]);
      ftt_cell_neighbors (child.c[i], &FTT_ROOT_CELL (newc)->neighbors);
      g_return_if_fail (!FTT_CELL_IS_LEAF (newc));
      newc->children->parent = newc;
      children->c[i] = newc;

      neighbor = FTT_ROOT_CELL (newc)->neighbors;
      for (d = 0; d < FTT_NEIGHBORS; d++)
	if (neighbor.c[d]) {
	  FttDirection od = FTT_OPPOSITE_DIRECTION (d);

	  if (FTT_CELL_IS_ROOT (neighbor.c[d])) {
	    g_assert (FTT_ROOT_CELL (neighbor.c[d])->neighbors.c[od] 
		      == child.c[i]);
	    FTT_ROOT_CELL (neighbor.c[d])->neighbors.c[od] = newc;
	  }
	  if (!FTT_CELL_IS_LEAF (neighbor.c[d])) {
	    g_assert (neighbor.c[d]->children->neighbors.c[od] == child.c[i]);
	    neighbor.c[d]->children->neighbors.c[od] = newc;
	  }
	}
    }
    else
      children->c[i] = NULL;

  g_free (root->children);
  g_free (root);
}

/**
 * ftt_cell_flatten:
 * @root: the root of the cell tree to flatten.
 * @d: the direction in which to flatten.
 * @cleanup: a #FttCellCleanupFunc to call before destroying a cell.
 * @data: user data to pass to @cleanup.
 *
 * Recursively destroys all the cells of the tree defined by @root
 * which do not form the boundary in direction @d. The resulting cell
 * tree is in effect a domain "flattened" in direction @d.
 *
 * The resulting domain is always one-cell thick in direction @d.  
 */
void ftt_cell_flatten (FttCell * root,
		       FttDirection d,
		       FttCellCleanupFunc cleanup,
		       gpointer data)
{
  g_return_if_fail (root != NULL);
  g_return_if_fail (d < FTT_NEIGHBORS);

  if (!FTT_CELL_IS_LEAF (root)) {
    struct _FttOct * oct;
    guint i;
#if FTT_2D
    static gint index[FTT_NEIGHBORS_2D][FTT_CELLS/2] =
    {{1, 3},
     {0, 2},
     {0, 1},
     {2, 3}};
#else  /* FTT_3D */
    static gint index[FTT_NEIGHBORS][FTT_CELLS/2] =
    {{1, 3, 5, 7},
     {0, 2, 4, 6},
     {0, 1, 4, 5},
     {2, 3, 6, 7},
     {0, 1, 2, 3},
     {4, 5, 6, 7}};
#endif /* FTT_3D */
    FttDirection od = FTT_OPPOSITE_DIRECTION (d);

    oct = root->children;
    for (i = 0; i < FTT_CELLS/2; i++) {
      FttCell * c = &(oct->cell[index[od][i]]);
      if (!FTT_CELL_IS_DESTROYED (c))
	ftt_cell_destroy (c, cleanup, data);
    }
    if (!FTT_CELL_IS_LEAF (root))
      for (i = 0; i < FTT_CELLS/2; i++)
	if (!FTT_CELL_IS_DESTROYED (&(oct->cell[index[d][i]])))
	  ftt_cell_flatten (&(oct->cell[index[d][i]]), d, cleanup, data);
  }
}

/**
 * ftt_cell_locate:
 * @root: a #FttCell.
 * @target: position of the point to look for.
 * @max_depth: maximum depth to consider (-1 means no restriction, see below for -2).
 *
 * Locates the cell of the tree defined by @root containing
 * @target. This is done efficiently in log(n) operations by using the
 * topology of the tree.
 *
 * If @max_depth is set to -2, the finest cell containing @target is
 * returned. This cell is not necessarily a leaf-cell in contrast to
 * the case where @max_depth is set to -1.
 *
 * Returns: a #FttCell of the tree defined by @root and
 * containing (boundary included) the point defined by @target or
 * %NULL if @target is not contained in any cell of @root.  
 */
FttCell * ftt_cell_locate (FttCell * root, 
			   FttVector target,
			   gint max_depth)
{
  FttVector pos;
  gdouble size;

  g_return_val_if_fail (root != NULL, NULL);

  ftt_cell_pos (root, &pos);
  size = ftt_cell_size (root)/2.;

  if (target.x > pos.x + size || target.x < pos.x - size ||
      target.y > pos.y + size || target.y < pos.y - size
#if !FTT_2D
      || target.z > pos.z + size || target.z < pos.z - size
#endif
      )
    return NULL;

  do {
    if (FTT_CELL_IS_LEAF (root) || ftt_cell_level (root) == max_depth)
      return root;
#if FTT_2D
    static guint index[2][2] = {{2,3},{0,1}};
    guint n = index[target.y > pos.y][target.x > pos.x];
#else  /* 3D */
    static guint index[2][2][2] = {{{6,7},{4,5}},{{2,3},{0,1}}};
    guint n = index[target.z > pos.z][target.y > pos.y][target.x > pos.x];
#endif /* 3D */
    root = &(root->children->cell[n]);
    size /= 2.;
    pos.x += coords[n][0]*size;
    pos.y += coords[n][1]*size;
#if !FTT_2D
    pos.z += coords[n][2]*size;
#endif /* 3D */
  } while (!FTT_CELL_IS_DESTROYED (root));
  return max_depth == -2 ? ftt_cell_parent (root) : NULL;
}

static void bubble_sort (FttCellChildren * child, gdouble * d)
{
  guint i, j;

  for (i = 0; i < FTT_CELLS - 1; i++)
    for (j = 0; j < FTT_CELLS - 1 - i; j++)
      if (d[j+1] < d[j]) {
	gdouble tmp = d[j];
	FttCell * cell = child->c[j];
	d[j] = d[j+1];
	d[j+1] = tmp;
	child->c[j] = child->c[j+1];
	child->c[j+1] = cell;
      }
}

/**
 * ftt_cell_point_distance2_min:
 * @cell: a #FttCell.
 * @p: a #GtsPoint.
 * 
 * Returns: the square of the minimum distance between @cell and @p.
 */
gdouble ftt_cell_point_distance2_min (FttCell * cell, GtsPoint * p)
{
  GtsBBox bb;
  gdouble dmin, xd1, xd2, yd1, yd2, zd1, zd2;
    
  g_return_val_if_fail (cell != NULL, G_MAXDOUBLE);
  g_return_val_if_fail (p != NULL, G_MAXDOUBLE);

  ftt_cell_bbox (cell, &bb);

  xd1 = (bb.x1 - p->x)*(bb.x1 - p->x);
  xd2 = (p->x - bb.x2)*(p->x - bb.x2);
  yd1 = (bb.y1 - p->y)*(bb.y1 - p->y);
  yd2 = (p->y - bb.y2)*(p->y - bb.y2);
  zd1 = (bb.z1 - p->z)*(bb.z1 - p->z);
  zd2 = (p->z - bb.z2)*(p->z - bb.z2);
  
  dmin = p->x < bb.x1 ? xd1 : p->x > bb.x2 ? xd2 : 0.0;
  dmin += p->y < bb.y1 ? yd1 : p->y > bb.y2 ? yd2 : 0.0;
  dmin += p->z < bb.z1 ? zd1 : p->z > bb.z2 ? zd2 : 0.0;

  return dmin;
}

void ftt_cell_point_distance2_internal (FttCell * root,
					GtsPoint * p,
					gdouble d,
					gdouble (* distance2) (FttCell *, GtsPoint *, gpointer),
					gpointer data,
					FttCell ** closest,
					gdouble * dmin)
{
  if (FTT_CELL_IS_LEAF (root)) {
    if (d < *dmin) {
      *dmin = d;
      if (closest)
	*closest = root;
    }
  }
  else {
    FttCellChildren child;
    gdouble dc[FTT_CELLS];
    guint i;

    ftt_cell_children (root, &child);
    for (i = 0; i < FTT_CELLS; i++)
      dc[i] = child.c[i] ? (* distance2) (child.c[i], p, data) : G_MAXDOUBLE;
    bubble_sort (&child, dc);
    for (i = 0; i < FTT_CELLS; i++)
      if (dc[i] < *dmin)
	ftt_cell_point_distance2_internal (child.c[i], p, dc[i], distance2, data, closest, dmin);
  }
}

/**
 * ftt_cell_point_distance2:
 * @root: a #FttCell.
 * @p: a #GtsPoint.
 * @distance2: the squared distance function.
 * @data: user data to pass to @distance2.
 * @closest: where to return the closest cell or %NULL.
 *
 * For non-leafs cells @distance2 must return a lower-bound for the
 * minimum distance (using for example ftt_cell_point_distance2_min()).
 *
 * Returns: the square of the minimum distance measured according to
 * @distance2 between @p and a leaf cell of @root.
 */
gdouble ftt_cell_point_distance2 (FttCell * root,
				  GtsPoint * p,
				  gdouble (* distance2) (FttCell *, GtsPoint *, gpointer),
				  gpointer data,
				  FttCell ** closest)
{
  gdouble d, dmin = G_MAXDOUBLE;

  g_return_val_if_fail (root != NULL, dmin);
  g_return_val_if_fail (p != NULL, dmin);
  g_return_val_if_fail (distance2 != NULL, dmin);

  if (closest)
    *closest = NULL;
  d = (* distance2) (root, p, data);
  if (d < dmin)
    ftt_cell_point_distance2_internal (root, p, d, distance2, data, closest, &dmin);
  return dmin;
}

/**
 * ftt_cell_depth:
 * @root: a #FttCell.
 *
 * Returns: the depth of the tree starting at @root, i.e. the maximum
 * level of any cell descendant of @root.  
 */
guint ftt_cell_depth (const FttCell * root)
{
  guint depth;

  g_return_val_if_fail (root != NULL, 0);

  depth = ftt_cell_level (root);
  if (root->children) {
    FttOct * oct = root->children;
    guint n;
    
    for (n = 0; n < FTT_CELLS; n++) 
      if (!FTT_CELL_IS_DESTROYED (&(oct->cell[n]))) {
	guint d = ftt_cell_depth (&(oct->cell[n]));
	if (d > depth)
	  depth = d;
      }
  }
  return depth;
}

/**
 * ftt_cell_write:
 * @root: a #FttCell.
 * @max_depth: the maximum depth at which to stop writing (-1 means no limit).
 * @fp: a file pointer.
 * @write: a #FttCellWriteFunc function or %NULL.
 * @data: user data to pass to @write.
 *
 * Writes in the file pointed to by @fp a text representation of the
 * cell tree starting at @root. If not %NULL, the user-defined
 * function @write is used to write the extra user data associated
 * with each cell.  
 */
void ftt_cell_write (const FttCell * root,
		     gint max_depth,
		     FILE * fp,
		     FttCellWriteFunc write,
		     gpointer data)
{
  guint flags;

  g_return_if_fail (root != NULL);
  g_return_if_fail (fp != NULL);

  flags = root->flags;
  if (FTT_CELL_IS_LEAF (root) || ftt_cell_level (root) == max_depth)
    flags |= FTT_FLAG_LEAF;

  fprintf (fp, "%u", flags);
  if (write && !FTT_CELL_IS_DESTROYED (root))
    (* write) (root, fp, data);
  fputc ('\n', fp);

  if ((flags & FTT_FLAG_LEAF) == 0) {
    FttOct * oct;
    guint i;

    oct = root->children;
    for (i = 0; i < FTT_CELLS; i++)
      ftt_cell_write (&(oct->cell[i]), max_depth, fp, write, data);
  }
}

/**
 * ftt_cell_write_binary:
 * @root: a #FttCell.
 * @max_depth: the maximum depth at which to stop writing (-1 means no limit).
 * @fp: a file pointer.
 * @write: a #FttCellWriteFunc function or %NULL.
 * @data: user data to pass to @write.
 *
 * Writes in the file pointed to by @fp a binary representation of the
 * cell tree starting at @root. If not %NULL, the user-defined
 * function @write is used to write the extra user data associated
 * with each cell.  
 */
void ftt_cell_write_binary (const FttCell * root,
			    gint max_depth,
			    FILE * fp,
			    FttCellWriteFunc write,
			    gpointer data)
{
  guint flags;

  g_return_if_fail (root != NULL);
  g_return_if_fail (fp != NULL);

  flags = root->flags;
  if (FTT_CELL_IS_LEAF (root) || ftt_cell_level (root) == max_depth)
    flags |= FTT_FLAG_LEAF;

  fwrite (&flags, sizeof (guint), 1, fp);
  if (write && !FTT_CELL_IS_DESTROYED (root))
    (* write) (root, fp, data);

  if ((flags & FTT_FLAG_LEAF) == 0) {
    FttOct * oct;
    guint i;

    oct = root->children;
    for (i = 0; i < FTT_CELLS; i++)
      ftt_cell_write_binary (&(oct->cell[i]), max_depth, fp, write, data);
  }
}

#define FTT_CELL_IS_FLAGGED_LEAF(cell) (((cell)->flags & FTT_FLAG_LEAF) != 0)

static gboolean oct_read (FttCell * parent, 
			  GtsFile * fp,
			  FttCellReadFunc read,
			  gpointer data);

static gboolean cell_read (FttCell * cell, 
			   GtsFile * fp,
			   FttCellReadFunc read,
			   gpointer data)
{
  guint flags;

  if (fp->type != GTS_INT) {
    gts_file_error (fp, "expecting an integer (flags)");
    return FALSE;
  }
  flags = atoi (fp->token->str);
  if (FTT_CELL_ID (cell) != (flags & FTT_FLAG_ID)) {
    gts_file_error (fp, 
		    "FTT_CELL_ID (cell) `%d' != (flags & FTT_FLAG_ID) `%d'\n"
		    "Make sure the file has %d spatial dimensions",
		    FTT_CELL_ID (cell), (flags & FTT_FLAG_ID), FTT_DIMENSION);
    return FALSE;
  }
  cell->flags = flags;

  gts_file_next_token (fp);
  if (fp->type != '\n' && read && !FTT_CELL_IS_DESTROYED (cell))
    (* read) (cell, fp, data);
  if (fp->type == GTS_ERROR)
    return FALSE;
  gts_file_first_token_after (fp, '\n');

  if (!FTT_CELL_IS_DESTROYED (cell) && !FTT_CELL_IS_FLAGGED_LEAF (cell))
    return oct_read (cell, fp, read, data);

  cell->flags &= ~FTT_FLAG_LEAF;
  return TRUE;
}

static gboolean oct_read (FttCell * parent,
			  GtsFile * fp,
			  FttCellReadFunc read,
			  gpointer data)
{
  FttOct * oct;
  guint n;

  oct = g_malloc0 (sizeof (FttOct));
  oct->level = ftt_cell_level (parent);
  oct->parent = parent;
  parent->children = oct;
  ftt_cell_pos (parent, &(oct->pos));
  
  for (n = 0; n < FTT_CELLS; n++) {
    oct->cell[n].parent = oct;
    oct->cell[n].flags = n;
  }

  for (n = 0; n < FTT_CELLS; n++)
    if (!cell_read (&(oct->cell[n]), fp, read, data))
      return FALSE;
  
  return TRUE;
}

static void set_neighbors (FttCell * cell)
{
  ftt_cell_neighbors (cell, &(cell->children->neighbors));
}

/**
 * ftt_cell_read:
 * @fp: a #GtsFile.
 * @read: a #FttCellReadFunc function or %NULL.
 * @data: user data to pass to @read.
 *
 * If an error occurs (i.e. corrupted file or file format incorrect),
 * the @error field of @fp is set. A possibly incomplete tree is then
 * returned.
 *
 * Returns: the root cell of the tree contained in the file pointed to
 * by @fp. If not %NULL, the user-defined function @read is used to
 * read the extra user data associated with each cell.  
 */
FttCell * ftt_cell_read (GtsFile * fp,
			 FttCellReadFunc read,
			 gpointer data)
{
  FttCell * root;
  guint l, depth;

  g_return_val_if_fail (fp != NULL, NULL);

  root = ftt_cell_new (NULL, NULL);
  cell_read (root, fp, read, data);

  depth = ftt_cell_depth (root);
  for (l = 0; l < depth; l++)
    ftt_cell_traverse (root, FTT_PRE_ORDER, 
		       FTT_TRAVERSE_LEVEL|FTT_TRAVERSE_NON_LEAFS, l, 
		       (FttCellTraverseFunc) set_neighbors, NULL);

  return root;
}

static gboolean oct_read_binary (FttCell * parent, 
				 GtsFile * fp,
				 FttCellReadFunc read,
				 gpointer data);

static gboolean cell_read_binary (FttCell * cell, 
				  GtsFile * fp,
				  FttCellReadFunc read,
				  gpointer data)
{
  guint flags;

  if (gts_file_read (fp, &flags, sizeof (guint), 1) != 1) {
    gts_file_error (fp, "expecting an integer (flags)");
    return FALSE;
  }
  if (FTT_CELL_ID (cell) != (flags & FTT_FLAG_ID)) {
    gts_file_error (fp, 
		    "FTT_CELL_ID (cell) `%d' != (flags & FTT_FLAG_ID) `%d'\n"
		    "Make sure the file has %d spatial dimensions",
		    FTT_CELL_ID (cell), (flags & FTT_FLAG_ID), FTT_DIMENSION);
    return FALSE;
  }
  cell->flags = flags;

  if (read && !FTT_CELL_IS_DESTROYED (cell))
    (* read) (cell, fp, data);
  if (fp->type == GTS_ERROR)
    return FALSE;

  if (!FTT_CELL_IS_DESTROYED (cell) && !FTT_CELL_IS_FLAGGED_LEAF (cell))
    return oct_read_binary (cell, fp, read, data);

  cell->flags &= ~FTT_FLAG_LEAF;
  return TRUE;
}

static gboolean oct_read_binary (FttCell * parent,
				 GtsFile * fp,
				 FttCellReadFunc read,
				 gpointer data)
{
  FttOct * oct;
  guint n;

  oct = g_malloc0 (sizeof (FttOct));
  oct->level = ftt_cell_level (parent);
  oct->parent = parent;
  parent->children = oct;
  ftt_cell_pos (parent, &(oct->pos));
  
  for (n = 0; n < FTT_CELLS; n++) {
    oct->cell[n].parent = oct;
    oct->cell[n].flags = n;
  }

  for (n = 0; n < FTT_CELLS; n++)
    if (!cell_read_binary (&(oct->cell[n]), fp, read, data))
      return FALSE;
  
  return TRUE;
}

/**
 * ftt_cell_read_binary:
 * @fp: a #GtsFile.
 * @read: a #FttCellReadFunc function or %NULL.
 * @data: user data to pass to @read.
 *
 * If an error occurs (i.e. corrupted file or file format incorrect),
 * the @error field of @fp is set. A possibly incomplete tree is then
 * returned.
 *
 * Returns: the root cell of the tree contained in the file pointed to
 * by @fp. If not %NULL, the user-defined function @read is used to
 * read the extra user data associated with each cell.  
 */
FttCell * ftt_cell_read_binary (GtsFile * fp,
				FttCellReadFunc read,
				gpointer data)
{
  FttCell * root;
  guint l, depth;

  g_return_val_if_fail (fp != NULL, NULL);

  root = ftt_cell_new (NULL, NULL);
  cell_read_binary (root, fp, read, data);

  depth = ftt_cell_depth (root);
  for (l = 0; l < depth; l++)
    ftt_cell_traverse (root, FTT_PRE_ORDER, 
		       FTT_TRAVERSE_LEVEL|FTT_TRAVERSE_NON_LEAFS, l, 
		       (FttCellTraverseFunc) set_neighbors, NULL);

  return root;
}

/**
 * ftt_refine_corner:
 * @cell: a #FttCell.
 *
 * Returns: %TRUE if any "corner" neighbors of @cell are more than one
 * level more refined, %FALSE otherwise (see figure topology.fig).
 */
gboolean ftt_refine_corner (const FttCell * cell)
{
  FttCellNeighbors neighbor;
  guint i;

  g_return_val_if_fail (cell != NULL, FALSE);

  ftt_cell_neighbors (cell, &neighbor);
  for (i = 0; i < FTT_NEIGHBORS; i++) {
    FttCell * n = neighbor.c[i];

    if (n && !FTT_CELL_IS_LEAF (n)) {
      FttCellChildren child;
      guint j, k;

      k = ftt_cell_children_direction (n, FTT_OPPOSITE_DIRECTION (i), &child);
      for (j = 0; j < k; j++) {
	FttCell * c = child.c[j];

	if (c) {
#if FTT_2D
	  static guint perpendicular[FTT_NEIGHBORS_2D][FTT_CELLS/2] =
	  {{2,3},
	   {2,3},
	   {1,0},
	   {1,0}};
	  FttCell * nc = ftt_cell_neighbor (c, perpendicular[i][j]);

	  if (nc && !FTT_CELL_IS_LEAF (nc))
	    return TRUE;
#else  /* FTT_3D */
	  static guint perpendicular[FTT_NEIGHBORS][FTT_CELLS/2][2] =
	  {{{4,2},{4,3},{5,2},{5,3}},
	   {{4,2},{4,3},{5,2},{5,3}},
	   {{4,1},{4,0},{5,1},{5,0}},
	   {{4,1},{4,0},{5,1},{5,0}},
	   {{2,1},{2,0},{3,1},{3,0}},
	   {{2,1},{2,0},{3,1},{3,0}}};
	  FttCell * nc0, * nc1;

	  nc0 = ftt_cell_neighbor (c, perpendicular[i][j][0]);
	  if (nc0 && !FTT_CELL_IS_LEAF (nc0))
	    return TRUE;
	  nc1 = ftt_cell_neighbor (c, perpendicular[i][j][1]);
	  if (nc1 && !FTT_CELL_IS_LEAF (nc1))
	    return TRUE;
#endif /* FTT_3D */
	  if (!FTT_CELL_IS_LEAF (c)) {
	    FttCellChildren child;
	    guint j, k;

	    k = ftt_cell_children_direction (c, FTT_OPPOSITE_DIRECTION (i), &child);
	    for (j = 0; j < k; j++)
	      if (child.c[j])
		return TRUE;
	  }
	}
      }	
    }
  }

  return FALSE;
}

static void copy_cell (const FttCell * from,
		       FttCell * to,
		       FttCellCopyFunc copy,
		       gpointer data)
{
  to->flags = from->flags;

  if (!FTT_CELL_IS_DESTROYED (from)) {
    if (copy)
      (* copy) (from, to, data);
    
    if (!FTT_CELL_IS_LEAF (from)) {
      FttOct * oct_from = from->children;
      FttOct * oct_to;
      guint n;
      
      oct_new (to, FALSE, NULL, NULL);
      oct_to = to->children;
      for (n = 0; n < FTT_CELLS; n++)
	copy_cell (&(oct_from->cell[n]), &(oct_to->cell[n]), copy, data);
    }
  }
}

/**
 * ftt_cell_copy:
 * @root: the root of the cell tree to copy.
 * @copy: a #FttCellCopyFunc or %NULL.
 * @data: user data to pass to @copy.
 *
 * Returns: a new #FttCell root of the cell tree copy of @root. The
 * attributes of the cells are copied using the user-defined @copy
 * function.
 */
FttCell * ftt_cell_copy (const FttCell * root,
			 FttCellCopyFunc copy,
			 gpointer data)
{
  FttCell * root_copy;

  g_return_val_if_fail (root != NULL, NULL);

  root_copy = ftt_cell_new (NULL, NULL);
  ftt_cell_neighbors (root, &FTT_ROOT_CELL (root_copy)->neighbors);
  ftt_cell_pos (root, &FTT_ROOT_CELL (root_copy)->pos);
  FTT_ROOT_CELL (root_copy)->level = ftt_cell_level (root);
							   
  copy_cell (root, root_copy, copy, data);

  return root_copy;
}

#include "ftt_internal.c"

/**
 * ftt_face_traverse:
 * @root: the root #FttCell of the tree to traverse.
 * @c: only the faces orthogonal to this component will be traversed - one of
 * %FTT_X, %FTT_Y, (%FTT_Z), %FTT_XYZ.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children and faces are to be visited.
 * @max_depth: the maximum depth of the traversal. Cells below this
 * depth will not be traversed. If @max_depth is -1 all cells in the
 * tree are visited.
 * @func: the function to call for each visited #FttCellFace.
 * @data: user data to pass to @func.
 *
 * Traverses a cell tree starting at the given root #FttCell. Calls
 * the given function for each face of the cells of the tree.
 *
 * If %FTT_TRAVERSE_BOUNDARY_FACES is not set in @flags, only
 * "double-sided" faces are traversed i.e. the @neighbor field of the
 * face is never %NULL.  
 */
void ftt_face_traverse (FttCell * root,
			FttComponent c,
			FttTraverseType order,
			FttTraverseFlags flags,
			gint max_depth,
			FttFaceTraverseFunc func,
			gpointer data)
{
  FttDirection d;
  gpointer datum[6];
  gboolean check = FALSE;
  gboolean boundary_faces;

  g_return_if_fail (root != NULL);
  g_return_if_fail (c >= FTT_X && c <= FTT_XYZ);
  g_return_if_fail (func != NULL);

  boundary_faces = ((flags & FTT_TRAVERSE_BOUNDARY_FACES) != 0);
  datum[1] = &max_depth;
  datum[2] = func;
  datum[3] = data;
  datum[4] = &check;
  datum[5] = &boundary_faces;
  if (c == FTT_XYZ) {
    if (boundary_faces) {
      check = TRUE;
      ftt_cell_traverse (root, order, flags, max_depth, 
			 (FttCellTraverseFunc) traverse_all_faces, 
			 datum);
    }
    else {
      ftt_cell_traverse (root, order, flags, max_depth, 
			 (FttCellTraverseFunc) traverse_all_direct_faces, 
			 datum);
      check = TRUE;
      datum[0] = &d;
      for (d = 1; d < FTT_NEIGHBORS; d += 2)
	ftt_cell_traverse_boundary (root, d, order, flags, max_depth, 
				    (FttCellTraverseFunc) traverse_face, 
				    datum);
    }
  }
  else {
    if (boundary_faces) {
      check = TRUE;
      datum[0] = &c;
      ftt_cell_traverse (root, order, flags, max_depth, 
			 (FttCellTraverseFunc) traverse_face_component,
			 datum);
    }
    else {
      d = 2*c;
      datum[0] = &d;
      ftt_cell_traverse (root, order, flags, max_depth, 
			 (FttCellTraverseFunc) traverse_face_direction, datum);
      d = 2*c + 1;
      check = TRUE;
      ftt_cell_traverse_boundary (root, d, order, flags, max_depth, 
				  (FttCellTraverseFunc) traverse_face, datum);
    }
  }
  ftt_cell_traverse (root, order, flags, max_depth, 
		     (FttCellTraverseFunc) reset_flag, NULL);
}

static void traverse_face_boundary (FttCell * cell, gpointer * datum) 
{
  FttDirection * d = datum[0];
  FttFaceTraverseFunc func = (FttFaceTraverseFunc) datum[1];
  gpointer data = datum[2];
  FttCellFace face;
  
  face.d = *d;
  face.cell = cell;
  face.neighbor = ftt_cell_neighbor (cell, face.d);
  (* func) (&face, data);
}

/**
 * ftt_face_traverse_boundary:
 * @root: the root #FttCell of the tree to traverse.
 * @d: the direction of the boundary to visit.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @max_depth: the maximum depth of the traversal. Cells below this
 * depth will not be traversed. If @max_depth is -1 all cells in the
 * tree are visited.
 * @func: the function to call for each visited #FttCellFace.
 * @data: user data to pass to @func.
 *
 * Traverses a cell tree starting at the given root #FttCell. Calls
 * the given function for each face of the cell tree forming the
 * boundary of the domain in direction @d.  
 */
void ftt_face_traverse_boundary (FttCell * root,
				 FttDirection d,
				 FttTraverseType order,
				 FttTraverseFlags flags,
				 gint max_depth,
				 FttFaceTraverseFunc func,
				 gpointer data)
{
  gpointer datum[3];

  g_return_if_fail (root != NULL);
  g_return_if_fail (d < FTT_NEIGHBORS);
  g_return_if_fail (func != NULL);

  datum[0] = &d;
  datum[1] = func;
  datum[2] = data;
  ftt_cell_traverse_boundary (root, d, order, flags, max_depth, 
			      (FttCellTraverseFunc) traverse_face_boundary, 
			      datum);
}

/**
 * ftt_cell_coarsen:
 * @root: a #FttCell root of a cell tree to coarsen.
 * @coarsen: a #FttCellCoarsenFunc.
 * @coarsen_data: user data to pass to @coarsen.
 * @cleanup: a #FttCellCleanupFunc to call before destroying a cell or %NULL.
 * @cleanup_data: user data to pass to @cleanup.
 *
 * Coarsens the cell tree defined by @root according to @coarsen.
 *
 * Returns: %TRUE if @root has been coarsened (i.e. @root is now a
 * leaf cell), %FALSE otherwise.
 */
gboolean ftt_cell_coarsen (FttCell * root,
			   FttCellCoarsenFunc coarsen,
			   gpointer coarsen_data,
			   FttCellCleanupFunc cleanup,
			   gpointer cleanup_data)
{
  guint i, n;
  gboolean coarsenable = TRUE;
  
  g_return_val_if_fail (root != NULL, FALSE);
  g_return_val_if_fail (coarsen != NULL, FALSE);

  if (FTT_CELL_IS_LEAF (root))
    return (* coarsen) (root, coarsen_data);

  for (i = 0; i < FTT_CELLS; i++)
    if (!FTT_CELL_IS_DESTROYED (&(root->children->cell[i])))
      coarsenable &= ftt_cell_coarsen (&(root->children->cell[i]), 
				       coarsen, coarsen_data, 
				       cleanup, cleanup_data);
  if (!coarsenable || !(* coarsen) (root, coarsen_data))
    return FALSE;

  {
    FttDirection d;

    for (d = 0; d < FTT_NEIGHBORS; d++) {
      FttCellChildren child;

      n = ftt_cell_children_direction (root, d, &child);
      for (i = 0; i < n; i++) {
	FttCell * neighbor;

	if (child.c[i] && (neighbor = ftt_cell_neighbor (child.c[i], d)) &&
	    !FTT_CELL_IS_LEAF (neighbor)) {
	  FttCellChildren child1;
	  guint j, k;
	  gboolean empty = TRUE;

	  k = ftt_cell_children_direction (neighbor, FTT_OPPOSITE_DIRECTION (d), &child1);
	  for (j = 0; j < k && empty; j++)
	    if (child1.c[j])
	      empty = FALSE;
	  if (!empty && !ftt_cell_coarsen (neighbor, coarsen, coarsen_data, 
					   cleanup, cleanup_data))
	    return FALSE;
	  if (!FTT_CELL_IS_LEAF (neighbor))
	    neighbor->children->neighbors.c[FTT_OPPOSITE_DIRECTION (d)] = NULL;
	}
      }
    }
  }

  if (cleanup)
    for (i = 0; i < FTT_CELLS; i++)
      if (!FTT_CELL_IS_DESTROYED (&(root->children->cell[i])))
	(* cleanup) (&(root->children->cell[i]), cleanup_data);
  g_free (root->children);
  root->children = NULL;

  return TRUE;
}

/**
 * ftt_direction_from_name:
 * @name: a direction name.
 *
 * Returns: the index of the direction @name or %FTT_NEIGHBORS if
 * @name is not a valid direction name.  
 */
FttDirection ftt_direction_from_name (const gchar * name)
{
  FttDirection d = 0;

  g_return_val_if_fail (name != NULL, FTT_NEIGHBORS);

  while (d < FTT_NEIGHBORS && strcmp (name, ftt_direction_name[d]))
    d++;
  return d;
}

static void cell_traverse_add (FttCell * cell, GPtrArray * a)
{
  g_ptr_array_add (a, cell);
}

/**
 * ftt_cell_traverse_new:
 * @root: the root #FttCell of the tree to traverse.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @max_depth: the maximum depth of the traversal. Cells below this
 * depth will not be traversed. If @max_depth is -1 all cells in the
 * tree are visited.
 *
 * Returns: a new #FttCellTraverse.
 */
FttCellTraverse * ftt_cell_traverse_new (FttCell * root,
					 FttTraverseType order,
					 FttTraverseFlags flags,
					 gint max_depth)
{
  FttCellTraverse * t;
  GPtrArray * a;

  g_return_val_if_fail (root != NULL, NULL);

  a = g_ptr_array_new ();
  ftt_cell_traverse (root, order, flags, max_depth,
		     (FttCellTraverseFunc) cell_traverse_add, a);
  g_ptr_array_add (a, NULL);
  t = g_malloc (sizeof (FttCellTraverse));
  t->current = t->cells = (FttCell **) a->pdata;
  g_ptr_array_free (a, FALSE);
  return t;
}

/**
 * ftt_cell_traverse_rewind:
 * @t: a #FttCellTraverse.
 *
 * Sets @t at the begining of the traversal.
 */
void ftt_cell_traverse_rewind (FttCellTraverse * t)
{
  g_return_if_fail (t != NULL);

  t->current = t->cells;
}

/**
 * ftt_cell_traverse_destroy:
 * @t: a #FttCellTraverse.
 *
 * Frees all the memory associated with @t.
 */
void ftt_cell_traverse_destroy (FttCellTraverse * t)
{
  g_return_if_fail (t != NULL);

  g_free (t->cells);
  g_free (t);
}

