/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2008 National Institute of Water and Atmospheric
 * Research
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
 * \brief Conversion to unstructured mesh formats.
 */

#include "unstructured.h"
#include "variable.h"
#include "config.h"
#include "version.h"
#include "graphic.h"
#include "solid.h"

#define NV (4*(FTT_DIMENSION - 1))

static void reset_pointers (FttCell * cell, GfsVariable ** v)
{
  guint i;
  for (i = 0; i < NV; i++)
    GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, v[i])) = NULL;
}

typedef struct {
  FttCell * cell;
  guint i, index;
} Vertex;

/* Using VTK convention */
static FttDirection d[NV][FTT_DIMENSION] = {
#if FTT_2D
  {FTT_LEFT,FTT_BOTTOM}, {FTT_RIGHT,FTT_BOTTOM}, {FTT_LEFT,FTT_TOP}, {FTT_RIGHT,FTT_TOP},
#else /* 3D */
  {FTT_LEFT,FTT_BOTTOM,FTT_BACK}, {FTT_RIGHT,FTT_BOTTOM,FTT_BACK}, 
  {FTT_LEFT,FTT_TOP,FTT_BACK}, {FTT_RIGHT,FTT_TOP,FTT_BACK},
  {FTT_LEFT,FTT_BOTTOM,FTT_FRONT}, {FTT_RIGHT,FTT_BOTTOM,FTT_FRONT}, 
  {FTT_LEFT,FTT_TOP,FTT_FRONT}, {FTT_RIGHT,FTT_TOP,FTT_FRONT}
#endif /* 3D */
};

static void vertex_pos (Vertex * v, FttVector * p, GfsSimulation * sim)
{
  ftt_corner_pos (v->cell, d[v->i], p);
  gfs_simulation_map_inverse (sim, p);
}

static float vertex_value (Vertex * vertex, GfsVariable * v, gint max_depth)
{
  return gfs_dimensional_value (v, gfs_cell_corner_value (vertex->cell, d[vertex->i], 
							  v, max_depth));
}

typedef struct {
  GfsVariable ** v;
  GfsDomain * domain;
  GSList * vertices;
  gint max_depth;
  guint size, index;
} AllocParams;

static void allocate_vertices (FttCell * cell, AllocParams * par)
{
  static gint dx[NV][FTT_DIMENSION] = {
#if FTT_2D
    {-1,-1}, {1,-1}, {-1,1}, {1,1},
#else /* 3D */
    {-1,-1,-1}, {1,-1,-1}, {-1,1,-1}, {1,1,-1},
    {-1,-1,1},  {1,-1,1},  {-1,1,1},  {1,1,1}
#endif /* 3D */
  };

  gdouble h = ftt_cell_size (cell)/128.;
  guint i;
  for (i = 0; i < NV; i++)
    if (GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, par->v[i])) == NULL) {
      Vertex * vertex = g_malloc (par->size);
      vertex->i = i;
      vertex->cell = cell;
      vertex->index = par->index++;
      GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, par->v[i])) = vertex;
      par->vertices = g_slist_prepend (par->vertices, vertex);
      
      FttVector p;
      ftt_corner_pos (cell, d[i], &p);
      FttCell * neighbor[NV];
      guint j;
      for (j = 0; j < NV; j++)
	if (i != j) {
	  FttVector q;
	  FttComponent c;
	  for (c = 0; c < FTT_DIMENSION; c++)
	    (&q.x)[c] = (&p.x)[c] - dx[j][c]*h;
	  FttCell * n = gfs_domain_locate (par->domain, q, par->max_depth, NULL);
	  if (n) {
	    guint k;
	    for (k = 0; k < j && n; k++)
	      if (n == neighbor[k]) {
		/* T-junction */
#if DEBUG
		fprintf (stderr, "tj: %g %g %g %g %g %g\n", p.x, p.y, p.z, q.x, q.y, q.z);
#endif
		neighbor[k] = n = NULL;
	      }
	  }
	  neighbor[j] = n;
	}
	else
	  neighbor[j] = NULL;
      for (j = 0; j < NV; j++)
	if (neighbor[j]) {
	  g_assert (GFS_DOUBLE_TO_POINTER (GFS_VALUE (neighbor[j], par->v[j])) == NULL);
	  GFS_DOUBLE_TO_POINTER (GFS_VALUE (neighbor[j], par->v[j])) = vertex;
	}
    }
}

static GSList * allocate_domain_vertices (GfsDomain * domain, 
					  gint max_depth, 
					  GfsVariable * v[NV],
					  guint size)
{
  g_return_val_if_fail (size >= sizeof (Vertex), NULL);

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, max_depth,
			    (FttCellTraverseFunc) reset_pointers, v);
  AllocParams par;
  par.v = v;
  par.domain = domain;
  par.max_depth = max_depth;
  par.size = size;  
  par.vertices = NULL;
  par.index = 0;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, max_depth,
			    (FttCellTraverseFunc) allocate_vertices, &par);
  return g_slist_reverse (par.vertices);
}

#if DEBUG
static void print_pos (Vertex * v)
{
  FttVector p;
  ftt_corner_pos (v->cell, d[v->i], &p);
  fprintf (stderr, "v: %g %g %g\n", p.x, p.y, p.z);
}
#endif /* DEBUG */

#if DEBUG
static void draw_vertices (FttCell * cell, GfsVariable ** v)
{
  guint i;
  FttVector c;
  ftt_cell_pos (cell, &c);
  for (i = 0; i < NV; i++) {
    Vertex * vertex = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, v[i]));
    FttVector p;
    ftt_corner_pos (vertex->cell, d[vertex->i], &p);
    fprintf (stderr, "vp: %g %g\nvp: %g %g\nvp: \n", c.x, c.y, p.x, p.y);
  }    
}
#endif /* DEBUG */

typedef struct {
  FILE * fp;
  GfsVariable ** v;
} WriteParams;

static void write_element (FttCell * cell, WriteParams * par)
{
  fprintf (par->fp, "%d", NV);
  guint i;
  for (i = 0; i < NV; i++) {
    Vertex * v = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, par->v[i]));
    fprintf (par->fp, " %d", v->index);
  }
  fputc ('\n', par->fp);
}

static void cell_count (FttCell * cell, guint * count)
{
  (*count)++; 
}

static guint local_domain_size (GfsDomain * domain, gint max_depth)
{
  guint n = 0;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, max_depth,
			    (FttCellTraverseFunc) cell_count, &n);
  return n;
}

/**
 * gfs_domain_write_vtk:
 * @domain: a #GfsDomain.
 * @max_depth: the maximum depth to consider.
 * @variables: a list of #GfsVariable to output.
 * @precision: the formatting string for converting float to ASCII.
 * @fp: a file pointer.
 *
 * Writes in @fp a VTK-formatted representation of @domain and of the
 * corresponding variables in the given list.
 */
void gfs_domain_write_vtk (GfsDomain * domain, gint max_depth, GSList * variables, 
			   const gchar * precision, FILE * fp)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (precision != NULL);
  g_return_if_fail (fp != NULL);

  GfsVariable * v[NV];
  guint i;
  for (i = 0; i < NV; i++)
    v[i] = gfs_temporary_variable (domain);

  GSList * vertices = allocate_domain_vertices (domain, max_depth, v, sizeof (Vertex));

  /* header */
  fprintf (fp, 
	   "# vtk DataFile Version 2.0\n"
	   "Gerris simulation version %s (%s)\n"
	   "ASCII\n"
	   "DATASET UNSTRUCTURED_GRID\n"
	   "\n", 
	   GFS_VERSION,
	   GFS_BUILD_VERSION);
  
  /* vertices */
  guint nv = g_slist_length (vertices);
  fprintf (fp, "POINTS %d float\n", nv);
  gchar * format = g_strdup_printf ("%s %s %s\n", precision, precision, precision);
  GSList * j = vertices;
  while (j) {
    FttVector p;
    vertex_pos (j->data, &p, GFS_SIMULATION (domain));
    fprintf (fp, format, p.x, p.y, p.z);
    j = j->next;
  }
  g_free (format);
  fputc ('\n', fp);

  /* elements */
  guint n_cells = local_domain_size (domain, max_depth);
  fprintf (fp, "CELLS %d %d\n", n_cells, n_cells*(NV + 1));
  WriteParams par;
  par.v = v;
  par.fp = fp;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, max_depth,
			    (FttCellTraverseFunc) write_element, &par);
  fprintf (fp, "\nCELL_TYPES %d\n",n_cells);
  for (i = 0; i < n_cells; i++) {
#if FTT_2D
    fputs ("8\n", fp);
#else
    fputs ("11\n", fp);
#endif
  }
  fputc ('\n', fp);

#if DEBUG
  fprintf (stderr, "vertices: %d\n", g_slist_length (vertices));
  g_slist_foreach (vertices, (GFunc) print_pos, NULL);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, max_depth,
			    (FttCellTraverseFunc) draw_vertices, v);
#endif /* DEBUG */
  
  /* write scalar fields */
  if (variables) {
    gchar * format = g_strdup_printf ("%s\n", precision);
    fprintf (fp, "POINT_DATA %d\n", nv);
    GSList * i = variables;
    while (i) {
      GfsVariable * v = i->data;
      fprintf (fp, "SCALARS %s float\nLOOKUP_TABLE default\n", v->name);
      GSList * j = vertices;
      while (j) {
	Vertex * vertex = j->data;
	fprintf (fp, format, vertex_value (vertex, v, max_depth));
	j = j->next;
      }
      fputc ('\n', fp);
      i = i->next;
    }
    g_free (format);
  }

  /* cleanup */
  g_slist_foreach (vertices, (GFunc) g_free, NULL);
  g_slist_free (vertices);
  for (i = 0; i < NV; i++)
    gts_object_destroy (GTS_OBJECT (v[i]));
}

static void write_tecplot_element (FttCell * cell, WriteParams * par)
{
  static guint tecplot_index[NV] = {
#if FTT_2D
    0, 1, 3, 2
#else /* 3D */
    0, 1, 3, 2,
    4, 5, 7, 6
#endif /* 3D */
  };
  guint i;
  for (i = 0; i < NV; i++) {
    Vertex * v = GFS_DOUBLE_TO_POINTER (GFS_VALUE (cell, par->v[tecplot_index[i]]));
    fprintf (par->fp, "%d ", v->index + 1);
  }
  fputc ('\n', par->fp);
}

/**
 * gfs_domain_write_tecplot:
 * @domain: a #GfsDomain.
 * @max_depth: the maximum depth to consider.
 * @variables: a list of #GfsVariable to output.
 * @precision: the formatting string for converting float to ASCII.
 * @fp: a file pointer.
 *
 * Writes in @fp a Tecplot-formatted representation of @domain and of the
 * corresponding variables in the given list.
 */
void gfs_domain_write_tecplot (GfsDomain * domain, gint max_depth, GSList * variables, 
			       const gchar * precision, FILE * fp)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (precision != NULL);
  g_return_if_fail (fp != NULL);

  GfsVariable * v[NV];
  guint i;
  for (i = 0; i < NV; i++)
    v[i] = gfs_temporary_variable (domain);

  GSList * vertices = allocate_domain_vertices (domain, max_depth, v, sizeof (Vertex));

  /* header */
  fprintf (fp,
	   " TITLE = \"Gerris simulation version %s (%s)\"\n",
	   GFS_VERSION,
	   GFS_BUILD_VERSION);

  fputs (FTT_DIMENSION == 2 ? " VARIABLES = \"X\", \"Y\"" : " VARIABLES = \"X\", \"Y\", \"Z\"", fp);
  GSList * j = variables;
  while (j) {
    GfsVariable * v = j->data;
    fprintf (fp, ", \"%s\"", v->name);
    j = j->next;
  }
  fputc ('\n', fp);

  guint nv = g_slist_length (vertices);
  guint n_cells = local_domain_size (domain, max_depth);
  fprintf (fp, " ZONE N=%i, E=%i, F=FEPOINT, ", nv, n_cells);
  fputs (FTT_DIMENSION == 2 ? "ET=QUADRILATERAL\n" : "ET=BRICK\n", fp);
  
  /* vertices and scalar data */
  gchar * xyzformat = 
#if FTT_2D
    g_strdup_printf ("%s %s", precision, precision);
#else
    g_strdup_printf ("%s %s %s", precision, precision, precision);
#endif
  gchar * format = g_strdup_printf (" %s", precision);
  j = vertices;
  while (j) {
    Vertex * vertex = j->data;
    FttVector p;
    vertex_pos (vertex, &p, GFS_SIMULATION (domain));
    fprintf (fp, xyzformat, p.x, p.y, p.z);
    GSList * k = variables;
    while (k) {
      fprintf (fp, format, vertex_value (vertex, k->data, max_depth));
      k = k->next;
    }
    fputc ('\n', fp);
    j = j->next;
  }
  g_free (format);
  g_free (xyzformat);

  /* elements */
  WriteParams par;
  par.v = v;
  par.fp = fp;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, max_depth,
			    (FttCellTraverseFunc) write_tecplot_element, &par);

  /* cleanup */
  g_slist_foreach (vertices, (GFunc) g_free, NULL);
  g_slist_free (vertices);
  for (i = 0; i < NV; i++)
    gts_object_destroy (GTS_OBJECT (v[i]));
}

#if !FTT_2D /* gfs_domain_write_tecplot_surface */

typedef struct {
  FttCell * cell;             /* cell */
  guint nnodes;               /* number of nodes (3 or 4) */
  FttVector node[4];          /* node coordinates */
} WallFace;

typedef struct {
  GfsDomain * domain;
  GSList * wallfaces;
  guint nnodes_total;
} AllocateFacesParams;

static void allocate_wallfaces (FttCell * cell, AllocateFacesParams * par)
{
  if (!GFS_IS_MIXED (cell)) return;

  GfsSolidVector * s = GFS_STATE (cell)->solid;
  FttVector normal; /* solid normal */
  gfs_solid_normal (cell, &normal);
  gts_vector_normalize (&normal.x);

  FttDirection d[12];
  /* array of node coordinates for a cut face */
  FttVector nodecutface[FTT_DIMENSION*(FTT_DIMENSION - 1) + 1];
  guint inode, inode2, jnode_max_sintheta = 0, nnodes = 
    gfs_cut_cube_vertices (cell, -1, &s->ca, &normal, nodecutface, d, NULL, NULL);
  g_assert (nnodes <= 6);
  
  guint nnodecutface = nnodes;

  if (nnodecutface > 3) {   /* reorder faces if necessary */
    /* Tecplot can think that opposite vertices of the quadrilateral surface 
       element are connected. This may result in ugly X-shaped surface elements. 
       Reorder the array of node if necessary, using bubble sort, 
       to maximize the sine of the angle between the vector from each 
       intersection point to the cell center and the vector from this point 
       to the next point in the array */

    FttVector facecenter = {s->ca.x, s->ca.y, s->ca.z};

    guint i_switchnodes = 0;        /* counter to avoid infinite loop */    
    gboolean switchnodes = FALSE;   /* logical variable used to reorder cut face nodes */

    do {
      i_switchnodes++;

      for (inode = 0; inode < nnodecutface; inode++) {
	FttVector node = nodecutface[inode];                /* face node coordinates */
	FttVector diff1 = {facecenter.x - node.x, 
			   facecenter.y - node.y, 
			   facecenter.z - node.z};
	gdouble length_diff1 = ftt_vector_norm (&diff1);

	gdouble max_sintheta = 0.;
	/* cycle through all other nodes (jnode) where cut face intersects cell edges */
	for (inode2 = 1; inode2 < nnodecutface; inode2++) {
	  gint jnode = (inode + inode2)%nnodecutface;
	  FttVector diff2 = {nodecutface[jnode].x - node.x, 
			     nodecutface[jnode].y - node.y,
			     nodecutface[jnode].z - node.z};
	  gdouble length_diff2 = ftt_vector_norm (&diff2);

	  gdouble sintheta = ((diff1.y*diff2.z - diff1.z*diff2.y)*normal.x + 
			      (diff1.z*diff2.x - diff1.x*diff2.z)*normal.y + 
			      (diff1.x*diff2.y - diff1.y*diff2.x)*normal.z)/
	    (length_diff1*length_diff2);

	  if (sintheta > max_sintheta) {
	    max_sintheta = sintheta;
	    jnode_max_sintheta = jnode;
	  }	  
	}
	
	/* terminate if cannot find positive angle between cut face nodes */
	g_assert (max_sintheta != 0.);
	
	inode2 = (inode + 1)%nnodecutface;
	if (jnode_max_sintheta != inode2) {
	  node = nodecutface[jnode_max_sintheta];
	  nodecutface[jnode_max_sintheta] = nodecutface[inode2];
	  nodecutface[inode2] = node;
	  
	  switchnodes = TRUE;
	}
      } /* inode-loop */
    } while (switchnodes && i_switchnodes < 1000);   /* avoid infinite loop */    
  } /* reorder faces if necessary */

  /* If there are more than 4 points in the array, 
     divide the cut face into a quadrilateral face and 
     either a quadrilateral or triangular face */

  /* assign data to nodeinfo array, increment number of wall faces and number of nodes */
  if (nnodecutface <= 4) {
    WallFace * face = g_malloc (sizeof (WallFace));

    for (inode = 0; inode < nnodecutface; inode++)
      face->node[inode] = nodecutface[inode];
    face->cell = cell;
    face->nnodes = nnodecutface;
    
    par->nnodes_total += nnodecutface;
    par->wallfaces = g_slist_append (par->wallfaces, face);
  }
  else { /* cut face must be divided into 2 quadrilateral/triangular faces */
    /* first face is quadrilateral */
    WallFace * face1 = g_malloc (sizeof (WallFace));
    
    for (inode = 0; inode < 4; inode++)
      face1->node[inode] = nodecutface[inode];
    face1->cell = cell;
    face1->nnodes = 4;

    par->wallfaces = g_slist_append (par->wallfaces, face1);
    par->nnodes_total += 4;

    /* second face is triangular if nnodecutface=5; otherwise quadrilateral */
    WallFace * face2 = g_malloc (sizeof (WallFace));

    nnodecutface -= 2;
    for (inode = 0; inode < nnodecutface; inode++)
      face2->node[inode] = 
	nodecutface[(inode + 3)%(nnodecutface + 2)];
    face2->cell = cell;
    face2->nnodes = nnodecutface;

    par->wallfaces = g_slist_append (par->wallfaces, face2);
    par->nnodes_total += nnodecutface;
  } /* cut face must be divided into 2 quadrilateral/triangular faces */
}

/**
 * gfs_domain_write_tecplot_surface:
 * @domain: a #GfsDomain.
 * @max_depth: the maximum depth to consider.
 * @variables: a list of #GfsVariable to output.
 * @precision: the formatting string for converting float to ASCII.
 * @fp: a file pointer.
 *
 * Writes in @fp a Tecplot-formatted representation of the solid
 * surface of @domain and of the corresponding variables in the given
 * list.
 */
void gfs_domain_write_tecplot_surface (GfsDomain * domain, gint max_depth, GSList * variables, 
				       const gchar * precision, FILE * fp)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (precision != NULL);
  g_return_if_fail (fp != NULL);

  GSList * j;
  AllocateFacesParams par = {domain, NULL, 0};

  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, max_depth,
			    (FttCellTraverseFunc) allocate_wallfaces, &par);

#if 0 /* include header for standalone event */
  fprintf (fp,
	   " TITLE = \"Gerris simulation version %s (%s)\"\n",
	   GFS_VERSION,
	   GFS_BUILD_VERSION);

  fputs (FTT_DIMENSION == 2 ? " VARIABLES = 'X', 'Y'" : " VARIABLES = 'X', 'Y', 'Z'", fp);
  j = variables;
  while (j) {
    GfsVariable * v = j->data;
    fprintf (fp, ", '%s'", v->name);
    j = j->next;
  }
  fputc ('\n', fp);
#endif /* include header for standalone event */

  guint inode, nwallfaces = g_slist_length (par.wallfaces);
  fprintf (fp, " ZONE T=\"WALL DATA\", N=%i, E=%i, F=FEPOINT, ", par.nnodes_total, nwallfaces);
  fputs ("ET=QUADRILATERAL\n", fp);

  /* nodes and scalar data */
  gchar * xyzformat = 
#if FTT_2D
    g_strdup_printf ("%s %s", precision, precision);
#else
    g_strdup_printf ("%s %s %s", precision, precision, precision);
#endif
  gchar * format = g_strdup_printf (" %s", precision);
 
  j = par.wallfaces;
  while (j) {
    WallFace * face = j->data;
    FttCell * cell = face->cell;

    for (inode = 0; inode < face->nnodes; inode++) {
      FttVector * node = &face->node[inode];
      fprintf (fp, xyzformat, node->x, node->y, node->z);

      GSList * k = variables;
      while (k) {
	GfsVariable * v = k->data;
	/* write out cell values; can use interpolate for more precise data representation */
	fprintf (fp, format, GFS_VALUE (cell, v));
	k = k->next;
      }
      fputc ('\n', fp);
    }

    j = j->next;
  }
  g_free (format);
  g_free (xyzformat);

  /* output node connectivity information */
  inode = 1; 
  j = par.wallfaces;
  while (j) {
    WallFace * face = j->data;

    if (face->nnodes == 4)
      fprintf (fp, "%d %d %d %d\n", inode, inode + 1, inode + 2, inode + 3);
    else                   /* triangular face */
      fprintf (fp, "%d %d %d %d\n", inode, inode + 1, inode + 2, inode + 2);

    inode += face->nnodes;
    j = j->next;
  }

  /* cleanup */
  g_slist_foreach (par.wallfaces, (GFunc) g_free, NULL);
  g_slist_free (par.wallfaces);
}

#endif /* gfs_domain_write_tecplot_surface */
