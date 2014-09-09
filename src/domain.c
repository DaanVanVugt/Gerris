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
 * \brief Spatial domain.
 */

#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "domain.h"

#include "advection.h"
#include "source.h"
#include "solid.h"
#include "adaptive.h"
#include "mpi_boundary.h"
#include "metric.h"
#include "version.h"
#include "init.h"

#include "config.h"

/* GfsLocateArray: Object */

static void locate_index (FttVector * p, GfsLocateArray * a, gint i[FTT_DIMENSION])
{
  gint c;
  for (c = 0; c < FTT_DIMENSION; c++)
    i[c] = floor (((&p->x)[c] - a->min[c])/a->h);
}

static void root_bounds (FttCell * root, GfsLocateArray * a)
{
  FttVector p;
  ftt_cell_pos (root, &p);
  gint i;
  for (i = 0; i < FTT_DIMENSION; i++) {
    if ((&p.x)[i] + a->h/2. > a->max[i]) a->max[i] = (&p.x)[i] + a->h/2.;
    if ((&p.x)[i] - a->h/2. < a->min[i]) a->min[i] = (&p.x)[i] - a->h/2.;
  }
}

static void box_bounds (GfsBox * box, GfsLocateArray * a)
{
  root_bounds (box->root, a);
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d]))
      root_bounds (GFS_BOUNDARY (box->neighbor[d])->root, a);
}

static gint locate_linear_index (FttVector * p, GfsLocateArray * a)
{
  gint i[FTT_DIMENSION], index = 0, c;
  locate_index (p, a, i);
  for (c = 0; c < FTT_DIMENSION; c++) {
    if (i[c] < 0 || i[c] >= a->n[c])
      return -1;
    index = index*a->n[c] + i[c];
  }
  return index;
}

static void box_index (GfsBox * b, GfsLocateArray * a)
{
  FttVector p;
  ftt_cell_pos (b->root, &p);
  gint i = locate_linear_index (&p, a);
  g_assert (i >= 0);
  g_assert (!a->root[i]);
  a->root[i] = g_slist_prepend (NULL, b);
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (b->neighbor[d])) {
      GfsBoundary * boundary = GFS_BOUNDARY (b->neighbor[d]);
      ftt_cell_pos (boundary->root, &p);
      gint i = locate_linear_index (&p, a);
      g_assert (i >= 0);
      a->root[i] = g_slist_prepend (a->root[i], boundary);
    }
}

/**
 * @domain: a #GfsDomain.
 *
 * Creates a rectangular array for fast location of which GfsBox
 * contains a given point.
 *
 * Returns: a new #GfsLocateArray.
 */
GfsLocateArray * gfs_locate_array_new (GfsDomain * domain)
{
  g_return_val_if_fail (domain != NULL, NULL);
  
  GfsLocateArray * a = g_malloc (sizeof (GfsLocateArray));
  guint i;
  a->h = ftt_level_size (domain->rootlevel);
  for (i = 0; i < FTT_DIMENSION; i++) {
    a->min[i] = G_MAXDOUBLE;
    a->max[i] = - G_MAXDOUBLE;
  }
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_bounds, a);
  guint size = 1;
  for (i = 0; i < FTT_DIMENSION; i++) {
    g_assert (a->max[i] > a->min[i]);
    a->n[i] = ceil ((a->max[i] - a->min[i])/a->h - 0.5);
    size *= a->n[i];
  }
  a->root = g_malloc0 (size*sizeof (GSList *));
  a->size = size;
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_index, a);
  return a;
}

/**
 * @a: a #GfsLocateArray.
 * @p: a #FttVector.
 *
 * Returns: a list of objects containing @p or %NULL.
 */
GSList * gfs_locate_array_locate (GfsLocateArray * a, FttVector * p)
{
  g_return_val_if_fail (a != NULL, NULL);
  g_return_val_if_fail (p != NULL, NULL);
  gint i = locate_linear_index (p, a);
  return i < 0 ? NULL : a->root[i];
}

/**
 * @a: a #GfsLocateArray.
 *
 * Frees the memory allocated to @a.
 */
void gfs_locate_array_destroy (GfsLocateArray * a)
{
  if (a) {
    gint i;
    for (i = 0; i < a->size; i++)
      g_slist_free (a->root[i]);
    g_free (a->root);
    g_free (a);
  }
}

/**
 * Spatial domain.
 * \beginobject{GfsDomain}
 */

static void domain_write (GtsObject * o, FILE * fp)
{
  GfsDomain * domain = GFS_DOMAIN (o);

  if (GTS_OBJECT_CLASS (gfs_domain_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_domain_class ())->parent_class->write) (o, fp);

  fputs (" { ", fp);
  if (domain->rootlevel != 0)
    fprintf (fp, "rootlevel = %u ", domain->rootlevel);
  if (domain->refpos.x != 0.)
    fprintf (fp, "x = %g ", domain->refpos.x);
  if (domain->refpos.y != 0.)
    fprintf (fp, "y = %g ", domain->refpos.y);
  if (domain->refpos.z != 0.)
    fprintf (fp, "z = %g ", domain->refpos.z);
  if (domain->lambda.x != 1.)
    fprintf (fp, "lx = %g ", domain->lambda.x);
  if (domain->lambda.y != 1.)
    fprintf (fp, "ly = %g ", domain->lambda.y);
  if (domain->lambda.z != 1.)
    fprintf (fp, "lz = %g ", domain->lambda.z);
  fprintf (fp, "version = %d ", atoi (GFS_BUILD_VERSION));
  if (!domain->overlap)
    fputs ("overlap = 0 ", fp);
  if (domain->max_depth_write > -2) {
    GSList * i = domain->variables_io;

    if (i != NULL) {
      fprintf (fp, "variables = %s", GFS_VARIABLE (i->data)->name);
      i = i->next;
      while (i) {
	fprintf (fp, ",%s", GFS_VARIABLE (i->data)->name);
	i = i->next;
      }
      fputc (' ', fp);
    }
  }
  if (domain->binary != FALSE)
    fprintf (fp, "binary = 1 ");
  fputc ('}', fp);
}

static void domain_read (GtsObject ** o, GtsFile * fp)
{
  GfsDomain * domain = GFS_DOMAIN (*o);
  GtsFileVariable var[] = {
    {GTS_UINT,   "rootlevel", TRUE},
    {GTS_DOUBLE, "x",         TRUE},
    {GTS_DOUBLE, "y",         TRUE},
    {GTS_DOUBLE, "z",         TRUE},
    {GTS_DOUBLE, "lx",        TRUE},
    {GTS_DOUBLE, "ly",        TRUE},
    {GTS_DOUBLE, "lz",        TRUE},
    {GTS_STRING, "variables", TRUE},
    {GTS_INT,    "binary",    TRUE},
    {GTS_INT,    "version",   TRUE},
    {GTS_INT,    "overlap",   TRUE},
    {GTS_NONE}
  };
  gchar * variables = NULL;

  if (GTS_OBJECT_CLASS (gfs_domain_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_domain_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  domain->version = -1;
  var[0].data = &domain->rootlevel;
  var[1].data = &domain->refpos.x;
  var[2].data = &domain->refpos.y;
  var[3].data = &domain->refpos.z;
  var[4].data = &domain->lambda.x;
  var[5].data = &domain->lambda.y;
  var[6].data = &domain->lambda.z;
  var[7].data = &variables;
  var[8].data = &domain->binary;
  var[9].data = &domain->version;
  var[10].data = &domain->overlap;
  gts_file_assign_variables (fp, var);
  if (fp->type == GTS_ERROR) {
    g_free (variables);
    return;
  }

  if (var[4].set || var[5].set || var[6].set)
    g_warning ("the (lx,ly,lz) parameters are obsolete, please use GfsMetricStretch instead");

#if FTT_2D
  if (var[3].set) {
    gts_file_variable_error (fp, var, "z", "unknown identifier `z'");
    return;
  }
  if (var[6].set) {
    gts_file_variable_error (fp, var, "lz", "unknown identifier `lz'");
    return;
  }
#endif

  if (var[4].set && domain->lambda.x <= 0.) {
    gts_file_variable_error (fp, var, "lx", "lx must be strictly positive");
    return;
  }
  if (var[5].set && domain->lambda.y <= 0.) {
    gts_file_variable_error (fp, var, "ly", "ly must be strictly positive");
    return;
  }
  if (var[6].set && domain->lambda.z <= 0.) {
    gts_file_variable_error (fp, var, "lz", "lz must be strictly positive");
    return;
  }

  if (variables != NULL) {
    gchar * variables1, * s;

    variables1 = g_strdup (variables);
    s = strtok (variables1, ",");
    while (s) {
      gfs_domain_add_variable (domain, s, NULL);
      s = strtok (NULL, ",");
    }
    g_free (variables1);
    domain->variables_io = gfs_variables_from_list (domain->variables, variables, &s);
    g_free (variables);
  } 
}

static void box_set_pos (GfsBox * box, FttVector * pos, 
			 FttDirection dold)
{
  FttVector p;
  static FttVector rpos[FTT_NEIGHBORS] = {
#if FTT_2D
    {1.,0.,0.}, {-1.,0.,0.}, {0.,1.,0.}, {0.,-1.,0.}
#else  /* FTT_3D */
    {1.,0.,0.}, {-1.,0.,0.}, {0.,1.,0.}, {0.,-1.,0.}, {0.,0.,1.}, {0.,0.,-1.}
#endif /* FTT_3D */
  };  
  static FttDirection id[FTT_NEIGHBORS][FTT_NEIGHBORS] = {
#if FTT_2D
    {0,1,2,3},
    {1,0,3,2},
    {2,3,1,0},
    {3,2,0,1},
#else  /* 3D */
    {0,1,2,3,5,4},
    {1,0,3,2,4,5},
    {2,3,1,0,5,4},
    {3,2,0,1,4,5},
    {4,5,2,3,0,1},
    {5,4,3,2,1,0}
#endif /* 3D */
  };

  ftt_cell_pos (box->root, &p);
  if (p.x != G_MAXDOUBLE) /* position already set */
    return;

  FttDirection i;
  gdouble size;
  size = ftt_cell_size (box->root);
  ftt_cell_set_pos (box->root, pos);
  for (i = 0; i < FTT_NEIGHBORS; i++) {
    FttDirection d = id[dold][i];
    
    p.x = pos->x + rpos[d].x*size;
    p.y = pos->y + rpos[d].y*size;
    p.z = pos->z + rpos[d].z*size;
    if (GFS_IS_BOX (box->neighbor[d]))
      box_set_pos (GFS_BOX (box->neighbor[d]), &p, d);
    else if (GFS_IS_BOUNDARY (box->neighbor[d]))
      ftt_cell_set_pos (GFS_BOUNDARY (box->neighbor[d])->root, &p);
  }
}

static void set_ref_pos (GfsBox * box, FttVector * pos)
{
  if (box->id == 1)
    box_set_pos (box, pos, FTT_RIGHT);
}

static void pid_max (GfsBox * box, gint * np)
{
  if (box->pid > *np)
    *np = box->pid;  
}

typedef struct {
  GSList * removed;
  gint pid;
} RemovedData;

static void removed_list (GfsBox * box, RemovedData * p)
{
  if (box->pid != p->pid)
    p->removed = g_slist_prepend (p->removed, box);
  else {
    FttDirection d;
    GfsBox * matching;

    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (GFS_IS_BOUNDARY_PERIODIC (box->neighbor[d]) &&
	  !GFS_IS_BOUNDARY_MPI (box->neighbor[d]) &&
	  (matching = GFS_BOUNDARY_PERIODIC (box->neighbor[d])->matching)->pid != p->pid) {
	GfsBoundaryPeriodic * b = GFS_BOUNDARY_PERIODIC (box->neighbor[d]);
	FttDirection rotate = b->d;
	gdouble orientation = b->rotate;
	gts_object_destroy (GTS_OBJECT (b));
	b = GFS_BOUNDARY_PERIODIC (gfs_boundary_mpi_new (gfs_boundary_mpi_class (), 
							 box, d, matching->pid, matching->id));
	if (orientation != 0.)
	  gfs_boundary_periodic_rotate (b, rotate, orientation);
      }
  }
}

static void mpi_links (GfsBox * box, GfsDomain * domain)
{
  FttDirection d;
  GtsObject * neighbor[FTT_NEIGHBORS];
  gint pid = box->pid;
  gint id = box->id;

  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOX (box->neighbor[d]) && GFS_BOX (box->neighbor[d])->pid == domain->pid)
      neighbor[d] = box->neighbor[d];
    else
      neighbor[d] = NULL;
  gts_object_destroy (GTS_OBJECT (box));

  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (neighbor[d])
      gfs_boundary_mpi_new (gfs_boundary_mpi_class (),
			    GFS_BOX (neighbor[d]), 
			    FTT_OPPOSITE_DIRECTION (d), 
			    pid, id);
}

static void add_id (GfsBox * box, GPtrArray * ids)
{
  if (box->id > ids->len)
    g_ptr_array_set_size (ids, box->id);
  g_ptr_array_index (ids, box->id - 1) = box;
}

static GPtrArray * box_ids (GfsDomain * domain)
{
  GPtrArray * ids = g_ptr_array_new ();
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) add_id, ids);
  return ids;
}

static void convert_boundary_mpi_into_edges (GfsBox * box, GPtrArray * ids)
{
  gint pid = gfs_box_domain (box)->pid;
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY_MPI (box->neighbor[d])) {
      GfsBoundaryMpi * b = GFS_BOUNDARY_MPI (box->neighbor[d]);
      GfsBox * nbox;
      if (b->id >= 0 && b->id <= ids->len && (nbox = g_ptr_array_index (ids, b->id - 1))) {
	FttDirection nd = GFS_BOUNDARY_PERIODIC (b)->d;
	if (!GFS_IS_BOUNDARY_MPI (nbox->neighbor[nd]))
	  g_warning ("!GFS_IS_BOUNDARY_MPI (nbox->neighbor[nd])");
	else {
	  GfsBoundaryMpi * nb = GFS_BOUNDARY_MPI (nbox->neighbor[nd]);
	  if (box->id != nb->id)
	    g_warning ("box->id != nb->id");
	  else {
	    int rotate = GFS_BOUNDARY_PERIODIC (b)->rotate;
	    gts_object_destroy (GTS_OBJECT (b));
	    gts_object_destroy (GTS_OBJECT (nb));
	    GfsGEdge * edge;
	    if (nd == FTT_OPPOSITE_DIRECTION (d))
	      /* standard edge */
	      edge = gfs_gedge_new (gfs_gedge_class (), box, nbox, d);
	    else {
	      /* rotated edge */
	      g_assert (rotate);
	      if (rotate > 0) {
		edge = gfs_gedge_new (gfs_gedge_class (), box, nbox, d);
		edge->rotate = nd;		
	      }
	      else {
		edge = gfs_gedge_new (gfs_gedge_class (), nbox, box, nd);
		edge->rotate = d;
	      }
	    }
	    gfs_gedge_link_boxes (edge);
	  }
	}
      }
    }
  if (pid >= 0)
    box->pid = pid;
}

static void domain_post_read (GfsDomain * domain, GtsFile * fp)
{
  gts_graph_foreach_edge (GTS_GRAPH (domain), (GtsFunc) gfs_gedge_link_boxes, NULL);

  domain->np = 0;
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) pid_max, &domain->np);
#ifdef HAVE_MPI
  if (domain->pid >= 0) { 
    /* Multiple PEs, make sure we have the max pid over all the boxes,
       in case each process loads a different file */
    int npmax;
    MPI_Allreduce (&domain->np, &npmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    domain->np = npmax;
  }
#endif /* HAVE_MPI */
  domain->np++; /* number of PEs according to pids */

  if (domain->np > 1 && domain->pid >= 0) { /* Multiple PEs */
    RemovedData p = { NULL, domain->pid };
    
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) set_ref_pos, &domain->refpos);
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) removed_list, &p);
#ifdef HAVE_MPI
    int comm_size;
    MPI_Comm_size (MPI_COMM_WORLD, &comm_size);
    if (domain->np != comm_size) {
      g_slist_free (p.removed);
      gts_file_error (fp, "it would be valid if one or %d PE were used", domain->np);
      return;
    }
#endif /* HAVE_MPI */
    g_slist_foreach (p.removed, (GFunc) mpi_links, domain);
    g_slist_free (p.removed);
  }
  else { /* Single PE */
    /* Create array for fast linking of ids to GfsBox pointers */
    GPtrArray * ids = box_ids (domain);
    
    /* Convert GfsBoundaryMpi into graph edges */
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) convert_boundary_mpi_into_edges, ids);

    g_ptr_array_free (ids, TRUE);

    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) set_ref_pos, &domain->refpos);
  }

  gfs_domain_match (domain);

  gfs_locate_array_destroy (domain->array);
  domain->array = gfs_locate_array_new (domain);

  domain->version = atoi (GFS_BUILD_VERSION);
}

static void free_pair (gpointer key, gpointer value)
{
  g_free (key);
  g_free (value);
}

static void cleanup_each_box (GfsBox * box, GfsDomain * domain)
{
  /* this is a necessary check when using graph partitioning */
  if (g_slist_length (GTS_SLIST_CONTAINEE (box)->containers) == 1) {
    ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
		       (FttCellTraverseFunc) gfs_cell_cleanup, domain);
    FttDirection d;
    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (GFS_IS_BOUNDARY (box->neighbor[d]))
	ftt_cell_traverse (GFS_BOUNDARY (box->neighbor[d])->root, 
			   FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			   (FttCellTraverseFunc) gfs_cell_cleanup, domain);
  }
}

static void domain_destroy (GtsObject * o)
{
  GfsDomain * domain = GFS_DOMAIN (o);
  GSList * i;

  gfs_clock_destroy (domain->timer);
  g_timer_destroy (domain->clock);

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) cleanup_each_box, domain);

  i = domain->variables;
  while (i) {
    GSList * next = i->next;
    gts_object_destroy (i->data);
    i = next;
  }
  g_assert (domain->variables == NULL);

  g_slist_foreach (domain->derived_variables, (GFunc) gts_object_destroy, NULL);
  g_slist_free (domain->derived_variables);
  domain->derived_variables = NULL;

  g_array_free (domain->allocated, TRUE);

  g_hash_table_foreach (domain->timers, (GHFunc) free_pair, NULL);
  g_hash_table_destroy (domain->timers);

  g_slist_free (domain->variables_io);

  gfs_locate_array_destroy (domain->array);
  domain->array = NULL;

  g_hash_table_destroy (domain->objects);

  g_ptr_array_free (domain->sorted, TRUE);
  domain->sorted = NULL;

  (* GTS_OBJECT_CLASS (gfs_domain_class ())->parent_class->destroy) (o);
}

static void add_item (gpointer item, GPtrArray * a)
{
  g_ptr_array_add (a, item);
}

static int compare_boxes (const void * p1, const void * p2)
{
  GfsBox * b1 = *(GfsBox **)p1;
  GfsBox * b2 = *(GfsBox **)p2;
  /* the check below is necessary when using graph partitioning */  
  if (GFS_IS_BOX (b1) && GFS_IS_BOX (b2))
    return b1->id < b2->id ? -1 : 1;
  else
    return 0;
}

static void domain_foreach (GtsContainer * c, 
			    GtsFunc func, 
			    gpointer data)
{
  GPtrArray * a = GFS_DOMAIN (c)->sorted;
  if (a == NULL) /* domain is being destroyed */
    (* GTS_CONTAINER_CLASS (GTS_OBJECT_CLASS (gfs_domain_class ())->parent_class)->foreach)
      (c, func, data);
  else {
    if (GFS_DOMAIN (c)->dirty) {
      g_ptr_array_set_size (a, 0);
      (* GTS_CONTAINER_CLASS (GTS_OBJECT_CLASS (gfs_domain_class ())->parent_class)->foreach)
	(c, (GtsFunc) add_item, a);
      qsort (a->pdata, a->len, sizeof (gpointer), compare_boxes);
      GFS_DOMAIN (c)->dirty = FALSE;
    }
    guint i;
    for (i = 0; i < a->len; i++)
      (* func) (a->pdata[i], data);
  }
}

static void domain_add (GtsContainer * c, GtsContainee * i)
{
  (* GTS_CONTAINER_CLASS (GTS_OBJECT_CLASS (gfs_domain_class ())->parent_class)->add) (c, i);
  GFS_DOMAIN (c)->dirty = TRUE;
}

static void domain_remove (GtsContainer * c, GtsContainee * i)
{
  (* GTS_CONTAINER_CLASS (GTS_OBJECT_CLASS (gfs_domain_class ())->parent_class)->remove) (c, i);
  GFS_DOMAIN (c)->dirty = TRUE;
}

static void domain_class_init (GfsDomainClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = domain_read;
  GTS_OBJECT_CLASS (klass)->write = domain_write;
  GTS_OBJECT_CLASS (klass)->destroy = domain_destroy;

  GTS_CONTAINER_CLASS (klass)->foreach = domain_foreach;
  GTS_CONTAINER_CLASS (klass)->add = domain_add;
  GTS_CONTAINER_CLASS (klass)->remove = domain_remove;

  klass->post_read = domain_post_read;
}

static void domain_init (GfsDomain * domain)
{
  domain->pid = -1;

#ifdef HAVE_MPI
  int size;

  MPI_Comm_size (MPI_COMM_WORLD, &size);
  if (size > 1)
    MPI_Comm_rank (MPI_COMM_WORLD, &domain->pid);
#endif /* HAVE_MPI */

  domain->clock = g_timer_new ();
  domain->timer = gfs_clock_new ();
  domain->timers = g_hash_table_new (g_str_hash, g_str_equal);

  gts_range_init (&domain->size);

  domain->profile_bc = FALSE;

  gts_range_init (&domain->mpi_messages);
  gts_range_init (&domain->mpi_wait);

  domain->rootlevel = 0;
  domain->refpos.x = domain->refpos.y = domain->refpos.z = 0.;
  domain->lambda.x = domain->lambda.y = domain->lambda.z = 1.;

  domain->allocated = g_array_new (FALSE, TRUE, sizeof (gboolean));
  domain->variables = NULL;

  domain->variables_io = NULL;
  domain->max_depth_write = -1;

  domain->cell_init = (FttCellInitFunc) gfs_cell_fine_init;
  domain->cell_init_data = domain;

  domain->version = atoi (GFS_BUILD_VERSION);

  domain->overlap = TRUE;

  domain->objects = g_hash_table_new (g_str_hash, g_str_equal);

  domain->np = 1;

  domain->sorted = g_ptr_array_new ();
  domain->dirty = TRUE;
  
  domain->projections = NULL;
}

GfsDomainClass * gfs_domain_class (void)
{
  static GfsDomainClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_domain_info = {
      "GfsDomain",
      sizeof (GfsDomain),
      sizeof (GfsDomainClass),
      (GtsObjectClassInitFunc) domain_class_init,
      (GtsObjectInitFunc) domain_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_wgraph_class ()),
				  &gfs_domain_info);
  }

  return klass;
}

typedef struct {
  FttTraverseFlags flags;
  gint max_depth;
  GfsVariable * v, * v1;
  FttComponent c;
  GfsLinearProblem * lp;
} BcData;

static void box_bc (GfsBox * box, BcData * p)
{
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      GfsBoundary * b = GFS_BOUNDARY (box->neighbor[d]);
      GfsBc * bc = gfs_boundary_lookup_bc (b, p->v);

      if (bc) {
	b->v = p->v1;
  	b->type = GFS_BOUNDARY_CENTER_VARIABLE;
	gfs_boundary_update (b);
	bc->v = p->v1;
	ftt_face_traverse_boundary (b->root, b->d,
				    FTT_PRE_ORDER, p->flags, p->max_depth,
				    bc->bc, bc);
	bc->v = p->v;
	gfs_boundary_send (b);
      }
    }
}

static void direction_face_bc (GtsObject * neighbor,
			       GfsVariable * v)
{
  if (GFS_IS_BOUNDARY (neighbor)) {
    GfsBoundary * b = GFS_BOUNDARY (neighbor);
    GfsBc * bc = gfs_boundary_lookup_bc (b, v);

    if (bc) {
      b->v = v;
      b->type = GFS_BOUNDARY_CENTER_VARIABLE;
      ftt_face_traverse_boundary (b->root, b->d,
				  FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				  bc->face_bc, bc);
      b->type = GFS_BOUNDARY_FACE_VARIABLE;
      gfs_boundary_send (b);
    }
  }
}

static void box_face_bc (GfsBox * box, BcData * p)
{
  if (p->c == FTT_XYZ) {
    FttDirection d;
    
    for (d = 0; d < FTT_NEIGHBORS; d++)
      direction_face_bc (box->neighbor[d], p->v);
  }
  else {
    direction_face_bc (box->neighbor[2*p->c], p->v);
    direction_face_bc (box->neighbor[2*p->c + 1], p->v);
  }
}

static void box_receive_bc (GfsBox * box, BcData * r)
{
  if (r->c == FTT_XYZ) {
    FttDirection d;
    
    for (d = 0; d < FTT_NEIGHBORS; d++) {
      FttDirection od = FTT_OPPOSITE_DIRECTION (d);

      if (GFS_IS_BOUNDARY (box->neighbor[od]))
	gfs_boundary_receive (GFS_BOUNDARY (box->neighbor[od]), r->flags, r->max_depth);
    }
  }
  else {
    if (GFS_IS_BOUNDARY (box->neighbor[2*r->c + 1]))
      gfs_boundary_receive (GFS_BOUNDARY (box->neighbor[2*r->c + 1]), r->flags, r->max_depth);
    if (GFS_IS_BOUNDARY (box->neighbor[2*r->c]))
      gfs_boundary_receive (GFS_BOUNDARY (box->neighbor[2*r->c]), r->flags, r->max_depth);
  }
}

static void box_match (GfsBox * box)
{
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      GfsBoundary * boundary = GFS_BOUNDARY (box->neighbor[d]);

      g_assert (GFS_BOUNDARY_CLASS (box->neighbor[d]->klass)->match);
      boundary->type = GFS_BOUNDARY_MATCH_VARIABLE;
      (* GFS_BOUNDARY_CLASS (box->neighbor[d]->klass)->match) (boundary);
      if (!boundary->root) {
	gts_object_destroy (GTS_OBJECT (boundary));
	box->neighbor[d] = NULL;
      }
      else
	gfs_boundary_send (boundary);
    }
}

static void box_synchronize (GfsBox * box, FttComponent * c)
{
  if (*c == FTT_XYZ) {
    FttDirection d;
    
    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (GFS_IS_BOUNDARY (box->neighbor[d]))
	gfs_boundary_synchronize (GFS_BOUNDARY (box->neighbor[d]));
  }
  else {
    if (GFS_IS_BOUNDARY (box->neighbor[2*(*c)]))
      gfs_boundary_synchronize (GFS_BOUNDARY (box->neighbor[2*(*c)]));
    if (GFS_IS_BOUNDARY (box->neighbor[2*(*c) + 1]))
      gfs_boundary_synchronize (GFS_BOUNDARY (box->neighbor[2*(*c) + 1]));
  }
}

/**
 * gfs_domain_copy_bc:
 * @domain: a #GfsDomain.
 * @flags: the traversal flags.
 * @max_depth: the maximum depth of the traversal.
 * @v: a #GfsVariable.
 * @v1: another #GfsVariable.
 *
 * Apply the boundary conditions of variable @v in @domain to variable @v1.
 */
void gfs_domain_copy_bc (GfsDomain * domain,
			 FttTraverseFlags flags,
			 gint max_depth,
			 GfsVariable * v,
			 GfsVariable * v1)
{
  BcData b = { flags, max_depth, v, v1, FTT_XYZ};

  g_return_if_fail (domain != NULL);
  g_return_if_fail (v != NULL);
  g_return_if_fail (v1 != NULL);

  if (domain->profile_bc)
    gfs_domain_timer_start (domain, "bc");

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_bc, &b);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_receive_bc, &b);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_synchronize, &b.c);

  if (domain->profile_bc)
    gfs_domain_timer_stop (domain, "bc");
}

/**
 * gfs_domain_bc:
 * @domain: a #GfsDomain.
 * @flags: the traversal flags.
 * @max_depth: the maximum depth of the traversal.
 * @v: a #GfsVariable.
 *
 * Apply the boundary conditions in @domain for variable @v.
 */
void gfs_domain_bc (GfsDomain * domain,
		    FttTraverseFlags flags,
		    gint max_depth,
		    GfsVariable * v)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (v != NULL);

  gfs_domain_copy_bc (domain, flags, max_depth, v, v);
}

static void box_homogeneous_bc (GfsBox * box, BcData * p)
{
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++) 
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      GfsBoundary * b = GFS_BOUNDARY (box->neighbor[d]);
      GfsBc * bc = gfs_boundary_lookup_bc (b, p->v);

      if (bc) {
	b->v = p->v1;
	bc->v = p->v1;
	b->type = GFS_BOUNDARY_CENTER_VARIABLE;
	ftt_face_traverse_boundary (b->root, b->d,
				    FTT_PRE_ORDER, p->flags, p->max_depth,
				    bc->homogeneous_bc, bc);
	bc->v = p->v;
	gfs_boundary_send (b);
      }
    }
}

static void box_homogeneous_bc_stencil (GfsBox * box, BcData * p)
{
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++) 
    if (GFS_IS_BOUNDARY (box->neighbor[d]) && !GFS_IS_BOUNDARY_PERIODIC (box->neighbor[d])) {
      GfsBoundary * b = GFS_BOUNDARY (box->neighbor[d]);
      GfsBc * bc = gfs_boundary_lookup_bc (b, p->v);
      
      if (bc) {
	b->v = p->v1;
	bc->v = p->v1;
	b->type = GFS_BOUNDARY_CENTER_VARIABLE;
	bc->lp = p->lp;
	ftt_face_traverse_boundary (b->root, b->d,
				    FTT_PRE_ORDER, p->flags, p->max_depth,
				    bc->homogeneous_bc_stencil, bc);
	bc->v = p->v;
	gfs_boundary_send (b);
      }
    }
}

/**
 * gfs_domain_homogeneous_bc:
 * @domain: a #GfsDomain.
 * @flags: the traversal flags.
 * @max_depth: the maximum depth of the traversal.
 * @ov: a #GfsVariable.
 * @v: a #GfsVariable of which @ov is an homogeneous version.
 *
 * Apply the boundary conditions in @domain for variable @ov using the
 * homogeneous version of the boundary condititons for @v.
 */
void gfs_domain_homogeneous_bc (GfsDomain * domain,
				FttTraverseFlags flags,
				gint max_depth,
				GfsVariable * ov,
				GfsVariable * v)
{
  BcData b = { flags, max_depth, v, ov, FTT_XYZ};

  g_return_if_fail (domain != NULL);
  g_return_if_fail (ov != NULL);
  g_return_if_fail (v != NULL);

  if (domain->profile_bc)
    gfs_domain_timer_start (domain, "bc");

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_homogeneous_bc, &b);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_receive_bc, &b);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_synchronize, &b.c);

  if (domain->profile_bc)
    gfs_domain_timer_stop (domain, "bc");
}

/**
 * gfs_domain_homogeneous_bc_stencil:
 * @domain: a #GfsDomain.
 * @flags: the traversal flags.
 * @max_depth: the maximum depth of the traversal.
 * @ov: a #GfsVariable.
 * @v: a #GfsVariable of which @ov is an homogeneous version.
 * @lp: a #GfsLinearProblem in which to store the stencil.
 *
 * Gets the stencils corresponding to the homogeneous boundary
 * conditions in @domain for variable @v.
 */
void gfs_domain_homogeneous_bc_stencil (GfsDomain * domain,
					FttTraverseFlags flags,
					gint max_depth,
					GfsVariable * ov,
					GfsVariable * v,
					GfsLinearProblem * lp)
{
  BcData b = { flags, max_depth, v, ov, FTT_XYZ, lp };

  g_return_if_fail (domain != NULL);
  g_return_if_fail (v != NULL);

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_homogeneous_bc_stencil, &b);
}

typedef struct {
  FttCellTraverseFunc func;
  gpointer data;
  FttTraverseType order;
  FttTraverseFlags flags;
  gint max_depth;
} TraverseData;

typedef struct {
  TraverseData t;
  BcData b;
} TraverseBcData;

static void update_mpi_cell (FttCell * cell, TraverseData * p)
{
  if ((cell->flags & GFS_FLAG_USED) == 0) {
    (* p->func) (cell, p->data);
    cell->flags |= GFS_FLAG_USED;
  }
}

static void update_other_cell (FttCell * cell, TraverseData * p)
{
  if ((cell->flags & GFS_FLAG_USED) != 0)
    cell->flags &= ~GFS_FLAG_USED;
  else
    (* p->func) (cell, p->data);
}

static void update_mpi_boundaries (GfsBox * box, TraverseBcData * p)
{
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++) 
    if (GFS_IS_BOUNDARY_MPI (box->neighbor[d])) {
      GfsBoundary * b = GFS_BOUNDARY (box->neighbor[d]);
      GfsBc * bc = gfs_boundary_lookup_bc (b, p->b.v);

      if (bc) {
	ftt_cell_traverse_boundary (box->root, d, p->t.order, p->t.flags, p->t.max_depth,
				    (FttCellTraverseFunc) update_mpi_cell, p);
	b->v = p->b.v1;
	bc->v = p->b.v1;
	b->type = GFS_BOUNDARY_CENTER_VARIABLE;
	ftt_face_traverse_boundary (b->root, b->d,
				    FTT_PRE_ORDER, p->b.flags, p->b.max_depth,
				    bc->bc, bc);
	bc->v = p->b.v;
	gfs_boundary_send (b);
      }
    }
}

static void update_other_homogeneous_boundaries (GfsBox * box, BcData * p)
{
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d]) &&
	!GFS_IS_BOUNDARY_MPI (box->neighbor[d])) {
      GfsBoundary * b = GFS_BOUNDARY (box->neighbor[d]);
      GfsBc * bc = gfs_boundary_lookup_bc (b, p->v);

      if (bc) {
	b->v = p->v1;
	bc->v = p->v1;
	b->type = GFS_BOUNDARY_CENTER_VARIABLE;
	ftt_face_traverse_boundary (b->root, b->d,
				    FTT_PRE_ORDER, p->flags, p->max_depth,
				    bc->homogeneous_bc, bc);
	bc->v = p->v;
	gfs_boundary_send (b);
      }
    }
}

/**
 * gfs_traverse_and_homogeneous_bc:
 * @domain: a #GfsDomain.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @max_depth: the maximum depth of the traversal. Cells below this
 * depth will not be traversed. If @max_depth is -1 all cells in the
 * tree are visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 * @ov: a #GfsVariable.
 * @v: a #GfsVariable of which @ov is an homogeneous version.
 *
 * For serial runs, this is identical to calling:
 *
 * gfs_domain_cell_traverse (domain, order, flags, max_depth, func, data);
 * gfs_domain_homogeneous_bc (domain, flags, max_depth, ov, v);
 *
 * For parallel runs, the communications needed to apply the boundary
 * conditions are overlapped with the calls to @func in the bulk of
 * the domain.
 */
void gfs_traverse_and_homogeneous_bc (GfsDomain * domain,
				      FttTraverseType order,
				      FttTraverseFlags flags,
				      gint max_depth,
				      FttCellTraverseFunc func,
				      gpointer data,
				      GfsVariable * ov,
				      GfsVariable * v)
{
  g_return_if_fail (domain != NULL);

  if (domain->pid < 0 || !domain->overlap) {
    gfs_domain_cell_traverse (domain, order, flags, max_depth, func, data);
    gfs_domain_homogeneous_bc (domain, flags, max_depth, ov, v);
  }
  else {
    TraverseBcData d = {
      { func, data, order, flags, max_depth }, 
      { flags, max_depth, v, ov, FTT_XYZ }
    };
    /* Update and send MPI boundary values */
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) update_mpi_boundaries, &d);
    /* Update bulk of domain and other boundaries */
    gfs_domain_cell_traverse (domain, order, flags, max_depth, 
			      (FttCellTraverseFunc) update_other_cell, &d);
    /* Apply homogeneous BC on other boundaries */
    gts_container_foreach (GTS_CONTAINER (domain), 
			   (GtsFunc) update_other_homogeneous_boundaries, &d.b);
    /* Receive and synchronize */
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_receive_bc, &d.b);
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_synchronize, &d.b.c);
  }
}

static void update_other_boundaries (GfsBox * box, BcData * p)
{
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d]) &&
	!GFS_IS_BOUNDARY_MPI (box->neighbor[d])) {
      GfsBoundary * b = GFS_BOUNDARY (box->neighbor[d]);
      GfsBc * bc = gfs_boundary_lookup_bc (b, p->v);

      if (bc) {
	b->v = p->v1;
	bc->v = p->v1;
	b->type = GFS_BOUNDARY_CENTER_VARIABLE;
	ftt_face_traverse_boundary (b->root, b->d,
				    FTT_PRE_ORDER, p->flags, p->max_depth,
				    bc->bc, bc);
	bc->v = p->v;
	gfs_boundary_send (b);
      }
    }
}

/**
 * gfs_traverse_and_bc:
 * @domain: a #GfsDomain.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @max_depth: the maximum depth of the traversal. Cells below this
 * depth will not be traversed. If @max_depth is -1 all cells in the
 * tree are visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 * @v: a #GfsVariable.
 * @v1: another #GfsVariable.
 *
 * For serial runs, this is identical to calling:
 *
 * gfs_domain_cell_traverse (domain, order, flags, max_depth, func, data);
 * gfs_domain_copy_bc (domain, flags, max_depth, v, v1);
 *
 * For parallel runs, the communications needed to apply the boundary
 * conditions are overlapped with the calls to @func in the bulk of
 * the domain.
 */
void gfs_traverse_and_bc (GfsDomain * domain,
			  FttTraverseType order,
			  FttTraverseFlags flags,
			  gint max_depth,
			  FttCellTraverseFunc func,
			  gpointer data,
			  GfsVariable * v,
			  GfsVariable * v1)
{
  g_return_if_fail (domain != NULL);

  if (domain->pid < 0 || !domain->overlap) {
    gfs_domain_cell_traverse (domain, order, flags, max_depth, func, data);
    gfs_domain_copy_bc (domain, flags, max_depth, v, v1);
  }
  else {
    TraverseBcData d = {
      { func, data, order, flags, max_depth },
      { flags, max_depth, v, v1, FTT_XYZ }
    };
    /* Update and send MPI boundary values */
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) update_mpi_boundaries, &d);
    /* Update bulk of domain and other boundaries */
    gfs_domain_cell_traverse (domain, order, flags, max_depth, 
    			      (FttCellTraverseFunc) update_other_cell, &d);
    /* Apply BC on other boundaries */
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) update_other_boundaries, &d.b);
    /* Receive and synchronize */
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_receive_bc, &d.b);
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_synchronize, &d.b.c);
  }
}

/**
 * gfs_domain_face_bc:
 * @domain: a #GfsDomain.
 * @c: a component.
 * @v: a #GfsVariable.
 *
 * Apply the boundary conditions on the faces of @domain for variable @v.
 */
void gfs_domain_face_bc (GfsDomain * domain,
			 FttComponent c,
			 GfsVariable * v)
{
  BcData b = { FTT_TRAVERSE_LEAFS, -1, v, v, c };

  g_return_if_fail (domain != NULL);
  g_return_if_fail (c == FTT_XYZ || (c >= 0 && c < FTT_DIMENSION));
  g_return_if_fail (v != NULL);

  if (domain->profile_bc)
    gfs_domain_timer_start (domain, "face_bc");

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_face_bc, &b);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_receive_bc, &b);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_synchronize, &b.c);

  if (domain->profile_bc)
    gfs_domain_timer_stop (domain, "face_bc");
}

static void box_changed (GfsBox * box, gboolean * changed)
{
  FttDirection d;

  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d]))
      *changed |= GFS_BOUNDARY (box->neighbor[d])->changed;
}

static void refine_cell_corner (FttCell * cell, GfsDomain * domain)
{
  if (FTT_CELL_IS_LEAF (cell) && ftt_refine_corner (cell))
    ftt_cell_refine_single (cell, domain->cell_init, domain->cell_init_data);
}

static void box_depth (GfsBox * box, guint * depth)
{
  guint d = ftt_cell_depth (box->root);

  if (d > *depth)
    *depth = d;
}

static gboolean domain_match (GfsDomain * domain)
{
  BcData b = { FTT_TRAVERSE_LEAFS, -1, NULL, NULL, FTT_XYZ };
  gboolean changed = FALSE;
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_match, NULL);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_receive_bc, &b);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_synchronize, &b.c);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_changed, &changed);
  if (changed) {
    gint l;
    guint depth = 0;
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_depth, &depth);
    for (l = depth - 2; l >= 0; l--)
      gfs_domain_cell_traverse (domain,
				FTT_PRE_ORDER, FTT_TRAVERSE_LEVEL, l,
				(FttCellTraverseFunc) refine_cell_corner, domain);
  }
  gfs_all_reduce (domain, changed, MPI_INT, MPI_MAX);
  return changed;
}

/**
 * gfs_domain_match:
 * @domain: a #GfsDomain.
 *
 * Match the boundaries of @domain.
 */
void gfs_domain_match (GfsDomain * domain)
{
  g_return_if_fail (domain != NULL);

  if (domain->profile_bc)
    gfs_domain_timer_start (domain, "match");

  while (domain_match (domain));

  if (domain->profile_bc)
    gfs_domain_timer_stop (domain, "match");
}

/**
 * gfs_domain_forget_boundary:
 * @domain: a #GfsDomain.
 * @boundary: a #GfsBoundary belonging to @domain.
 *
 * Makes @domain permanently ignore @boundary when performing
 * "locate()" queries.
 */
void gfs_domain_forget_boundary (GfsDomain * domain, GfsBoundary * boundary)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (boundary != NULL);
  g_return_if_fail (gfs_box_domain (boundary->box) == domain);

  GfsLocateArray * a = domain->array;
  if (a == NULL)
    return;

  gint i;
  for (i = 0; i < a->size; i++)
    a->root[i] = g_slist_remove (a->root[i], boundary);
}

static void dirichlet_bc (FttCell * cell)
{
  cell->flags |= GFS_FLAG_DIRICHLET;
  GFS_STATE (cell)->solid->fv = 0.;
}

static void neumann_bc (FttCell * cell)
{
  cell->flags &= ~GFS_FLAG_DIRICHLET;
  GFS_STATE (cell)->solid->fv = 0.;
}

static gboolean is_velocity (GfsVariable * v, GfsDomain * domain)
{
  FttComponent c;
  GfsVariable ** u = gfs_domain_velocity (domain);

  for (c = 0; c < FTT_DIMENSION; c++)
    if (v == u[c])
      return TRUE;
  return FALSE;
}

/**
 * gfs_domain_surface_bc:
 * @domain: a #GfsDomain.
 * @v: a #GfsVariable.
 *
 * Apply boundary conditions for variable @v on embedded surfaces. 
 */
void gfs_domain_surface_bc (GfsDomain * domain,
			    GfsVariable * v)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (v != NULL);

  if (v->surface_bc)
    gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL,
      (FttCellTraverseFunc) GFS_SURFACE_GENERIC_BC_CLASS (GTS_OBJECT (v->surface_bc)->klass)->bc, 
			       v->surface_bc);
  else if (is_velocity (v, domain))
    gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL,
			       (FttCellTraverseFunc) dirichlet_bc, NULL);
  else
    gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL,
			       (FttCellTraverseFunc) neumann_bc, NULL);
}

static void box_traverse (GfsBox * box, TraverseData * d)
{
  ftt_cell_traverse (box->root, d->order, d->flags, d->max_depth, d->func, d->data);
}

/**
 * gfs_domain_cell_traverse:
 * @domain: a #GfsDomain.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @max_depth: the maximum depth of the traversal. Cells below this
 * depth will not be traversed. If @max_depth is -1 all cells in the
 * tree are visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 *
 * Traverses the cell trees of @domain. Calls the given function for
 * each cell visited.  
 */
void gfs_domain_cell_traverse (GfsDomain * domain,
			       FttTraverseType order,
			       FttTraverseFlags flags,
			       gint max_depth,
			       FttCellTraverseFunc func,
			       gpointer data)
{
  TraverseData d = { func, data, order, flags, max_depth };

  g_return_if_fail (domain != NULL);
  g_return_if_fail (func != NULL);

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_traverse, &d);
}

static void cell_traverse_add (FttCell * cell, GPtrArray * a)
{
  g_ptr_array_add (a, cell);
}

/**
 * gfs_domain_cell_traverse_new:
 * @domain: a #GfsDomain.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @max_depth: the maximum depth of the traversal. Cells below this
 * depth will not be traversed. If @max_depth is -1 all cells in the
 * tree are visited.
 *
 * Returns: a new #FttCellTraverse.
 */
FttCellTraverse * gfs_domain_cell_traverse_new (GfsDomain * domain,
						FttTraverseType order,
						FttTraverseFlags flags,
						gint max_depth)
{
  g_return_val_if_fail (domain != NULL, NULL);

  GPtrArray * a = g_ptr_array_new ();
  gfs_domain_cell_traverse (domain, order, flags, max_depth,
			    (FttCellTraverseFunc) cell_traverse_add, a);
  g_ptr_array_add (a, NULL);
  FttCellTraverse * t = g_malloc (sizeof (FttCellTraverse));
  t->current = t->cells = (FttCell **) a->pdata;
  g_ptr_array_free (a, FALSE);
  return t;
}

/**
 * gfs_domain_traverse_layers:
 * @domain: a #GfsDomain.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 *
 * Traverses the leaf cells of @domain for each layer. Calls the given
 * function for each cell visited.
 *
 * By default it is identical to gfs_domain_traverse_leaves() but can
 * be overloaded for specific (layered) domains.
 */
void gfs_domain_traverse_layers (GfsDomain * domain,
				 FttCellTraverseFunc func,
				 gpointer data)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (func != NULL);

  if (domain->traverse_layers)
    (* domain->traverse_layers) (domain, func, data);
  else
    gfs_domain_traverse_leaves (domain, func, data);
}

static void box_traverse_box (GfsBox * box, gpointer * datum)
{
  FttTraverseType * order = datum[0];
  FttTraverseFlags * flags = datum[1];
  gint * max_depth = datum[2];
  FttCellTraverseFunc func = (FttCellTraverseFunc) datum[3];
  gpointer data = datum[4];
  GtsBBox * bb = datum[5];

  ftt_cell_traverse_box (box->root, bb, 
			 *order, *flags, *max_depth, func, data);
}

/**
 * gfs_domain_cell_traverse_box:
 * @domain: a #GfsDomain.
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
 * Traverses the cell trees of @domain. Calls the given function for
 * each cell visited. Only the cells overlapping with @box are visited.
 */
void gfs_domain_cell_traverse_box (GfsDomain * domain,
				   GtsBBox * box,
				   FttTraverseType order,
				   FttTraverseFlags flags,
				   gint max_depth,
				   FttCellTraverseFunc func,
				   gpointer data)
{
  gpointer datum[6];

  datum[0] = &order;
  datum[1] = &flags;
  datum[2] = &max_depth;
  datum[3] = func;
  datum[4] = data;
  datum[5] = box;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (box != NULL);
  g_return_if_fail (func != NULL);

  gts_container_foreach (GTS_CONTAINER (domain), 
			 (GtsFunc) box_traverse_box, datum);
}

static void box_traverse_condition (GfsBox * box, gpointer * datum)
{
  FttTraverseType * order = datum[0];
  FttTraverseFlags * flags = datum[1];
  gint * max_depth = datum[2];
  FttCellTraverseFunc func = (FttCellTraverseFunc) datum[3];
  gpointer data = datum[4];
  gboolean (* condition) (FttCell *, gpointer) = datum[5];
  gpointer cdata = datum[6];

  ftt_cell_traverse_condition (box->root, *order, *flags, *max_depth, func, data,
			       condition, cdata);
}

/**
 * gfs_domain_cell_traverse_condition:
 * @domain: a #GfsDomain.
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
 * Traverses the cell trees of @domain. Calls the given function for
 * each cell visited.
 *
 * Traversal of any branch of the tree is stopped whenever @condition
 * is not verified.
 */
void gfs_domain_cell_traverse_condition (GfsDomain * domain,
					 FttTraverseType order,
					 FttTraverseFlags flags,
					 gint max_depth,
					 FttCellTraverseFunc func,
					 gpointer data,
					 gboolean (* condition) (FttCell *, gpointer),
					 gpointer cdata)
{
  gpointer datum[7];

  datum[0] = &order;
  datum[1] = &flags;
  datum[2] = &max_depth;
  datum[3] = func;
  datum[4] = data;
  datum[5] = condition;
  datum[6] = cdata;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (func != NULL);
  g_return_if_fail (condition != NULL);

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_traverse_condition, datum);
}

static void traverse_mixed (GfsBox * box, TraverseData * d)
{
  gfs_cell_traverse_mixed (box->root, d->order, d->flags, d->func, d->data);
}

/**
 * gfs_domain_traverse_mixed:
 * @domain: a #GfsDomain.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 *
 * Calls @func for each mixed cell of @domain.
 */
void gfs_domain_traverse_mixed (GfsDomain * domain,
				FttTraverseType order,
				FttTraverseFlags flags,
				FttCellTraverseFunc func,
				gpointer data)
{
  TraverseData d = { func, data, order, flags, -1 };

  g_return_if_fail (domain != NULL);
  g_return_if_fail (func != NULL);

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) traverse_mixed, &d);
}

typedef struct {
  FttCellTraverseCutFunc func;
  gpointer data;
  FttTraverseType order;
  FttTraverseFlags flags;
  GfsGenericSurface * s;
} TraverseCut;

static void traverse_cut (GfsBox * box, TraverseCut * p)
{
  gfs_cell_traverse_cut (box->root, p->s, p->order, p->flags, p->func, p->data);
}

/**
 * gfs_domain_traverse_cut:
 * @domain: a #GfsDomain.
 * @s: a #GfsGenericSurface.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 *
 * Calls @func for each cell of @domain cut by @s.
 */
void gfs_domain_traverse_cut (GfsDomain * domain,
			      GfsGenericSurface * s,
			      FttTraverseType order,
			      FttTraverseFlags flags,
			      FttCellTraverseCutFunc func,
			      gpointer data)
{
  TraverseCut p = { func, data, order, flags, s };

  g_return_if_fail (domain != NULL);
  g_return_if_fail (s != NULL);
  g_return_if_fail (func != NULL);

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) traverse_cut, &p);
}

static void traverse_cut_2D (GfsBox * box, TraverseCut * p)
{
  gfs_cell_traverse_cut_2D (box->root, p->s, p->order, p->flags, p->func, p->data);
}

/**
 * gfs_domain_traverse_cut_2D:
 * @domain: a #GfsDomain.
 * @s: a #GfsGenericSurface.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 *
 * Calls @func for each cell of @domain cut by @s.
 *
 * The cells are flattened in the z-direction.
 */
void gfs_domain_traverse_cut_2D (GfsDomain * domain,
				 GfsGenericSurface * s,
				 FttTraverseType order,
				 FttTraverseFlags flags,
				 FttCellTraverseCutFunc func,
				 gpointer data)
{
  TraverseCut p = { func, data, order, flags, s };

  g_return_if_fail (domain != NULL);
  g_return_if_fail (s != NULL);
  g_return_if_fail (func != NULL);

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) traverse_cut_2D, &p);
}

/**
 * gfs_domain_depth:
 * @domain: a #GfsDomain.
 *
 * Returns: the maximum depth of the cell trees of @domain. This
 * function is global i.e. it returns the maximum depth over all the
 * processes (for parallel execution).
 */
guint gfs_domain_depth (GfsDomain * domain)
{
  guint depth = 0;

  g_return_val_if_fail (domain != NULL, 0);

  gts_container_foreach (GTS_CONTAINER (domain),
			 (GtsFunc) box_depth, &depth);
  gfs_all_reduce (domain, depth, MPI_UNSIGNED, MPI_MAX);
  return depth;
}

#include "ftt_internal.c"

/**
 * gfs_domain_face_traverse:
 * @domain: a #GfsDomain.
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
 * Traverses a @domain. Calls the given function for each face
 * of the cells of the domain.
 *
 * If %FTT_TRAVERSE_BOUNDARY_FACES is not set in @flags, only
 * "double-sided" faces are traversed i.e. the @neighbor field of the
 * face is never %NULL.  
 */
void gfs_domain_face_traverse (GfsDomain * domain,
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
  
  g_return_if_fail (domain != NULL);
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
      gfs_domain_cell_traverse (domain, order, flags, max_depth, 
	  (FttCellTraverseFunc) traverse_all_faces, 
				datum);
    }
    else {
      gfs_domain_cell_traverse (domain, order, flags, max_depth, 
	  (FttCellTraverseFunc) traverse_all_direct_faces, 
				datum);
      datum[0] = &d;
      check = TRUE;
      for (d = 1; d < FTT_NEIGHBORS; d += 2)
	gfs_domain_cell_traverse_boundary (domain, 
					   d, order, flags, max_depth, 
					   (FttCellTraverseFunc) traverse_face, datum);
    }
  }
  else if (c == FTT_XY) {
    gfs_domain_face_traverse (domain, FTT_X, order, flags, max_depth, func, data);
    gfs_domain_face_traverse (domain, FTT_Y, order, flags, max_depth, func, data);
  }
  else {
    if (boundary_faces) {
      check = TRUE;
      datum[0] = &c;
      gfs_domain_cell_traverse (domain, order, flags, max_depth, 
				(FttCellTraverseFunc) traverse_face_component,
				datum);
    }
    else {
      d = 2*c;
      datum[0] = &d;
      gfs_domain_cell_traverse (domain, order, flags, max_depth, 
				(FttCellTraverseFunc) traverse_face_direction, 
				datum);
      d = 2*c + 1;
      check = TRUE;
      gfs_domain_cell_traverse_boundary (domain, d, order, flags, max_depth, 
					 (FttCellTraverseFunc) traverse_face, datum);
    }
  }
  gfs_domain_cell_traverse (domain, order, flags, max_depth, 
			    (FttCellTraverseFunc) reset_flag, NULL);
}

static void cell_traverse_boundary (GfsBox * box, gpointer * datum)
{
  FttDirection * d = datum[0];

  if (!GFS_IS_BOX (box->neighbor[*d])) {
    FttTraverseType * order = datum[1];
    FttTraverseFlags * flags = datum[2];
    gint * max_depth = datum[3];
    FttCellTraverseFunc func = (FttCellTraverseFunc) datum[4];
    gpointer data = datum[5];

    ftt_cell_traverse_boundary (box->root, 
				*d, *order, *flags, *max_depth, func, data);
  }
}

/**
 * gfs_domain_cell_traverse_boundary:
 * @domain: a #GfsDomain.
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
 * Traverses the boundary of a domain in direction @d. Calls the given
 * function for each cell visited.  
 */
void gfs_domain_cell_traverse_boundary (GfsDomain * domain,
					FttDirection d,
					FttTraverseType order,
					FttTraverseFlags flags,
					gint max_depth,
					FttCellTraverseFunc func,
					gpointer data)
{
  gpointer datum[6];
  
  datum[0] = &d;
  datum[1] = &order;
  datum[2] = &flags;
  datum[3] = &max_depth;
  datum[4] = func;
  datum[5] = data;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (d < FTT_NEIGHBORS);
  g_return_if_fail (func != NULL);

  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) cell_traverse_boundary, datum);
}

static void add_stats (const FttCell * cell, gpointer * data)
{
  GtsRange * s = data[0];
  gdouble v = GFS_VALUE (cell, GFS_VARIABLE (data[1]));

  if (v != GFS_NODATA)
    gts_range_add_value (s, v);
}

#ifdef HAVE_MPI
static void range_reduce (void * i, void * o, 
			  int * len,
			  MPI_Datatype * type)
{
  gdouble * in = (gdouble *) i;
  gdouble * inout = (gdouble *) o;
  g_assert (*len == 5);
  
  if (in[0] < inout[0]) /* min */
    inout[0] = in[0];
  if (in[1] > inout[1]) /* max */
    inout[1] = in[1];
  inout[2] += in[2];    /* sum */
  inout[3] += in[3];    /* sum2 */
  inout[4] += in[4];    /* n */
}

static void domain_range_reduce (GfsDomain * domain, GtsRange * s)
{
  if (domain->pid >= 0) {
    double in[5];
    double out[5] = { G_MAXDOUBLE, - G_MAXDOUBLE, 0., 0., 0. };
    MPI_Op op;
    
    MPI_Op_create (range_reduce, TRUE, &op);
    in[0] = s->min; in[1] = s->max; in[2] = s->sum; in[3] = s->sum2;
    in[4] = s->n;
    MPI_Allreduce (in, out, 5, MPI_DOUBLE, op, MPI_COMM_WORLD);
    MPI_Op_free (&op);
    s->min = out[0]; s->max = out[1]; s->sum = out[2]; s->sum2 = out[3];
    s->n = out[4];
  }
}
#else /* not HAVE_MPI */
static void domain_range_reduce (GfsDomain * domain, GtsRange * s)
{
}
#endif /* HAVE_MPI */

/**
 * gfs_domain_stats_variable:
 * @domain: the domain to obtain statistics from.
 * @v: a #GfsVariable.
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 * @condition: a condition or %NULL.
 * @cdata: user data to pass to @condition.
 *
 * Traverses the domain defined by @domain using
 * gfs_domain_cell_traverse() and gathers statistics about variable
 * @v.
 *
 * Only cells veryfing @condition are taken into account (if
 * @condition is not %NULL). See also
 * gfs_domain_cell_traverse_condition().
 *
 * Returns: a #GtsRange containing the statistics about @v.
 */
GtsRange gfs_domain_stats_variable (GfsDomain * domain,
				    GfsVariable * v,
				    FttTraverseFlags flags,
				    gint max_depth,
				    gboolean (* condition) (FttCell *, gpointer),
				    gpointer cdata)
{
  GtsRange s;
  gpointer data[2];

  g_return_val_if_fail (domain != NULL, s);
  g_return_val_if_fail (v != NULL, s);

  gts_range_init (&s);
  data[0] = &s;
  data[1] = v;
  if (condition)
    gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, flags, max_depth, 
					(FttCellTraverseFunc) add_stats, data,
					condition, cdata);
  else
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth, 
			      (FttCellTraverseFunc) add_stats, data);
  domain_range_reduce (domain, &s);
  gts_range_update (&s);

  return s;
}

static void add_stats_solid (FttCell * cell, GtsRange * s)
{
  gts_range_add_value (s, GFS_STATE (cell)->solid->a);
}

/**
 * gfs_domain_stats_solid:
 * @domain: the domain to obtain statistics from.
 *
 * Traverses the domain defined by @domain using gfs_domain_traverse_mixed()
 * and gathers statistics about the solid volume fraction in mixed cells.
 *
 * Returns: statistics about the solid volume fraction @a in mixed cells.
 */
GtsRange gfs_domain_stats_solid (GfsDomain * domain)
{
  GtsRange s;

  g_return_val_if_fail (domain != NULL, s);

  gts_range_init (&s);
  gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			    (FttCellTraverseFunc) add_stats_solid, &s);
  domain_range_reduce (domain, &s);
  gts_range_update (&s);

  return s;
}

static void add_stats_merged (GSList * m, gpointer * data)
{
  GtsRange * solid =  data[0];
  GtsRange * number = data[1];
  gdouble a = 0.;
  guint n = 0;

  while (m) {
    FttCell * c = m->data;

    a += GFS_IS_MIXED (c) ? GFS_STATE (c)->solid->a : 1.;
    n++;
    m = m->next;
  }
  if (n > 1 || a < 1.)
    gts_range_add_value (solid, a);
  if (n > 1)
    gts_range_add_value (number, n);
}

/**
 * gfs_domain_stats_merged:
 * @domain: the domain to obtain statistics from.
 * @solid: #GtsRange in which to return stats for the total solid
 * volume fraction of merged cells. 
 * @number: #GtsRange in which to return stats for the number of cells
 * used per merged cell.
 *
 * Traverses the domain defined by @domain using
 * gfs_domain_traverse_merged() and gathers statistics about the total
 * solid volume fraction of merged cells and the number of cells used
 * per merged cell.
 */
void gfs_domain_stats_merged (GfsDomain * domain,
			     GtsRange * solid,
			     GtsRange * number)
{
  gpointer data[2];

  g_return_if_fail (domain != NULL);
  g_return_if_fail (solid != NULL);
  g_return_if_fail (number != NULL);

  gts_range_init (solid);
  gts_range_init (number);
  data[0] = solid;
  data[1] = number;
  gfs_domain_traverse_merged (domain,
			     (GfsMergedTraverseFunc) add_stats_merged, data);
  domain_range_reduce (domain, solid);
  domain_range_reduce (domain, number);
  gts_range_update (solid);
  gts_range_update (number);
}

static void cell_count (FttCell * cell, guint * count)
{
  (*count)++;
}

#define BPID(b) (gfs_box_domain (b)->pid >= 0 && (b)->pid > 0 ? (b)->pid : 0)

static void box_count (GfsBox * b, GArray * a)
{
  guint count = 0, pid = BPID(b);
  ftt_cell_traverse (b->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
		     (FttCellTraverseFunc) cell_count, &count);
  if (pid >= a->len)
    g_array_set_size (a, pid + 1);
  g_array_index (a, guint, pid) += count;
}

static void boundary_size (GfsBox * box, GArray * a)
{
  FttDirection d;
  guint count = 0;

  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY_MPI (box->neighbor[d]) ||
	(GFS_IS_BOX (box->neighbor[d]) && GFS_BOX (box->neighbor[d])->pid != box->pid))
      ftt_cell_traverse_boundary (box->root, d, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				  (FttCellTraverseFunc) cell_count, &count);
  g_array_index (a, guint, BPID (box)) += count;
}

/**
 * gfs_domain_stats_balance:
 * @domain: the domain to obtain statistics from.
 * @size: #GtsRange in which to return stats for the total size of the domain.
 * @boundary: #GtsRange in which to return stats for the size of the parallel 
 * boundaries of the domain.
 * @mpiwait:  #GtsRange in which to return stats for the average time spend
 * waiting for MPI calls in each PE.
 *
 * Gathers statistics about the sizes of the domains, their parallel
 * boundaries and the execution time on each PE.  
 */
void gfs_domain_stats_balance (GfsDomain * domain,
			       GtsRange * size,
			       GtsRange * boundary,
			       GtsRange * mpiwait)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (size != NULL);
  g_return_if_fail (boundary != NULL);
  g_return_if_fail (mpiwait != NULL);

  gts_range_init (size);
  gts_range_init (boundary);
  gts_range_init (mpiwait);

  if (domain->timestep.n > 0)
    gts_range_add_value (mpiwait, domain->mpi_wait.sum/domain->timestep.n);

  GArray * a = g_array_new (FALSE, TRUE, sizeof (guint));
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_count, a);
  guint i;
  for (i = 0; i < a->len; i++) {
    guint v = g_array_index (a, guint, i);
    if (v > 0) {
      gts_range_add_value (size, v);
      g_array_index (a, guint, i) = 0;
    }
  }
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) boundary_size, a);
  for (i = 0; i < a->len; i++) {
    guint v = g_array_index (a, guint, i);
    if (v > 0)
      gts_range_add_value (boundary, v);
  }
  domain_range_reduce (domain, size);
  domain_range_reduce (domain, boundary);
  domain_range_reduce (domain, mpiwait);
  g_array_free (a, TRUE);
  gts_range_update (size);
  gts_range_update (boundary);
  gts_range_update (mpiwait);
}

static void add_norm (const FttCell * cell, gpointer * data)
{
  GfsNorm * n = data[0];
  GfsVariable * v = data[1];

  gfs_norm_add (n, GFS_VALUE (cell, v), gfs_cell_volume (cell, v->domain));
}

static void add_norm_weighted (FttCell * cell, gpointer * data)
{
  GfsNorm * n = data[0];
  GfsVariable * v = data[1];
  GfsFunction * w = data[2];

  gfs_norm_add (n, GFS_VALUE (cell, v), gfs_cell_volume (cell, v->domain)*
		gfs_function_value (w, cell));
}

#ifdef HAVE_MPI
static void norm_reduce (void * i, void * o, 
			 int * len,
			 MPI_Datatype * type)
{
  gdouble * in = (gdouble *) i;
  gdouble * inout = (gdouble *) o;
  g_assert (*len == 5);
  
  inout[0] += in[0];    /* bias */
  inout[1] += in[1];    /* first */
  inout[2] += in[2];    /* second */
  if (in[3] > inout[3]) /* infty */
    inout[3] = in[3];    
  inout[4] += in[4];    /* w */
}

static void domain_norm_reduce (GfsDomain * domain, GfsNorm * n)
{
  if (domain->pid >= 0) {
    double in[5];
    double out[5] = { 0., 0., 0., - G_MAXDOUBLE, 0. };
    MPI_Op op;

    MPI_Op_create (norm_reduce, TRUE, &op);
    in[0] = n->bias; in[1] = n->first; in[2] = n->second; in[3] = n->infty;
    in[4] = n->w;
    MPI_Allreduce (in, out, 5, MPI_DOUBLE, op, MPI_COMM_WORLD);
    MPI_Op_free (&op);
    n->bias = out[0]; n->first = out[1]; n->second = out[2]; n->infty = out[3];
    n->w = out[4];
  }
}
#else /* not HAVE_MPI */
static void domain_norm_reduce (GfsDomain * domain, GfsNorm * n)
{
}
#endif /* not HAVE_MPI */

/**
 * gfs_domain_norm_variable:
 * @domain: the domain to obtain norm from.
 * @v: a #GfsVariable.
 * @w: a #GfsFunction or %NULL.
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 * @condition: a condition or %NULL.
 * @cdata: user data to pass to @condition.
 *
 * Traverses the domain defined by @domain using gfs_domain_cell_traverse()
 * and gathers norm statistics about variable @v.
 *
 * The norm is weighted by the volume of each cell times the value of
 * function @w (if @w is not %NULL).
 *
 * Only cells veryfing @condition are taken into account (if
 * @condition is not %NULL). See also
 * gfs_domain_cell_traverse_condition().
 *
 * Returns: a #GfsNorm containing the norm statistics about @v.
 */
GfsNorm gfs_domain_norm_variable (GfsDomain * domain,
				  GfsVariable * v,
				  GfsFunction * w,
				  FttTraverseFlags flags,
				  gint max_depth,
				  gboolean (* condition) (FttCell *, gpointer),
				  gpointer cdata)
{
  GfsNorm n;
  gpointer data[3];

  g_return_val_if_fail (domain != NULL, n);
  g_return_val_if_fail (v != NULL, n);
  
  gfs_norm_init (&n);
  data[0] = &n;
  data[1] = v;
  data[2] = w;
  FttCellTraverseFunc func = w != NULL ?
    (FttCellTraverseFunc) add_norm_weighted : 
    (FttCellTraverseFunc) add_norm;
  if (w)
    gfs_catch_floating_point_exceptions ();
  if (condition)
    gfs_domain_cell_traverse_condition (domain, FTT_PRE_ORDER, flags, max_depth, 
					(FttCellTraverseFunc) func, data,
					condition, cdata);
  else
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth, 
			      (FttCellTraverseFunc) func, data);
  if (w)
    gfs_restore_fpe_for_function (w);
  domain_norm_reduce (domain, &n);
  gfs_norm_update (&n);

  return n;
}

typedef struct {
  GfsVariable * res;
  gdouble bias;
  GfsNorm n;
} ResData;

static void add_norm_residual (const FttCell * cell, ResData * p)
{
  GfsDomain * domain = p->res->domain;
  gdouble a = domain->cell_metric ? (* domain->cell_metric) (domain, cell) : 1.;
  gdouble size = ftt_cell_size (cell);
  gfs_norm_add (&p->n, GFS_VALUE (cell, p->res)/(a*size*size), 1.);
  p->bias += GFS_VALUE (cell, p->res);
}

/**
 * gfs_domain_norm_residual:
 * @domain: the domain to obtain the norm from.
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 * @dt: the time step.
 * @res: the residual.
 *
 * Traverses the domain defined by @domain using gfs_domain_cell_traverse()
 * and gathers norm statistics about the volume weighted relative residual
 * (i.e. the sum of the residual over the volume defined by each cell
 * divided by the total volume of the cell).
 *
 * Returns: a #GfsNorm containing the norm statistics about the volume
 * weighted relative residual.  
 */
GfsNorm gfs_domain_norm_residual (GfsDomain * domain,
				  FttTraverseFlags flags,
				  gint max_depth,
				  gdouble dt,
				  GfsVariable * res)
{
  ResData p = { res, 0. };

  g_return_val_if_fail (domain != NULL, p.n);
  g_return_val_if_fail (res != NULL, p.n);
  
  gfs_norm_init (&p.n);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth, 
			   (FttCellTraverseFunc) add_norm_residual, &p);
  domain_norm_reduce (domain, &p.n);
  gfs_all_reduce (domain, p.bias, MPI_DOUBLE, MPI_SUM);
  gfs_norm_update (&p.n);

  dt *= dt;
  p.n.bias = p.bias*dt;
  p.n.first *= dt;
  p.n.second *= dt;
  p.n.infty *= dt;
  return p.n;
}

/**
 * gfs_domain_velocity:
 * @domain: a #GfsDomain.
 *
 * Returns: the components of the velocity vector for @domain or %NULL.
 */
GfsVariable ** gfs_domain_velocity (GfsDomain * domain)
{
  FttComponent c;
  static gchar name[][2] = {"U","V","W"};

  g_return_val_if_fail (domain != NULL, NULL);
  
  for (c = 0; c < FTT_DIMENSION; c++) {
    GfsVariable * v = gfs_variable_from_name (domain->variables, name[c]);
    if (v == NULL)
      return NULL;
    domain->velocity[c] = v;
  }
  return domain->velocity;
}

static void add_norm_velocity (FttCell * cell, gpointer * data)
{
  GfsVariable ** u = data[0];
  GfsNorm * n = data[1];
  
  gfs_norm_add (n, gfs_vector_norm (cell, u), gfs_cell_volume (cell, u[0]->domain));
}

/**
 * gfs_domain_norm_velocity:
 * @domain: the domain to obtain the norm from.
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 *
 * Traverses the domain defined by @domain using gfs_domain_cell_traverse()
 * and gathers norm statistics about velocity.
 *
 * Returns: a #GfsNorm containing the norm statistics about the velocity.
 */
GfsNorm gfs_domain_norm_velocity (GfsDomain * domain,
				  FttTraverseFlags flags,
				  gint max_depth)
{
  GfsNorm n;
  gpointer data[2];

  g_return_val_if_fail (domain != NULL, n);
  
  gfs_norm_init (&n);
  data[0] = gfs_domain_velocity (domain);
  data[1] = &n;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth, 
			   (FttCellTraverseFunc) add_norm_velocity, data);
  domain_norm_reduce (domain, &n);
  gfs_norm_update (&n);

  return n;
}

/**
 * gfs_domain_read:
 * @fp: a #GtsFile.
 *
 * Reads the graph nodes (#GfsBox) and edges and the
 * corresponding boundaries (#GfsBoundaryMpi if necessary) defined in
 * @fp.
 *
 * Returns: the #GfsDomain or %NULL if an error occured, in which case
 * the corresponding @fp fields (@pos and @error) are set.
 */
GfsDomain * gfs_domain_read (GtsFile * fp)
{
  GfsDomain * domain;

  g_return_val_if_fail (fp != NULL, NULL);
						 
  if (!(domain = GFS_DOMAIN (gts_graph_read (fp))))
    return NULL;

  (* GFS_DOMAIN_CLASS (GTS_OBJECT (domain)->klass)->post_read) (domain, fp);
  if (fp->type == GTS_ERROR) {
    gts_object_destroy (GTS_OBJECT (domain));
    return NULL;
  }

  return domain;
}

typedef struct {
  GSList * boxlist;
  guint bid;
  gboolean one_box_per_pe;
  gint pid;
  GfsVariable * newboxp;
  GfsDomain * domain;
} SplitPar;

static void box_split (GfsBox * box, SplitPar * p)
{
  guint refid = FTT_DIMENSION == 2 ? 2 : 6;
  FttCellChildren child;
  FttDirection d;
  guint i;
  GfsDomain * domain = gfs_box_domain (box);

  p->boxlist = g_slist_prepend (p->boxlist, box);

  if (FTT_CELL_IS_LEAF (box->root))
    ftt_cell_refine_single (box->root, (FttCellInitFunc) gfs_cell_init, domain);

  ftt_cell_children (box->root, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i]) {
      GfsBox * newbox = GFS_BOX (gts_object_new (GTS_OBJECT (box)->klass));

      GTS_OBJECT (newbox)->reserved = domain;
      if (p->one_box_per_pe)
	newbox->pid = (p->pid)++;
      else
	newbox->pid = box->pid;
      if (box->id == 1 && i == refid)
	newbox->id = 1;
      else
	newbox->id = (p->bid)++;

      GFS_DOUBLE_TO_POINTER (GFS_VALUE (child.c[i], p->newboxp)) = newbox;

      if (FTT_CELL_IS_LEAF (child.c[i]))
	ftt_cell_refine_single (child.c[i], (FttCellInitFunc) gfs_cell_init, domain);
    }

  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d])) {
      GfsBoundary * boundary = GFS_BOUNDARY (box->neighbor[d]);
      GfsBoundaryClass * klass = GFS_BOUNDARY_CLASS (GTS_OBJECT (boundary)->klass);

      ftt_cell_destroy (boundary->root, (FttCellCleanupFunc) gfs_cell_cleanup, domain);
      boundary->root = NULL;
      
      ftt_cell_children_direction (box->root, d, &child);
      for (i = 0; i < FTT_CELLS/2; i++)
	if (child.c[i]) {
	  GfsBox * newbox = GFS_DOUBLE_TO_POINTER (GFS_VALUE (child.c[i], p->newboxp));
	  GtsObject * newboundary = GTS_OBJECT (gfs_boundary_new (klass, newbox, d));

	  if (GFS_IS_BOUNDARY_PERIODIC (newboundary)) {
	    GFS_BOUNDARY_PERIODIC (newboundary)->matching = 
	      GFS_BOUNDARY_PERIODIC (boundary)->matching;
	    GFS_BOUNDARY_PERIODIC (newboundary)->d = GFS_BOUNDARY_PERIODIC (boundary)->d;
	  }
	  else
	    gfs_object_clone (GTS_OBJECT (boundary), GTS_OBJECT (newboundary));
	}
      gts_object_destroy (GTS_OBJECT (boundary));
      box->neighbor[d] = NULL;
    }
}

static GtsGEdge * node_is_linked (GtsGNode * n1, GtsGNode * n2, FttDirection d)
{
  GSList * i = GTS_SLIST_CONTAINER (n1)->items;
  while (i) {
    if (GTS_GNODE_NEIGHBOR (n1, i->data) == n2 &&
	GFS_GEDGE (i->data)->d == d)
      return i->data;
    i = i->next;
  }
  return NULL;
}

static void box_link (GfsBox * box, SplitPar * p)
{
  FttCellChildren child;
  guint i;

  ftt_cell_children (box->root, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i]) {
       GfsBox * newbox = GFS_DOUBLE_TO_POINTER (GFS_VALUE (child.c[i], p->newboxp));
       FttDirection d;
       
       g_assert (newbox);
       gts_container_add (GTS_CONTAINER (p->domain), GTS_CONTAINEE (newbox));

       for (d = 0; d < FTT_NEIGHBORS; d++)
	 if (newbox->neighbor[d] != NULL && GFS_IS_BOUNDARY_PERIODIC (newbox->neighbor[d])) {
	   GfsBox * matching =  GFS_BOUNDARY_PERIODIC (newbox->neighbor[d])->matching;
	   static gint match[FTT_CELLS][FTT_NEIGHBORS] = {
#if FTT_2D
	     { -1, 1, 2, -1 },
	     { 0, -1, 3, -1 },
	     { -1, 3, -1, 0 },
	     { 2, -1, -1, 1 }
#else /* 3D */
	     { -1, 1, 2, -1, 4, -1 },
	     { 0, -1, 3, -1, 5, -1 },
	     { -1, 3, -1, 0, 6, -1 },
	     { 2, -1, -1, 1, 7, -1 },
	     { -1, 5, 6, -1, -1, 0 },
	     { 4, -1, 7, -1, -1, 1 },
	     { -1, 7, -1, 4, -1, 2 },
	     { 6, -1, -1, 5, -1, 3 },	     
#endif /* 3D */
	   };
	   gint ci = match[FTT_CELL_ID (child.c[i])][d];
	   g_assert (ci >= 0);
	   FttCellChildren neighbors;
	   ftt_cell_children (matching->root, &neighbors);
	   FttCell * neighbor = neighbors.c[ci];
	   g_assert (neighbor);
	   GfsBox * newbox1 = GFS_DOUBLE_TO_POINTER (GFS_VALUE (neighbor, p->newboxp));
	   g_assert (newbox1);
	   GFS_BOUNDARY_PERIODIC (newbox->neighbor[d])->matching = newbox1;
	   if (!node_is_linked (GTS_GNODE (newbox1), GTS_GNODE (newbox), 
				FTT_OPPOSITE_DIRECTION (d))) {
	     GfsGEdge * edge = GFS_GEDGE (gts_gedge_new (GTS_GRAPH (p->domain)->edge_class,
							 GTS_GNODE (newbox), 
							 GTS_GNODE (newbox1)));
	     edge->d = d;
	   }
	 }
	 else if (newbox->neighbor[d] == NULL) {
	   FttCell * neighbor = ftt_cell_neighbor (child.c[i], d);

	   if (neighbor) {
	     GfsBox * newbox1 = GFS_DOUBLE_TO_POINTER (GFS_VALUE (neighbor, p->newboxp));
	     FttDirection od = FTT_OPPOSITE_DIRECTION (d);
	     GfsGEdge * edge;

	     g_assert (newbox1);
	     newbox->neighbor[d] = GTS_OBJECT (newbox1);
	     g_assert (newbox1->neighbor[od] == NULL);
	     newbox1->neighbor[od] = GTS_OBJECT (newbox);
	     edge = GFS_GEDGE (gts_gedge_new (GTS_GRAPH (p->domain)->edge_class,
					      GTS_GNODE (newbox), 
					      GTS_GNODE (newbox1)));
	     edge->d = d;
	   }
	 }
    }
}

static void box_destroy (GfsBox * box, GfsVariable * newboxp)
{
  GfsBox * newbox[FTT_CELLS];
  FttCellChildren child;
  guint i;

  ftt_cell_children (box->root, &child);
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i])
      newbox[i] = GFS_DOUBLE_TO_POINTER (GFS_VALUE (child.c[i], newboxp));
    else
      newbox[i] = NULL;

  ftt_cell_destroy_root (box->root, &child, (FttCellCleanupFunc) gfs_cell_cleanup, 
			 gfs_box_domain (box));
  box->root = NULL;
  for (i = 0; i < FTT_CELLS; i++)
    if (child.c[i]) {
      newbox[i]->root = child.c[i];
      FTT_ROOT_CELL (newbox[i]->root)->parent = newbox[i];
    }

  gts_object_destroy (GTS_OBJECT (box));
}

static void get_ref_pos (GfsBox * box, FttVector * pos)
{
  if (box->id == 1)
    ftt_cell_pos (box->root, pos);
}

/**
 * gfs_domain_split:
 * @domain: a #GfsDomain.
 * @one_box_per_pe: if %TRUE each new box created is assigned to a
 * different process, otherwise the newly created box inherits the pid
 * of its parent.
 *
 * Splits each box of @domain into its (4 in 2D, 8 in 3D)
 * children. The corresponding newly created boxes are added to the
 * graph and the parent boxes are destroyed.
 */
void gfs_domain_split (GfsDomain * domain, gboolean one_box_per_pe)
{
  SplitPar p;

  g_return_if_fail (domain != NULL);

  p.newboxp = gfs_temporary_variable (domain);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, 1,
  			   (FttCellTraverseFunc) gfs_cell_reset, p.newboxp);
  p.boxlist = NULL;
  p.bid = 2;
  p.pid = 0;
  p.one_box_per_pe = one_box_per_pe;
  p.domain = domain;
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_split, &p);
  g_slist_foreach (p.boxlist, (GFunc) box_link, &p);
  g_slist_foreach (p.boxlist, (GFunc) box_destroy, p.newboxp);
  g_slist_free (p.boxlist);
  gts_object_destroy (GTS_OBJECT (p.newboxp));

  gfs_domain_match (domain);
  domain->rootlevel++;
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) get_ref_pos, &domain->refpos);
}

/**
 * gfs_domain_locate:
 * @domain: a #GfsDomain.
 * @target: position of the point to look for.
 * @max_depth: maximum depth to consider (-1 means no restriction, see below for -2).
 * @where: a pointer to a #GfsBox or %NULL.
 *
 * Locates the cell of @domain containing @target. This is done
 * efficiently in log(n) operations by using the topology of the cell
 * trees.
 *
 * If @max_depth is set to -2, the finest cell containing @target is
 * returned. This cell is not necessarily a leaf-cell in contrast to
 * the case where @max_depth is set to -1.
 *
 * If @where is not %NULL it is filled with the #GfsBox containing the
 * cell.
 *
 * Returns: a #FttCell of @domain containing (boundary included) the
 * point defined by @target or %NULL if @target is not contained in
 * any cell of @domain.  
 */
FttCell * gfs_domain_locate (GfsDomain * domain,
			     FttVector target,
			     gint max_depth,
			     GfsBox ** where)
{
  g_return_val_if_fail (domain != NULL, NULL);
  g_return_val_if_fail (domain->array != NULL, NULL);

  GSList * b = gfs_locate_array_locate (domain->array, &target);
  if (b && GFS_IS_BOX (b->data)) {
    if (where)
      *where = b->data;
    return ftt_cell_locate (GFS_BOX (b->data)->root, target, max_depth);
  }
  return NULL;
}

/**
 * gfs_domain_boundary_locate:
 * @domain: a #GfsDomain.
 * @target: position of the point to look for.
 * @max_depth: maximum depth to consider (-1 means no restriction).
 * @where: a pointer to a #GtsObject.
 *
 * Locates the cell of @domain or of its boundary containing @target.
 *
 * If @where is not %NULL it is filled with the #GtsObject (either a
 * #GfsBox or a #GfsBoundary) containing the cell.
 *
 * Returns: a #FttCell of @domain or of its boundary containing the
 * point defined by @target or %NULL if @target is not contained in
 * any cell of @domain or of its boundary.
 */
FttCell * gfs_domain_boundary_locate (GfsDomain * domain,
				      FttVector target,
				      gint max_depth,
				      GtsObject ** where)
{
  g_return_val_if_fail (domain != NULL, NULL);
  g_return_val_if_fail (domain->array != NULL, NULL);

  GSList * b = gfs_locate_array_locate (domain->array, &target);
  if (!b)
    return NULL;
  if (GFS_IS_BOX (b->data)) {
    if (where)
      *where = b->data;
    return ftt_cell_locate (GFS_BOX (b->data)->root, target, max_depth);
  }
  else
    while (b) {
      g_assert (GFS_IS_BOUNDARY (b->data));
      FttCell * cell = ftt_cell_locate (GFS_BOUNDARY (b->data)->root, target, max_depth);
      if (cell && GFS_CELL_IS_BOUNDARY (cell)) {
	if (where)
	  *where = b->data;
	return cell;
      }
      b = b->next;
    }
  return NULL;
}

static void box_distance2 (GfsBox * box, GPtrArray * a)
{
  g_ptr_array_add (a, box);
}

static void bubble_sort (GPtrArray * a, gdouble * d)
{
  guint i, j;

  for (i = 0; i < a->len - 1; i++)
    for (j = 0; j < a->len - 1 - i; j++)
      if (d[j+1] < d[j]) {
	gdouble tmp = d[j];
	gpointer data = a->pdata[j];
	d[j] = d[j+1];
	d[j+1] = tmp;
	a->pdata[j] = a->pdata[j+1];
	a->pdata[j+1] = data;
      }
}

/**
 * gfs_domain_cell_point_distance2:
 * @domain: a #GfsDomain.
 * @p: a #GtsPoint.
 * @distance2: the squared distance function.
 * @data: user data to pass to @distance2.
 * @closest: where to return the closest cell or %NULL.
 *
 * For non-leafs cells @distance2 must return a lower-bound for the
 * minimum distance (using for example ftt_cell_point_distance2_min()).
 *
 * Returns: the square of the minimum distance measured according to
 * @distance2 between @p and a leaf cell of @domain.
 */
gdouble gfs_domain_cell_point_distance2 (GfsDomain * domain,
					 GtsPoint * p,
					 gdouble (* distance2) (FttCell *, GtsPoint *, gpointer),
					 gpointer data,
					 FttCell ** closest)
{
  gdouble dmin = G_MAXDOUBLE;
  GPtrArray * a;
  gdouble * d;
  guint i;

  g_return_val_if_fail (domain != NULL, dmin);
  g_return_val_if_fail (p != NULL, dmin);
  g_return_val_if_fail (distance2 != NULL, dmin);

  if (closest)
    *closest = NULL;
  a = g_ptr_array_new ();
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_distance2, a);
  d = g_malloc (sizeof (gdouble)*a->len);
  for (i = 0; i < a->len; i++)
    d[i] = (* distance2) (GFS_BOX (a->pdata[i])->root, p, data);
  bubble_sort (a, d);
  for (i = 0; i < a->len; i++)
    if (d[i] < dmin)
      ftt_cell_point_distance2_internal (GFS_BOX (a->pdata[i])->root, p, d[i],
					 distance2, data, closest, &dmin);
  g_free (d);
  g_ptr_array_free (a, TRUE);
  return dmin;
}

/**
 * gfs_domain_advect_point:
 * @domain: a #GfsDomain.
 * @p: a #FttVector.
 * @dt: the time step.
 *
 * Updates the coordinates of point @p at time t + @dt using the
 * velocity field defined by @domain.
 *
 * If @p is not contained within @domain, the coordinates are unchanged.
 */
void gfs_domain_advect_point (GfsDomain * domain, 
			      FttVector * p,
			      gdouble dt)
{
  FttCell * cell;
  FttVector p0, p1;
  FttComponent c;
  GfsVariable ** u;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (p != NULL);

  p0 = p1 = *p;
  cell = gfs_domain_locate (domain, p0, -1, NULL);
  if (cell == NULL)
    return;
  u = gfs_domain_velocity (domain);
  for (c = 0; c < FTT_DIMENSION; c++)
    (&p1.x)[c] += dt*gfs_interpolate (cell, p0, u[c])/2.;
  cell = gfs_domain_locate (domain, p1, -1, NULL);
  if (cell == NULL)
    return;
  for (c = 0; c < FTT_DIMENSION; c++)
    (&p->x)[c] += dt*gfs_interpolate (cell, p1, u[c]);
}

static void count (FttCell * cell, guint * n)
{
  (*n)++;
}

/**
 * gfs_domain_size:
 * @domain: a #GfsDomain.
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 *
 * Returns: the number of cells of @domain traversed using @flags and
 * @max_depth.
 */
guint gfs_domain_size (GfsDomain * domain,
		       FttTraverseFlags flags,
		       gint max_depth)
{
  guint n = 0;

  g_return_val_if_fail (domain != NULL, 0);
  
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth, 
			   (FttCellTraverseFunc) count, &n);
  gfs_all_reduce (domain, n, MPI_UNSIGNED, MPI_SUM);
  return n;
}

typedef struct {
  gdouble cfl;
  GfsVariable ** v;
  GfsDomain * domain;
} CflData;

static void minimum_mac_cfl (FttCellFace * face, CflData * p)
{
  gdouble un = GFS_STATE (face->cell)->f[face->d].un;
  gdouble length = ftt_cell_size (face->cell);
  if (p->domain->cell_metric) {
    gdouble fm = (* p->domain->face_metric) (p->domain, face);
    if (fm <= 0.) /* e.g. Axi metric on the axis */
      return;
    length *= (* p->domain->cell_metric) (p->domain, face->cell)/fm;
  }
  if (un != 0.) {
    gdouble cflu = length/fabs (un);
    if (cflu*cflu < p->cfl)
      p->cfl = cflu*cflu;
  }
  FttComponent c = face->d/2;
  if (p->v[c]->sources) {
    gdouble g = 0.;
    GSList * i = GTS_SLIST_CONTAINER (p->v[c]->sources)->items;
    while (i) {
      GfsSourceGeneric * s = i->data;
      if (s->face_value)
	g += (* s->face_value) (s, face, p->v[c]);
      i = i->next;
    }
    if (g != 0.) {
      gdouble cflg = 2.*length/fabs (g);
      if (cflg < p->cfl)
	p->cfl = cflg;
    }
  }
}

static void minimum_cfl (FttCell * cell, CflData * p)
{
  gdouble length = ftt_cell_size (cell);
  if (p->domain->cell_metric)
    length *= (* p->domain->cell_metric) (p->domain, cell);

  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++) {
    gdouble fm;
    if (p->domain->face_metric) {
      FttCellFace f;
      f.cell = cell; f.d = 2*c;
      gdouble fm1 = (* p->domain->face_metric) (p->domain, &f);
      f.d = 2*c + 1;
      gdouble fm2 = (* p->domain->face_metric) (p->domain, &f);
      fm = MAX (fm1, fm2);
    }
    else
      fm = 1.;
    if (GFS_VALUE (cell, p->v[c]) != 0.) {
      gdouble cflu = length/fabs (fm*GFS_VALUE (cell, p->v[c]));

      if (cflu*cflu < p->cfl)
	p->cfl = cflu*cflu;
    }
    if (p->v[c]->sources) {
      gdouble g = gfs_variable_mac_source (p->v[c], cell);

      if (g != 0.) {
	gdouble cflg = 2.*length/fabs (fm*g);

	if (cflg < p->cfl)
	  p->cfl = cflg;
      }
    }
  }
}

/**
 * gfs_domain_cfl:
 * @domain: a #GfsDomain.
 * @flags: which types of cells are to be visited.
 * @max_depth: maximum depth of the traversal.
 *
 * Returns: the minimum over the cells of @domain (traversed using
 * @flags and @max_depth) of the time scale defined by the size of the
 * cell and the norm of either the local velocity or the local
 * acceleration.
 */
gdouble gfs_domain_cfl (GfsDomain * domain,
			FttTraverseFlags flags,
			gint max_depth)
{
  CflData p;

  g_return_val_if_fail (domain != NULL, 0.);

  p.cfl = G_MAXDOUBLE;
  p.v = gfs_domain_velocity (domain);
  p.domain = domain;
  gfs_domain_face_traverse (domain, FTT_XYZ, FTT_PRE_ORDER, flags, max_depth, 
			    (FttFaceTraverseFunc) minimum_mac_cfl, &p);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, flags, max_depth, 
			    (FttCellTraverseFunc) minimum_cfl, &p);
  gfs_all_reduce (domain, p.cfl, MPI_DOUBLE, MPI_MIN);
  return sqrt (p.cfl);
}

/**
 * gfs_cell_init:
 * @cell: a #FttCell.
 * @domain: a #GfsDomain containing @cell.
 *
 * Allocates the memory for fluid state data associated to @cell or its children.
 */
void gfs_cell_init (FttCell * cell, GfsDomain * domain)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (domain != NULL);

  if (FTT_CELL_IS_LEAF (cell)) {
    g_return_if_fail (cell->data == NULL);
    cell->data = g_malloc0 (gfs_domain_variables_size (domain));
  }
  else {
    FttCellChildren child;
    guint n;

    ftt_cell_children (cell, &child);
    for (n = 0; n < FTT_CELLS; n++) {
      g_return_if_fail (child.c[n]->data == NULL);
      child.c[n]->data = g_malloc0 (gfs_domain_variables_size (domain));
    }
    if (GFS_CELL_IS_BOUNDARY (cell))
      for (n = 0; n < FTT_CELLS; n++)
	child.c[n]->flags |= GFS_FLAG_BOUNDARY;
  }
}

/**
 * gfs_cell_reinit:
 * @cell: a #FttCell.
 * @domain: a #GfsDomain containing @cell.
 *
 * Re-allocates the memory for fluid state data associated to @cell.
 */
void gfs_cell_reinit (FttCell * cell, GfsDomain * domain)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (cell->data != NULL);
  g_return_if_fail (domain != NULL);

  cell->data = g_realloc (cell->data, gfs_domain_variables_size (domain));
}

/**
 * gfs_cell_fine_init:
 * @parent: a #FttCell.
 * @domain: a #GfsDomain containing @parent.
 *
 * Initialises the children of @parent.
 */
void gfs_cell_fine_init (FttCell * parent, GfsDomain * domain)
{
  GSList * i;

  g_return_if_fail (parent != NULL);
  g_return_if_fail (!FTT_CELL_IS_LEAF (parent));
  g_return_if_fail (domain != NULL);

  gfs_cell_init (parent, domain);

  if (!GFS_CELL_IS_BOUNDARY (parent) && GFS_IS_MIXED (parent))
    gfs_solid_coarse_fine (parent, domain);

  /* metric is used by gfs_cell_coarse_fine(), make sure it is
     initialised first */
  i = domain->variables;
  while (i) {
    GfsVariable * v = i->data;
    if (GFS_IS_VARIABLE_METRIC (v))
      (* v->coarse_fine) (parent, v);
    i = i->next;
  }

  /* initialise remaining variables */
  i = domain->variables;
  while (i) {
    GfsVariable * v = i->data;
    if (!GFS_IS_VARIABLE_METRIC (v))
      (* v->coarse_fine) (parent, v);
    i = i->next;
  }
}

/**
 * gfs_cell_copy:
 * @from: a #FttCell to copy attributes from.
 * @to: a #FttCell to copy attributes to.
 * @domain: the #GfsDomain containing @from.
 *
 * Copies the attributes of the fluid cell @from to the fluid cell @to.
 */
void gfs_cell_copy (const FttCell * from, 
		    FttCell * to,
		    GfsDomain * domain)
{
  GfsSolidVector * solid;
  GfsStateVector * froms, * tos;

  g_return_if_fail (from != NULL);
  g_return_if_fail (to != NULL);
  g_return_if_fail (from != to);  
  g_return_if_fail (domain != NULL);

  froms = GFS_STATE (from);
  tos = GFS_STATE (to);
  if (froms != NULL) {
    if (tos == NULL) {
      gfs_cell_init (to, domain);
      tos = GFS_STATE (to);
    }
    solid = tos->solid;
    memcpy (to->data, from->data, gfs_domain_variables_size (domain));
    if (froms->solid == NULL) {
      if (solid)
	g_free (solid);
    }
    else {
      tos->solid = solid;
      *solid = *(froms->solid);
    }
  }
  else if (tos != NULL)
    gfs_cell_cleanup (to, domain);
}

/**
 * gfs_cell_write:
 * @cell: a #FttCell.
 * @fp: a file pointer.
 * @variables: the list of #GfsVariable to be written.
 *
 * Writes in @fp the fluid data associated with @cell and described by
 * @variables. This function is generally used in association with
 * ftt_cell_write().  
 */
void gfs_cell_write (const FttCell * cell, FILE * fp,
		     GSList * variables)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (fp != NULL);

  if (GFS_IS_MIXED (cell)) {
    GfsStateVector * s = GFS_STATE (cell);
    guint i;

    for (i = 0; i < FTT_NEIGHBORS; i++)
      fprintf (fp, " %g", s->solid->s[i]);
    fprintf (fp, " %g", s->solid->a);
    for (i = 0; i < FTT_DIMENSION; i++)
      fprintf (fp, " %g", (&s->solid->cm.x)[i]);
  }
  else
    fputs (" -1", fp);
  
  while (variables) {
    fprintf (fp, " %g", GFS_VALUE (cell, GFS_VARIABLE (variables->data)));
    variables = variables->next;
  }
}

/**
 * gfs_cell_read:
 * @cell: a #FttCell.
 * @fp: a #GtsFile.
 * @domain: the #GfsDomain containing @cell.
 *
 * Reads from @fp the fluid data associated with @cell and described
 * by @domain->variables_io. This function is generally used in
 * association with ftt_cell_read().  
 */
void gfs_cell_read (FttCell * cell, GtsFile * fp, GfsDomain * domain)
{
  gdouble s0;
  GfsStateVector * s;
  GSList * i;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (fp != NULL);
  g_return_if_fail (domain != NULL);

  if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
    gts_file_error (fp, "expecting a number (solid->s[0])");
    return;
  }
  s0 = atof (fp->token->str);
  if (s0 < 0. && s0 != -1.) {
    gts_file_error (fp, "solid->s[0] must be positive");
    return;
  }
  gts_file_next_token (fp);

  gfs_cell_init (cell, domain);
  s = cell->data;
  if (s0 >= 0.) {
    guint i;

    s->solid = g_malloc0 (sizeof (GfsSolidVector));
    s->solid->s[0] = s0;

    for (i = 1; i < FTT_NEIGHBORS; i++) {
      if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
	gts_file_error (fp, "expecting a number (solid->s[%d])", i);
	return;
      }
      s->solid->s[i] = atof (fp->token->str);
      gts_file_next_token (fp);
    }
    if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
      gts_file_error (fp, "expecting a number (solid->a)");
      return;
    }
    s->solid->a = atof (fp->token->str);
    gts_file_next_token (fp);
    for (i = 0; i < FTT_DIMENSION; i++) {
      if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
	gts_file_error (fp, "expecting a number (solid->cm[%d])", i);
	return;
      }
      (&s->solid->cm.x)[i] = atof (fp->token->str);
      gts_file_next_token (fp);
    }
  }

  i = domain->variables_io;
  while (i) {
    GfsVariable * v = i->data;

    if (fp->type != GTS_INT && fp->type != GTS_FLOAT) {
      gts_file_error (fp, "expecting a number (%s)", v->name);
      return;
    }
    GFS_VALUE (cell, v) = atof (fp->token->str);
    gts_file_next_token (fp);
    i = i->next;
  }
}

/**
 * gfs_cell_write_binary:
 * @cell: a #FttCell.
 * @fp: a file pointer.
 * @variables: the list of #GfsVariable to be written.
 *
 * Writes in @fp the fluid data associated with @cell and described by
 * @variables. This function is generally used in association with
 * ftt_cell_write_binary().
 */
void gfs_cell_write_binary (const FttCell * cell, FILE * fp,
			    GSList * variables)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (fp != NULL);

  if (GFS_IS_MIXED (cell)) {
    GfsStateVector * s = GFS_STATE (cell);

    fwrite (s->solid->s, sizeof (gdouble), FTT_NEIGHBORS, fp);
    fwrite (&s->solid->a, sizeof (gdouble), 1, fp);
    fwrite (&s->solid->cm.x, sizeof (gdouble), FTT_DIMENSION, fp);
    fwrite (&s->solid->ca.x, sizeof (gdouble), FTT_DIMENSION, fp);
  }
  else {
    gdouble a = -1.;
    fwrite (&a, sizeof (gdouble), 1, fp);
  }
  
  while (variables) {
    gdouble a = GFS_VALUE (cell, GFS_VARIABLE (variables->data));
    fwrite (&a, sizeof (gdouble), 1, fp);
    variables = variables->next;
  }
}

/**
 * gfs_cell_read_binary:
 * @cell: a #FttCell.
 * @fp: a #GtsFile.
 * @domain: the #GfsDomain containing @cell.
 *
 * Reads from @fp the fluid data associated with @cell and described
 * by @domain->variables_io. This function is generally used in
 * association with ftt_cell_read_binary().
 */
void gfs_cell_read_binary (FttCell * cell, GtsFile * fp, GfsDomain * domain)
{
  gdouble s0;
  GfsStateVector * s;
  GSList * i;

  g_return_if_fail (cell != NULL);
  g_return_if_fail (fp != NULL);
  g_return_if_fail (domain != NULL);

  if (gts_file_read (fp, &s0, sizeof (gdouble), 1) != 1) {
    gts_file_error (fp, "expecting a number (solid->s[0])");
    return;
  }
  if (s0 < 0. && s0 != -1.) {
    gts_file_error (fp, "solid->s[0] must be positive");
    return;
  }

  gfs_cell_init (cell, domain);
  s = cell->data;
  if (s0 >= 0.) {
    s->solid = g_malloc0 (sizeof (GfsSolidVector));
    s->solid->s[0] = s0;
    
    if (gts_file_read (fp, &s->solid->s[1], sizeof (gdouble), FTT_NEIGHBORS - 1) 
	!= FTT_NEIGHBORS - 1) {
      gts_file_error (fp, "expecting numbers (solid->s[1..%d])", FTT_NEIGHBORS - 1);
      return;
    }
    if (gts_file_read (fp, &s->solid->a, sizeof (gdouble), 1) != 1) {
      gts_file_error (fp, "expecting a number (solid->a)");
      return;
    }
    if (gts_file_read (fp, &s->solid->cm.x, sizeof (gdouble), FTT_DIMENSION) != FTT_DIMENSION) {
      gts_file_error (fp, "expecting numbers (solid->cm[0..%d])", FTT_DIMENSION - 1);
      return;
    }
    if (domain->version >= 90628 &&
	gts_file_read (fp, &s->solid->ca.x, sizeof (gdouble), FTT_DIMENSION) != FTT_DIMENSION) {
      gts_file_error (fp, "expecting numbers (solid->ca[0..%d])", FTT_DIMENSION - 1);
      return;
    }
  }

  i = domain->variables_io;
  while (i) {
    GfsVariable * v = i->data;
    gdouble a;

    if (gts_file_read (fp, &a, sizeof (gdouble), 1) != 1) {
      gts_file_error (fp, "expecting a number (%s)", v->name);
      return;
    }
    GFS_VALUE (cell, v) = a;
    i = i->next;
  }
}

static void box_realloc (GfsBox * box, GfsDomain * domain)
{
  FttDirection d;

  ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
		     (FttCellTraverseFunc) gfs_cell_reinit, domain);
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOUNDARY (box->neighbor[d]))
      ftt_cell_traverse (GFS_BOUNDARY (box->neighbor[d])->root, 
			 FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			 (FttCellTraverseFunc) gfs_cell_reinit, domain);
}

/**
 * gfs_domain_alloc:
 * @domain: a #GfsDomain.
 *
 * Returns: the index of a memory location newly allocated for each
 * cell of @domain.
 */
guint gfs_domain_alloc (GfsDomain * domain)
{
  guint i = 0;

  g_return_val_if_fail (domain != NULL, -1);

  while (i < domain->allocated->len && g_array_index (domain->allocated, gboolean, i))
    i++;
  if (i == domain->allocated->len) {
    g_array_set_size (domain->allocated, domain->allocated->len + 1);
    gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) box_realloc, domain);
  }
  g_array_index (domain->allocated, gboolean, i) = TRUE;
  return i;
}

/**
 * gfs_domain_free:
 * @domain: a #GfsDomain.
 * @i: a memory location index previously allocated using gfs_domain_alloc().
 *
 * Frees the memory location of @domain defined by @i.
 */
void gfs_domain_free (GfsDomain * domain, guint i)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (i < domain->allocated->len);
  g_return_if_fail (g_array_index (domain->allocated, gboolean, i));

  g_array_index (domain->allocated, gboolean, i) = FALSE;
}

/**
 * gfs_domain_add_variable:
 * @domain: a #GfsDomain.
 * @name: the name of the variable to add or %NULL.
 * @description: the variable description or %NULL.
 *
 * Adds a new variable @name to @domain.
 *
 * Returns: the new variable or %NULL if a variable with the same name
 * already exists.  
 */
GfsVariable * gfs_domain_add_variable (GfsDomain * domain,
				       const gchar * name,
				       const gchar * description)
{
  GfsVariable * v;

  g_return_val_if_fail (domain != NULL, NULL);

  if ((v = gfs_variable_new (gfs_variable_class (), domain, name, description)) == NULL)
    return NULL;
  domain->variables = g_slist_append (domain->variables, v);
  return v;
}

/**
 * gfs_domain_get_or_add_variable:
 * @domain: a #GfsDomain.
 * @name: the name of the variable to add or get.
 * @description: the variable description or %NULL.
 *
 * Adds a new variable @name to @domain or returns the variable of
 * @domain with the same name. In either case the description of the
 * variable name is set to @description (if it is not %NULL).
 *
 * Returns: the new or already existing variable or %NULL if @name is a
 * reserved variable name.
 */
GfsVariable * gfs_domain_get_or_add_variable (GfsDomain * domain,
					      const gchar * name,
					      const gchar * description)
{
  GfsVariable * v;

  g_return_val_if_fail (domain != NULL, NULL);
  g_return_val_if_fail (name != NULL, NULL);

  v = gfs_variable_from_name (domain->variables, name);
  if (v != NULL) {
    if (description) {
      if (v->description)
	g_free (v->description);
      v->description = g_strdup (description);
    }
  }
  else
    v = gfs_domain_add_variable (domain, name, description);
  return v;
}

typedef struct {
  gdouble * f, * m;
  GfsVariable * v;
  GfsFunction * weight;
  GfsSourceDiffusion * d;
} Force;

static void add_pressure_force (FttCell * cell, Force * f)
{
  gdouble weight = f->weight ? gfs_function_value (f->weight, cell) : 1.;

  if (weight != 0.) {
    gdouble * r = &GFS_STATE (cell)->solid->ca.x;
    FttVector ff, mm;
    FttComponent c;
    
    gfs_pressure_force (cell, f->v, &ff);
    gts_vector_cross (&mm.x, r, &ff.x);
    for (c = 0; c < 3; c++) {
      f->f[c] += weight*(&ff.x)[c];
      f->m[c] += weight*(&mm.x)[c];
    }
  }
}

static GfsSourceDiffusion * source_diffusion (GfsVariable * v)
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

static void add_viscous_force (FttCell * cell, Force * f)
{
  gdouble weight = f->weight ? gfs_function_value (f->weight, cell) : 1.;

  if (weight != 0.) {
    gdouble D;
    GfsSolidVector * s = GFS_STATE (cell)->solid;
    gdouble * r = &s->ca.x;
    FttVector ff, mm, n, g;
    FttComponent c;
    
    g_assert (((cell)->flags & GFS_FLAG_DIRICHLET) != 0);
    gfs_cell_dirichlet_gradient (cell, f->v->i, -1, s->fv, &g);
    
    D = - gfs_source_diffusion_cell (f->d, cell);
    n.x = s->s[1] - s->s[0];
    n.y = s->s[3] - s->s[2];
#if FTT_2D
    ff.z = 0.;
    switch (f->v->component) {
    case FTT_X:
      ff.x = D*(2.*g.x*n.x + g.y*n.y);
      ff.y = D*g.y*n.x;
      break;
    case FTT_Y:
      ff.x = D*g.x*n.y;
      ff.y = D*(2.*g.y*n.y + g.x*n.x);
      break;
    default:
      g_assert_not_reached ();
    }
#else /* 3D */
    n.z = s->s[5] - s->s[4];
    D *= ftt_cell_size (cell);
    switch (f->v->component) {
    case FTT_X:
      ff.x = D*(2.*g.x*n.x + g.y*n.y + g.z*n.z);
      ff.y = D*g.y*n.x;
      ff.z = D*g.z*n.x;
      break;
    case FTT_Y:
      ff.y = D*(2.*g.y*n.y + g.x*n.x + g.z*n.z);
      ff.x = D*g.x*n.y;
      ff.z = D*g.z*n.y;
      break;
    case FTT_Z:
      ff.z = D*(2.*g.z*n.z + g.x*n.x + g.y*n.y);
      ff.x = D*g.x*n.z;
      ff.y = D*g.y*n.z;
      break;
    default:
      g_assert_not_reached ();
    }
#endif /* 3D */
    gts_vector_cross (&mm.x, r, &ff.x);
    for (c = 0; c < 3; c++) {
      f->f[c] += weight*(&ff.x)[c];
      f->m[c] += weight*(&mm.x)[c];
    }
  }
}

/**
 * gfs_domain_solid_force:
 * @domain: a #GfsDomain.
 * @pf: a #FttVector.
 * @vf: a #FttVector.
 * @pm: a #FttVector.
 * @vm: a #FttVector.
 * @weight: an optional weight.
 *
 * Fills @pf and @vf (resp. @pm and @vm) with the components of the
 * net pressure and viscous forces (resp. pressure and viscous
 * moments) applied by the fluid on the solid surface embedded in
 * @domain.
 *
 * The reference point for the moments is the origin of the coordinate system.
 */
void gfs_domain_solid_force (GfsDomain * domain, 
			     FttVector * pf,
			     FttVector * vf,
			     FttVector * pm,
			     FttVector * vm,
			     GfsFunction * weight)
{
  FttComponent c;
  GfsVariable ** v;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (pf != NULL);
  g_return_if_fail (vf != NULL);
  g_return_if_fail (pm != NULL);
  g_return_if_fail (vm != NULL);

  if (GFS_IS_AXI (domain))
    g_assert_not_implemented ();

  pf->x = pf->y = pf->z = 0.;
  pm->x = pm->y = pm->z = 0.;
  Force f;
  f.f = (gdouble *) pf;
  f.m = (gdouble *) pm;
  f.v = gfs_variable_from_name (domain->variables, "P");
  f.weight = weight;
  if (weight)
    gfs_catch_floating_point_exceptions ();
  gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			     (FttCellTraverseFunc) add_pressure_force, &f);
  if (weight)
    gfs_restore_fpe_for_function (weight);
  vf->x = vf->y = vf->z = 0.;
  vm->x = vm->y = vm->z = 0.;
  v = gfs_domain_velocity (domain);
  for (c = 0; c < FTT_DIMENSION; c++) {
    GfsSourceDiffusion * D = source_diffusion (v[c]);

    if (D) {
      gfs_domain_surface_bc (domain, v[c]);
      f.f = (gdouble *) vf;
      f.m = (gdouble *) vm;
      f.v = v[c];
      f.d = D;
      f.weight = weight;
      if (weight)
	gfs_catch_floating_point_exceptions ();
      gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
				 (FttCellTraverseFunc) add_viscous_force, &f);
      if (weight)
	gfs_restore_fpe_for_function (weight);
    }
  }

  for (c = 0; c < 3; c++) {
    gfs_all_reduce (domain, (&pf->x)[c], MPI_DOUBLE, MPI_SUM);
    gfs_all_reduce (domain, (&vf->x)[c], MPI_DOUBLE, MPI_SUM);
    gfs_all_reduce (domain, (&pm->x)[c], MPI_DOUBLE, MPI_SUM);
    gfs_all_reduce (domain, (&vm->x)[c], MPI_DOUBLE, MPI_SUM);
  }
}

#define THRESHOLD 1e-4

static void tag_cell_fraction (GtsFifo * fifo,
			       FttCell * cell,
			       GfsVariable * c, GfsVariable * v,
			       guint tag)
{
  FttDirection d;
  FttCellNeighbors n;

  g_assert (FTT_CELL_IS_LEAF (cell));
  ftt_cell_neighbors (cell, &n);
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (n.c[d] && GFS_VALUE (n.c[d], v) == 0. && GFS_VALUE (n.c[d], c) > THRESHOLD) {
      if (FTT_CELL_IS_LEAF (n.c[d])) {
	GFS_VALUE (n.c[d], v) = tag;
	gts_fifo_push (fifo, n.c[d]);
      }
      else {
	FttCellChildren child;
	FttDirection od = FTT_OPPOSITE_DIRECTION (d);
	guint i;

	ftt_cell_children_direction (n.c[d], od, &child);
	for (i = 0; i < FTT_CELLS/2; i++)
	  if (child.c[i] && GFS_VALUE (child.c[i], v) == 0. &&
	      GFS_VALUE (child.c[i], c) > THRESHOLD) {
	    GFS_VALUE (child.c[i], v) = tag;
	    gts_fifo_push (fifo, child.c[i]);
	  }
      }
    }
}

typedef struct {
  GfsVariable * v, * c;
  FttDirection d;
  guint * touch, * tags, tag, tagshift;
} TagPar;

static void tag_new_fraction_region (FttCell * cell, TagPar * p)
{
  if (GFS_VALUE (cell, p->v) == 0. && GFS_VALUE (cell, p->c) > THRESHOLD) {
    GtsFifo * fifo = gts_fifo_new ();

    GFS_VALUE (cell, p->v) = ++p->tag;
    gts_fifo_push (fifo, cell);
    while ((cell = gts_fifo_pop (fifo)))
      tag_cell_fraction (fifo, cell, p->c, p->v, p->tag);
    gts_fifo_destroy (fifo);
  }
}

/* @touch defines the touching connectivity. This function updates
   @touch with the info that region tagged with @tag1 touches the
   region tagged with @tag2  */
static void touching_regions (guint tag1, guint tag2, guint * touch)
{
  if (tag2 < tag1) {
    guint tmp = tag1;
    tag1 = tag2;
    tag2 = tmp;
  }
  else if (tag2 == tag1)
    return;
  guint ntag = touch[tag2];
  if (ntag == tag1)
    return;
  if (ntag == 0)
    touch[tag2] = tag1;
  else {
    if (tag1 < ntag)
      touch[tag2] = tag1;
    touching_regions (tag1, ntag, touch);
  }
}

#ifdef HAVE_MPI
static void reduce_touching_regions (void * in, void * inout, int * len, MPI_Datatype * type)
{
  guint * ltouch = (guint *) in;
  guint * gtouch = (guint *) inout;
  guint i;

  for (i = 1; i < *len; i++)
    if (ltouch[i] > 0)
      touching_regions (i, ltouch[i], gtouch);
}

static void shift_tags (FttCell * cell, TagPar * p)
{
  if (GFS_VALUE (cell, p->v) > 0.)
    GFS_VALUE (cell, p->v) += p->tagshift;
}
#endif /* HAVE_MPI */

static void unify_tag_range (GfsDomain * domain, TagPar * p)
{
#ifdef HAVE_MPI
  if (domain->pid >= 0) {
    int gsize;
    guint * tags;
    MPI_Comm_size (MPI_COMM_WORLD, &gsize);
    tags = g_malloc (sizeof (guint)*gsize);
    MPI_Allgather (&p->tag, 1, MPI_UNSIGNED, tags, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
    /* tags[] now contains the p->tag value on each PE */
    guint i;
    p->tag = 0;
    for (i = 0; i < gsize; i++)
      p->tag += tags[i];
    /* shift tag values to get a single tag space across all PEs */
    if (domain->pid > 0) {
      p->tagshift = 0;
      for (i = 0; i < domain->pid; i++)
	p->tagshift += tags[i];
      gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				(FttCellTraverseFunc) shift_tags, p);
    }
    g_free (tags);
  }
#endif /* HAVE_MPI */
}

static void match_periodic_bc (FttCell * cell, TagPar * p)
{
  guint tag = GFS_VALUE (cell, p->v);
  if (tag > 0) {
    FttCell * neighbor = ftt_cell_neighbor (cell, p->d);
    guint ntag = GFS_VALUE (neighbor, p->v);
    if (ntag > 0)
      touching_regions (tag, ntag, p->touch);
  }
}

static void match_box_bc (GfsBox * box, TagPar * p)
{
  for (p->d = 0; p->d < FTT_NEIGHBORS; p->d++)
    if (GFS_IS_BOUNDARY_PERIODIC (box->neighbor[p->d]))
      ftt_cell_traverse_boundary (box->root, p->d, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
				  (FttCellTraverseFunc) match_periodic_bc, p);
}

static void fix_touching (FttCell * cell, TagPar * p)
{
  GFS_VALUE (cell, p->v) = p->tags[p->touch[(guint) GFS_VALUE (cell, p->v)]];
}

/**
 * gfs_domain_tag_droplets:
 * @domain: a #GfsDomain.
 * @c: the volume fraction.
 * @tag: a #GfsVariable.
 *
 * Fills the @tag variable of the cells of @domain with the (strictly
 * positive) index of the droplet they belong to. The cells belonging
 * to the background phase have an index of zero.
 *
 * Note that the volume fraction @c must be defined on all levels.
 *
 * Returns: the number of droplets.
 */
guint gfs_domain_tag_droplets (GfsDomain * domain,
			       GfsVariable * c,
			       GfsVariable * tag)
{
  g_return_val_if_fail (domain != NULL, 0);
  g_return_val_if_fail (c != NULL, 0);
  g_return_val_if_fail (tag != NULL, 0);

  TagPar p;
  gboolean touching = FALSE;
  p.c = c;
  p.v = tag;
  p.tag = 0;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, tag);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) tag_new_fraction_region, &p);

  /* the rest of the algorithm deals with periodic and parallel BCs */
  unify_tag_range (domain, &p);
  gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, tag);
  p.touch = g_malloc0 ((p.tag + 1)*sizeof (guint));
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) match_box_bc, &p);

#ifdef HAVE_MPI
  if (domain->pid >= 0) {
    guint * gtouch = g_malloc0 ((p.tag + 1)*sizeof (guint));
    MPI_Op op;    
    MPI_Op_create (reduce_touching_regions, FALSE, &op);
    MPI_Allreduce (p.touch, gtouch, p.tag + 1, MPI_UNSIGNED, op, MPI_COMM_WORLD);
    MPI_Op_free (&op);
    g_free (p.touch);
    p.touch = gtouch;
  }
#endif /* HAVE_MPI */
  
  /* find region with smallest tag touching each region i.e. reduces
     the chain of touching tags */
  guint i, maxtag = 0;
  for (i = 1; i <= p.tag; i++) {
    guint touch = p.touch[i];
    while (touch > 0) {
      p.touch[i] = touch;
      touch = p.touch[touch];
      touching = TRUE;
    }
    if (p.touch[i] == 0 && i > maxtag)
      maxtag = i;
  }
  
  /* fix touching regions */
  if (touching) {
    guint ntag = 0; /* fresh tag index */
    p.tags = g_malloc ((maxtag + 1)*sizeof (guint));
    p.tags[0] = 0;
    for (i = 1; i <= maxtag; i++)
      if (p.touch[i] == 0) { /* this region is not touching any other */
	p.touch[i] = i;
	p.tags[i] = ++ntag;
      }
    maxtag = ntag;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) fix_touching, &p);
    g_free (p.tags);
  }

  g_free (p.touch);
  return maxtag;
}

typedef struct {
  GfsVariable * tag, * c;
  guint * sizes;
  guint n, min;
  gdouble val;
} RemoveDropletsPar;

static void compute_droplet_size (FttCell * cell, RemoveDropletsPar * p)
{
  guint i = GFS_VALUE (cell, p->tag);
  if (i > 0)
    p->sizes[i - 1]++;
}

static void reset_small_fraction (FttCell * cell, RemoveDropletsPar * p)
{
  guint i = GFS_VALUE (cell, p->tag);
  if (i > 0 && p->sizes[i - 1] < p->min)
    GFS_VALUE (cell, p->c) = p->val;
}

static int greater (const void * a, const void * b)
{
  return *((guint *)a) > *((guint *)b) ? -1 : 1;
}

/**
 * gfs_domain_remove_droplets:
 * @domain: a #GfsDomain.
 * @c: a #GfsVariable.
 * @v: a #GfsVariable.
 * @min: the minimum size (in cells) of the droplets.
 * @val: the value used to reset @v.
 *
 * Resets the @v variable (using @val) of all the droplets (defined by
 * the @c variable) smaller than @min cells if @min is positive, or
 * all the droplets but the -$min largest ones if @min is negative.
 */
void gfs_domain_remove_droplets (GfsDomain * domain,
				 GfsVariable * c,
				 GfsVariable * v,
				 gint min,
				 gdouble val)
{
  RemoveDropletsPar p;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (c != NULL);
  g_return_if_fail (v != NULL);

  p.c = c;
  p.tag = gfs_temporary_variable (domain);
  p.n = gfs_domain_tag_droplets (domain, c, p.tag);
  if (p.n > 0 && -min < (gint) p.n) {
    p.sizes = g_malloc0 (p.n*sizeof (guint));
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) compute_droplet_size, &p);
#ifdef HAVE_MPI
    if (domain->pid >= 0) {
      guint * sizes = g_malloc0 (p.n*sizeof (guint));
      MPI_Allreduce (p.sizes, sizes, p.n, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
      g_free (p.sizes);
      p.sizes = sizes;
    }
#endif
    if (min >= 0)
      p.min = min;
    else {
      guint * tmp = g_malloc (p.n*sizeof (guint));
      memcpy (tmp, p.sizes, p.n*sizeof (guint));
      qsort (tmp, p.n, sizeof (guint), greater);
      g_assert (-1 - min < p.n);
      p.min = tmp[-1 - min];
      g_free (tmp);
    }
    p.c = v;
    p.val = val;
    gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			      (FttCellTraverseFunc) reset_small_fraction, &p);
    g_free (p.sizes);
  }
  gts_object_destroy (GTS_OBJECT (p.tag));
}

static void tag_cell (GtsFifo * fifo, FttCell * cell, GfsVariable * v, guint tag, guint * size)
{
  FttDirection d;
  FttCellNeighbors n;
  GfsSolidVector * solid = GFS_STATE (cell)->solid;

  g_assert (FTT_CELL_IS_LEAF (cell));
  (*size)++;
  ftt_cell_neighbors (cell, &n);
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (n.c[d] && GFS_VALUE (n.c[d], v) == 0. &&
	!GFS_CELL_IS_BOUNDARY (n.c[d]) &&
	(!solid || solid->s[d] > 0.)) {
      if (FTT_CELL_IS_LEAF (n.c[d])) {
	GFS_VALUE (n.c[d], v) = tag;
	gts_fifo_push (fifo, n.c[d]);
      }
      else {
	FttCellChildren child;
	FttDirection od = FTT_OPPOSITE_DIRECTION (d);
	guint i, j;
	
	j = ftt_cell_children_direction (n.c[d], od, &child);
	for (i = 0; i < j; i++)
	  if (child.c[i] && GFS_VALUE (child.c[i], v) == 0. &&
	      (!GFS_IS_MIXED (child.c[i]) || GFS_STATE (child.c[i])->solid->s[od] > 0.)) {
	    GFS_VALUE (child.c[i], v) = tag;
	    gts_fifo_push (fifo, child.c[i]);
	  }
      }
    }
}

static void tag_new_region (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];

  if (GFS_VALUE (cell, v) == 0.) {
    GArray * sizes = data[1];
    guint size = 0;
    GtsFifo * fifo = gts_fifo_new ();

    GFS_VALUE (cell, v) = sizes->len + 1;
    gts_fifo_push (fifo, cell);
    while ((cell = gts_fifo_pop (fifo)))
      tag_cell (fifo, cell, v, sizes->len + 1, &size);
    g_array_append_val (sizes, size);
    gts_fifo_destroy (fifo);
  }
}

static gboolean remove_small (FttCell * cell, gpointer * data)
{
  if (FTT_CELL_IS_LEAF (cell)) {
    GArray * sizes = data[0];
    GfsVariable * v = data[5];
    guint * min = data[1], i = GFS_VALUE (cell, v) - 1.;

    g_assert (GFS_VALUE (cell, v) > 0.);
    if (g_array_index (sizes, guint, i) < *min) {
      if (FTT_CELL_IS_ROOT (cell))
	g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR, "root cell belongs to a pond");
      else
	ftt_cell_destroy (cell, data[2], data[3]);
      return TRUE;
    }
    return FALSE;
  }
  else {
    FttCellChildren child;
    guint i;
    gboolean changed = FALSE;

    ftt_cell_children (cell, &child);
    for (i = 0; i < FTT_CELLS; i++)
      if (child.c[i] && remove_small (child.c[i], data))
	changed = TRUE;
    if (FTT_CELL_IS_LEAF (cell)) {
      /* all the children have been destroyed i.e. the cell belongs to a small pond */
      if (FTT_CELL_IS_ROOT (cell))
	g_log (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR, "root cell belongs to a pond");
      else
	ftt_cell_destroy (cell, data[2], data[3]);
    }
    else if (changed)
      gfs_cell_init_solid_fractions_from_children (cell);
    return changed;
  }
}

static void remove_small_box (GfsBox * box, gpointer * data)
{
  gboolean * changed = data[4];

  if (remove_small (box->root, data))
    *changed = TRUE;
}

/**
 * gfs_domain_remove_ponds:
 * @domain: a #GfsDomain.
 * @min: the minimum size (in cells) of the ponds.
 * @cleanup: a #FttCellCleanupFunc or %NULL.
 * @data: user data to pass to @cleanup.
 *
 * Removes all the fluid "ponds" of @domain smaller than @min cells
 * if @min is positive, or all the ponds but the -@min largest ones
 * if @min is negative.
 *
 * If the domain is modified its boundaries are re"matched" using
 * gfs_domain_match().
 */
void gfs_domain_remove_ponds (GfsDomain * domain, 
			      gint min,
			      FttCellCleanupFunc cleanup,
			      gpointer data)
{
  GArray * sizes;
  gpointer dat[6];
  guint minsize;
  gboolean changed = FALSE;
  GfsVariable * v;

  g_return_if_fail (domain != NULL);

  v = gfs_temporary_variable (domain);
  sizes = g_array_new (FALSE, FALSE, sizeof (guint));
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) gfs_cell_reset, v);
  dat[0] = v;
  dat[1] = sizes;
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) tag_new_region, dat);
  g_assert (sizes->len > 0);
  if (min >= 0)
    minsize = min;
  else if (-min >= sizes->len)
    minsize = 0;
  else {
    guint * tmp = g_malloc (sizes->len*sizeof (guint));
    memcpy (tmp, sizes->data, sizes->len*sizeof (guint));
    qsort (tmp, sizes->len, sizeof (guint), greater);
    minsize = tmp[-1 - min];
    g_free (tmp);
  }
  dat[0] = sizes;
  dat[1] = &minsize;
  dat[2] = cleanup;
  dat[3] = data;
  dat[4] = &changed;
  dat[5] = v;
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) remove_small_box, dat);
  g_array_free (sizes, TRUE);
  gts_object_destroy (GTS_OBJECT (v));
  if (changed)
    gfs_domain_match (domain);
}

static gboolean tag_speck (FttCell * cell, GfsVariable * v)
{
  if (GFS_VALUE (cell, v) == 0.) {
    FttDirection d;
    FttCellNeighbors n;
    GfsSolidVector * solid = GFS_STATE (cell)->solid;
    
    g_assert (FTT_CELL_IS_LEAF (cell));
    ftt_cell_neighbors (cell, &n);
    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (!n.c[d])
	return FALSE;
    GFS_VALUE (cell, v) = 1.;
    for (d = 0; d < FTT_NEIGHBORS; d++)
      if (GFS_VALUE (n.c[d], v) == 0. && 
	  !GFS_CELL_IS_BOUNDARY (n.c[d]) &&
	  solid->s[d] > 0. && solid->s[d] < 1.) {
	g_assert (GFS_IS_MIXED (n.c[d]));
	if (FTT_CELL_IS_LEAF (n.c[d])) {
	  if (!tag_speck (n.c[d], v)) {
	    GFS_VALUE (cell, v) = 0.;
	    return FALSE;
	  }
	}
	else {
	  FttCellChildren child;
	  FttDirection od = FTT_OPPOSITE_DIRECTION (d);
	  guint i;
	  
	  ftt_cell_children_direction (n.c[d], od, &child);
	  for (i = 0; i < FTT_CELLS/2; i++)
	    if (!child.c[i] || (GFS_VALUE (child.c[i], v) == 0 && 
				GFS_IS_MIXED (child.c[i]) &&
				!tag_speck (child.c[i], v))) {
	      GFS_VALUE (cell, v) = 0.;
	      return FALSE;
	    }
	}
      }
  }
  return TRUE;
}

static void fill_speck (FttCell * cell, gpointer * data)
{
  GfsVariable * v = data[0];

  if (GFS_VALUE (cell, v) == 1.) {
    gboolean * changed = data[1];
    g_free (GFS_STATE (cell)->solid);
    GFS_STATE (cell)->solid = NULL;
    *changed = TRUE;
  }
}

/**
 * gfs_domain_remove_specks:
 * @domain: a #GfsDomain.
 *
 * Removes all the solid "specks" of @domain. Solid specks are islands
 * which do not contain any empty cell.
 *
 * Note that the domain's boundaries are not "matched" automatically.
 */
void gfs_domain_remove_specks (GfsDomain * domain)
{
  gboolean changed = FALSE;
  GfsVariable * v;
  gpointer data[2];

  g_return_if_fail (domain != NULL);

  v = gfs_temporary_variable (domain);
  gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL,
			     (FttCellTraverseFunc) gfs_cell_reset, v);
  gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			     (FttCellTraverseFunc) tag_speck, v);
  data[0] = v;
  data[1] = &changed;
  gfs_domain_traverse_mixed (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS,
			     (FttCellTraverseFunc) fill_speck, data);
  gts_object_destroy (GTS_OBJECT (v));
  if (changed)
    gfs_domain_cell_traverse (domain, FTT_POST_ORDER, FTT_TRAVERSE_NON_LEAFS, -1,
			      (FttCellTraverseFunc) gfs_cell_init_solid_fractions_from_children, 
			      NULL);
    
}

/**
 * gfs_domain_timer_start:
 * @domain: a #GfsDomain.
 * @name: the name of the timer.
 *
 * Starts timer @name of @domain. If @name does not exist it is
 * created first.
 */
void gfs_domain_timer_start (GfsDomain * domain, const gchar * name)
{
  GfsTimer * t;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (name != NULL);

  t = g_hash_table_lookup (domain->timers, name);
  if (t == NULL) {
    t = g_malloc (sizeof (GfsTimer));
    gts_range_init (&t->r);
    g_hash_table_insert (domain->timers, g_strdup (name), t);
  }
  else
    g_return_if_fail (t->start < 0.);
  t->start = gfs_clock_elapsed (domain->timer);
  gfs_debug ("starting %s at %g", name, t->start);
}

/**
 * gfs_domain_timer_stop:
 * @domain: a #GfsDomain.
 * @name: the name of the timer.
 *
 * Stops timer @name of @domain. This function fails if @name is not a
 * timer of @domain.
 */
void gfs_domain_timer_stop (GfsDomain * domain, const gchar * name)
{
  GfsTimer * t;
  gdouble end;

  g_return_if_fail (domain != NULL);
  end = gfs_clock_elapsed (domain->timer);
  g_return_if_fail (name != NULL);

  t = g_hash_table_lookup (domain->timers, name);
  g_return_if_fail (t != NULL);
  g_return_if_fail (t->start >= 0.);

  gts_range_add_value (&t->r, end - t->start);
  gts_range_update (&t->r);
  gfs_debug ("stopping %s: elapsed: %g", name, end - t->start);
  t->start = -1.;
}

static void cell_combine_traverse (FttCell * cell,
				   FttCell * parent,
				   FttCellCombineTraverseFunc inside,
				   gpointer idata,
				   FttCellTraverseFunc outside,
				   gpointer odata)
{
  FttCell * locate;
  FttVector p;

  ftt_cell_pos (cell, &p);
  locate = ftt_cell_locate (parent, p, ftt_cell_level (cell));
  if (locate == NULL) {
    if (outside)
      ftt_cell_traverse (cell, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1, outside, odata);
  }
  else {
    if (FTT_CELL_IS_LEAF (cell))
      (* inside) (cell, locate, idata);
    else {
      FttCellChildren child;
      guint i;

      ftt_cell_children (cell, &child);
      for (i = 0; i < FTT_CELLS; i++)
	if (child.c[i])
	  cell_combine_traverse (child.c[i], locate, inside, idata, outside, odata);
    }
  }  
}

static void box_combine_traverse (GfsBox * box, gpointer * data)
{
  FttVector p;
  FttCell * locate;

  ftt_cell_pos (box->root, &p);
  locate = gfs_domain_locate (data[0], p, ftt_cell_level (box->root), NULL);
  if (locate == NULL) {
    if (data[3])
      ftt_cell_traverse (box->root, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1, data[3], data[4]);
  }
  else
    cell_combine_traverse (box->root, locate, data[1], data[2], data[3], data[4]);
}

/**
 * gfs_domain_combine_traverse:
 * @domain1: a #GfsDomain.
 * @domain2: another #GfsDomain.
 * @inside: function to call for each pair of cells.
 * @idata: user data to pass to @inside.
 * @outside: function to call for cells falling outside of @domain2 or
 * %NULL.
 * @odata: user data to pass to @outside.
 *
 * Calls @inside for each leaf cell of @domain1 contained in
 * @domain2. The second cell argument to @inside is set to the cell of
 * @domain2 containing the first cell argument.
 *
 * If @outside is not %NULL it is called for each leaf cell of
 * @domain1 which is outside of @domain2.
 */
void gfs_domain_combine_traverse (GfsDomain * domain1,
				  GfsDomain * domain2,
				  FttCellCombineTraverseFunc inside,
				  gpointer idata,
				  FttCellTraverseFunc outside,
				  gpointer odata)				  
{
  gpointer data[5];

  g_return_if_fail (domain1 != NULL);
  g_return_if_fail (domain2 != NULL);
  g_return_if_fail (inside != NULL);

  data[0] = domain2;
  data[1] = inside;
  data[2] = idata;
  data[3] = outside;
  data[4] = odata;

  gts_container_foreach (GTS_CONTAINER (domain1), (GtsFunc) box_combine_traverse, data);
}

/**
 * gfs_domain_add_derived_variable:
 * @domain: a #GfsDomain.
 * @info: the #GfsDerivedVariableInfo.
 *
 * Adds a derived variable described by @info to @domain.
 *
 * Returns: the #GfsDerivedVariable if the variable was successfully
 * added to @domain or %NULL if a variable with the same name already
 * exists.
 */
GfsDerivedVariable * gfs_domain_add_derived_variable (GfsDomain * domain, 
						      GfsDerivedVariableInfo info)
{
  GfsDerivedVariable * v;

  g_return_val_if_fail (domain != NULL, NULL);

  if (gfs_variable_from_name (domain->variables, info.name) ||
      gfs_derived_variable_from_name (domain->derived_variables, info.name))
    return NULL;
  v = GFS_DERIVED_VARIABLE (gts_object_new (GTS_OBJECT_CLASS (gfs_derived_variable_class ())));
  v->name = g_strdup (info.name);
  v->description = info.description ? g_strdup (info.description) : NULL;
  v->func = info.func;
  v->data = info.data;
  domain->derived_variables = g_slist_prepend (domain->derived_variables, v);
  GTS_OBJECT (v)->reserved = domain;
  return v;
}

/**
 * gfs_domain_remove_derived_variable:
 * @domain: a #GfsDomain.
 * @name: the name of a #GfsDerivedVariable.
 *
 * Removes derived variable @name from @domain.
 *
 * Returns: %TRUE if the variable was successfully removed from @domain or
 * %FALSE if a derived variable with the this name does not exist.
 */
gboolean gfs_domain_remove_derived_variable (GfsDomain * domain, const gchar * name)
{
  GSList * i;
  
  g_return_val_if_fail (domain != NULL, FALSE);
  g_return_val_if_fail (name != NULL, FALSE);

  i = domain->derived_variables;
  while (i) {
    GfsDerivedVariable * u = i->data;

    if (!strcmp (u->name, name)) {
      gts_object_destroy (GTS_OBJECT (u));
      domain->derived_variables = g_slist_remove_link (domain->derived_variables, i);
      g_slist_free (i);
      return TRUE;
    }
    i = i->next;
  }
  return FALSE;
}

typedef struct {
  FttDirection d;
  GfsFunction * f;
  GfsVariable * v;
} SumData;

static gdouble product (FttCell * cell, GfsFunction * f)
{
  GfsSolidVector * solid = GFS_STATE (cell)->solid;
  return ftt_cell_volume (cell)*(solid ? solid->a : 1.)*gfs_function_value (f, cell);
}

static void sum (FttCell * cell, SumData * data)
{
  FttCell * n = ftt_cell_neighbor (cell, data->d);
  GfsSolidVector * solid = GFS_STATE (cell)->solid;

  if (!n || GFS_CELL_IS_BOUNDARY (n) || (solid && solid->s[data->d] == 0.)) {
    gdouble s = 0.;

    n = cell;
    do {
      /* fixme: does not work if the resolution varies along data->d */
      g_assert (ftt_cell_level (n) == ftt_cell_level (cell));
      s += product (n, data->f);
      GFS_VALUE (n, data->v) = s;
      n = ftt_cell_neighbor (n, FTT_OPPOSITE_DIRECTION (data->d));
    } while (n && !GFS_CELL_IS_BOUNDARY (n) && 
	     (!GFS_IS_MIXED (n) || GFS_STATE (n)->solid->s[data->d] > 0.));
  }
}

/**
 * gfs_domain_sum:
 * @domain: a #GfsDomain.
 * @d: the #FttDirection.
 * @f: a #GfsFunction.
 * @v: a #GfsVariable.
 *
 * Fills variable @v of each cell of @domain with the sum in direction
 * @d of the volume-weighted function @f.
 */
void gfs_domain_sum (GfsDomain * domain, FttDirection d, GfsFunction * f, GfsVariable * v)
{
  SumData data;

  g_return_if_fail (domain != NULL);
  g_return_if_fail (d >= 0 && d < FTT_NEIGHBORS);
  g_return_if_fail (f != NULL);
  g_return_if_fail (v != NULL);

  data.d = d;
  data.f = f;
  data.v = v;
  gfs_catch_floating_point_exceptions ();
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) sum, &data);
  gfs_restore_fpe_for_function (f);
}

static void filter (FttCell * cell, gpointer * data)
{
  FttDirection d[4*(FTT_DIMENSION - 1)][FTT_DIMENSION] = {
#if FTT_2D
    {FTT_RIGHT, FTT_TOP}, {FTT_RIGHT, FTT_BOTTOM}, {FTT_LEFT, FTT_TOP}, {FTT_LEFT, FTT_BOTTOM}
#else
    {FTT_RIGHT, FTT_TOP, FTT_FRONT}, {FTT_RIGHT, FTT_BOTTOM, FTT_FRONT}, 
    {FTT_LEFT, FTT_TOP, FTT_FRONT}, {FTT_LEFT, FTT_BOTTOM, FTT_FRONT},
    {FTT_RIGHT, FTT_TOP, FTT_BACK}, {FTT_RIGHT, FTT_BOTTOM, FTT_BACK}, 
    {FTT_LEFT, FTT_TOP, FTT_BACK}, {FTT_LEFT, FTT_BOTTOM, FTT_BACK}
#endif
  };
  guint i;
  gdouble val = 0.;
  GfsVariable * a = data[0];
  GfsVariable * b = data[1];

  for (i = 0; i < 4*(FTT_DIMENSION - 1); i++)
    val += gfs_cell_corner_value (cell, d[i], a, -1);
  GFS_VALUE (cell, b) = val/(4*(FTT_DIMENSION - 1));
}

/**
 * gfs_domain_filter:
 * @domain: a #GfsDomain.
 * @v: a #GfsVariable.
 * @fv: the filtered variable or %NULL.
 *
 * Apply a "corner-averaging" filter to variable @v on all leaf cells
 * of @domain.
 *
 * If @fv is %NULL, @v is replaced by its filtered value.
 */
void gfs_domain_filter (GfsDomain * domain, GfsVariable * v, GfsVariable * fv)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (v != NULL);

  gpointer data[2];
  data[0] = v;
  data[1] = fv ? fv : gfs_temporary_variable (domain);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, -1,
			    (FttCellTraverseFunc) filter, data);
  if (fv == NULL) {
    gfs_variables_swap (data[0], data[1]);
    gts_object_destroy (data[1]);
    gfs_domain_bc (domain, FTT_TRAVERSE_LEAFS, -1, v);
  }
  else
    gfs_domain_copy_bc (domain, FTT_TRAVERSE_LEAFS, -1, v, fv);
}

struct _GfsRequest {
  void * buf;
#ifdef HAVE_MPI
  MPI_Request request[2];
#endif
};

/**
 * gfs_send_objects:
 * @list: a list of #GtsObject.
 * @dest: the rank of the destination PE.
 *
 * Sends the objects in @list to PE @dest of a parallel simulation.
 * This is a non-blocking operation which returns a handler which must
 * be cleared by calling gfs_wait().
 *
 * Note that this functions assumes that the write() method of the
 * #GtsObject sent begins by writing the object class name.
 *
 * Returns: a #GfsRequest.
 */
GfsRequest * gfs_send_objects (GSList * list, int dest)
{
#ifdef HAVE_MPI
  char * buf;
  size_t len;
  FILE * fp = open_memstream (&buf, &len);
  if (fp == NULL)
    g_error ("gfs_send_objects(): could not open_memstream:\n%s", strerror (errno));
  while (list) {
    GtsObject * object = list->data;
    g_assert (object->klass->write != NULL);
    (* object->klass->write) (object, fp);
    fputc ('\n', fp);
    list = list->next;
  }
  fclose (fp);
  GfsRequest * r = g_malloc0 (sizeof (GfsRequest));
  long length = len;
  MPI_Isend (&length, 1, MPI_LONG, dest, 0, MPI_COMM_WORLD, &r->request[0]);
  gfs_debug ("sending %ld bytes to PE %d", length, dest);
  if (length > 0) {
    r->buf = buf;
    MPI_Isend (r->buf, length, MPI_BYTE, dest, 1, MPI_COMM_WORLD, &r->request[1]);
  }
  return r;
#else  /* not HAVE_MPI */
  return NULL;
#endif /* HAVE_MPI */
}

/**
 * gfs_wait:
 * @r: a #GfsRequest.
 *
 * Waits for completion of and deallocates @r.
 */
void gfs_wait (GfsRequest * r)
{
#ifdef HAVE_MPI
  g_return_if_fail (r != NULL);

  MPI_Status status;
  MPI_Wait (&r->request[0], &status);
  if (r->buf) {
    MPI_Wait (&r->request[1], &status);
    free (r->buf);
  }
  g_free (r);
#endif /* HAVE_MPI */
}

/**
 * gfs_receive_objects:
 * @domain: a #GfsDomain.
 * @src: the rank of the source PE.
 *
 * Receives a list of #GtsObject from PE @src of a parallel simulation.
 *
 * Returns: a list of newly-allocated objects.
 */
GSList * gfs_receive_objects (GfsDomain * domain, int src)
{
  g_return_val_if_fail (domain != NULL, NULL);

#ifdef HAVE_MPI
  MPI_Status status;
  long length;
  MPI_Recv (&length, 1, MPI_LONG, src, 0, MPI_COMM_WORLD, &status);
  gfs_debug ("receiving %ld bytes from PE %d", length, src);
  if (length > 0) {
    char * buf = g_malloc (length);
    MPI_Recv (buf, length, MPI_BYTE, src, 1, MPI_COMM_WORLD, &status);
    GtsFile * fp = gts_file_new_from_buffer (buf, length);
    GSList * list = NULL;
    while (fp->type == GTS_STRING) {
      GtsObjectClass * klass = gfs_object_class_from_name (fp->token->str);
      if (klass == NULL)
	g_error ("gfs_receive_object():%d:%d: unknown class '%s'", 
		 fp->line, fp->pos, fp->token->str);
      GtsObject * object = gts_object_new (klass);
      gfs_object_simulation_set (object, domain);
      g_assert (klass->read);
      (* klass->read) (&object, fp);
      if (fp->type == GTS_ERROR)
	g_error ("gfs_receive_object():%d:%d: %s", fp->line, fp->pos, fp->error);
      list = g_slist_prepend (list, object);
      while (fp->type == '\n')
	gts_file_next_token (fp);
    }
    gts_file_destroy (fp);
    g_free (buf);
    return list;
  }
#endif /* HAVE_MPI */
  return NULL;
}

static void unlink_box (GfsBox * box, gint * dest)
{
  FttDirection d;
  for (d = 0; d < FTT_NEIGHBORS; d++)
    if (GFS_IS_BOX (box->neighbor[d])) {
      GfsBox * nbox = GFS_BOX (box->neighbor[d]);
      FttDirection od = FTT_OPPOSITE_DIRECTION (d);
      nbox->neighbor[od] = NULL;
      gfs_boundary_mpi_new (gfs_boundary_mpi_class (), nbox, od, *dest, box->id);
      box->neighbor[d] = NULL;
      gfs_boundary_mpi_new (gfs_boundary_mpi_class (), box, d, nbox->pid, nbox->id);
    }
    else if (GFS_IS_BOUNDARY_PERIODIC (box->neighbor[d]) && 
	     !GFS_IS_BOUNDARY_MPI (box->neighbor[d])) {
      GfsBoundaryPeriodic * boundary = GFS_BOUNDARY_PERIODIC (box->neighbor[d]);
      g_assert (boundary->matching);
      GfsBoundaryPeriodic * matching = 
	GFS_BOUNDARY_PERIODIC (boundary->matching->neighbor[boundary->d]);
      g_assert (GFS_IS_BOUNDARY_PERIODIC (matching));
      GfsBox * nbox = GFS_BOUNDARY (matching)->box;
      FttDirection od = FTT_OPPOSITE_DIRECTION (d);
      g_assert (nbox->neighbor[od] == GTS_OBJECT (matching));
      gts_object_destroy (GTS_OBJECT (matching));
      nbox->neighbor[od] = NULL;
      gfs_boundary_mpi_new (gfs_boundary_mpi_class (), nbox, od, *dest, box->id);
      gts_object_destroy (GTS_OBJECT (box->neighbor[d]));
      box->neighbor[d] = NULL;
      gfs_boundary_mpi_new (gfs_boundary_mpi_class (), box, d, nbox->pid, nbox->id);
    }
}

static void setup_binary_IO (GfsDomain * domain)
{
  /* make sure that all the variables are sent */
  g_slist_free (domain->variables_io);
  domain->variables_io = NULL;
  GSList * i = domain->variables;
  while (i) {
    if (GFS_VARIABLE (i->data)->name)
      domain->variables_io = g_slist_append (domain->variables_io, i->data);
    i = i->next;
  }
  domain->binary = TRUE;	
}

/**
 * gfs_send_boxes:
 * @domain: a #GfsDomain.
 * @boxes: a list of #GfsBox belonging to @domain.
 * @dest: the destination processor id.
 *
 * Send boxes to @dest and removes them from @domain.
 * This is a non-blocking operation.
 *
 * Returns: a #GfsRequest which must be cleared using gfs_wait().
 */
GfsRequest * gfs_send_boxes (GfsDomain * domain, GSList * boxes, int dest)
{
  g_return_val_if_fail (domain != NULL, NULL);
  g_return_val_if_fail (dest != domain->pid, NULL);

  g_slist_foreach (boxes, (GFunc) unlink_box, &dest);
  setup_binary_IO (domain);
  GfsRequest * r = gfs_send_objects (boxes, dest);
  g_slist_foreach (boxes, (GFunc) gts_object_destroy, NULL);
  gfs_locate_array_destroy (domain->array);
  domain->array = gfs_locate_array_new (domain);
  return r;
}

/**
 * gfs_receive_boxes:
 * @domain: a #GfsDomain.
 * @src: the source processor id.
 *
 * Receive boxes from @src and adds them to @domain.
 *
 * Returns: the list of boxes received.
 */
GSList * gfs_receive_boxes (GfsDomain * domain, int src)
{
  g_return_val_if_fail (domain != NULL, NULL);
  g_return_val_if_fail (src != domain->pid, NULL);

  setup_binary_IO (domain);
  GSList * boxes = gfs_receive_objects (domain, src);
  if (boxes) {
    /* Create array for fast linking of ids to GfsBox pointers */
    GPtrArray * ids = box_ids (domain);
	  
    /* Convert internal GfsBoundaryMpi into graph edges */
    g_slist_foreach (boxes, (GFunc) convert_boundary_mpi_into_edges, ids);
    g_ptr_array_free (ids, TRUE);

    /* Update GfsLocateArray */
    gfs_locate_array_destroy (domain->array);
    domain->array = gfs_locate_array_new (domain);
  }
  return boxes;
}

/**
 * gfs_object_from_name:
 * @domain: a #GfsDomain.
 * @name: the name.
 *
 * Returns: the object of @domain called @name or %NULL.
 */
GtsObject * gfs_object_from_name (GfsDomain * domain, const gchar * name)
{
  g_return_val_if_fail (domain != NULL, NULL);
  g_return_val_if_fail (name != NULL, NULL);

  return g_hash_table_lookup (domain->objects, name);
}

/** \endobject{GfsDomain} */

/**
 * Projection of a GfsDomain along a coordinate direction.
 * \beginobject{GfsDomainProjection}
 */

static void gfs_domain_projection_destroy (GtsObject * o)
{
  GfsDomainProjection * p = GFS_DOMAIN_PROJECTION (o);
  p->domain->projections = g_slist_remove (p->domain->projections, p);

  (* GTS_OBJECT_CLASS (gfs_domain_projection_class ())->parent_class->destroy) (o);
} 

static void gfs_domain_projection_class_init (GfsEventClass * klass)
{
  GTS_OBJECT_CLASS (klass)->destroy = gfs_domain_projection_destroy;
}

GfsDomainClass * gfs_domain_projection_class (void)
{
  static GfsDomainClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo info = {
      "GfsDomainProjection",
      sizeof (GfsDomainProjection),
      sizeof (GfsDomainClass),
      (GtsObjectClassInitFunc) gfs_domain_projection_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_domain_class ()), &info);
  }

  return klass;
}

typedef struct {
  GfsDomainProjection * proj;
  FttCell * cell;
  GfsVariable * maxlevel;
} ProjData;

static gboolean overlap (FttCell * cell1, gpointer data)
{
  ProjData * r = data;
  FttCell * cell2 = r->cell;
  if (ftt_cell_level (cell1) < ftt_cell_level (cell2)) {
    FttCell * tmp = cell2;
    cell2 = cell1;
    cell1 = tmp;
  }
  FttVector p, q;
  ftt_cell_pos (cell1, &p);
  ftt_cell_pos (cell2, &q);
  gdouble h = ftt_cell_size (cell2)/2.;
  FttComponent c;
  for (c = 0; c < FTT_DIMENSION; c++)
    if (c != r->proj->c && ((&p.x)[c] < (&q.x)[c] - h || (&p.x)[c] > (&q.x)[c] + h))
      return FALSE;
  return TRUE;
}

static void update_maxlevel (FttCell * cell, int * maxlevel)
{
  int level = ftt_cell_level (cell);
  if (level > *maxlevel)
    *maxlevel = level;
}

static void project_refine (FttCell * cell, ProjData * p)
{
  int level = ftt_cell_level (cell), maxlevel = 0;
  p->cell = cell;
  gfs_domain_cell_traverse_condition (p->proj->domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, level + 1,
				      (FttCellTraverseFunc) update_maxlevel, &maxlevel,
				      overlap, p);
  GFS_VALUE (cell, p->maxlevel) = maxlevel;
  if (FTT_CELL_IS_LEAF (cell)) {
    if (maxlevel > level) {
      ftt_cell_refine_single (cell, 
			      GFS_DOMAIN (p->proj)->cell_init, 
			      GFS_DOMAIN (p->proj)->cell_init_data);
      ftt_cell_flatten (cell, 2*p->proj->c, (FttCellCleanupFunc) gfs_cell_cleanup, p->proj);
    }
  }
  else
    ftt_cell_flatten (cell, 2*p->proj->c, (FttCellCleanupFunc) gfs_cell_cleanup, p->proj);
}

static gboolean finer (FttCell * cell, ProjData * p)
{
  int level = ftt_cell_level (cell);
  g_assert (level >= GFS_VALUE (cell, p->maxlevel));
  return (level > GFS_VALUE (cell, p->maxlevel));
}

static void project_coarsen_box (GfsBox * box, ProjData * p)
{
  ftt_cell_coarsen (box->root,
		    (FttCellCoarsenFunc) finer, p,
		    (FttCellCleanupFunc) gfs_cell_cleanup, p->proj);
}

/**
 * gfs_domain_projection_reshape:
 * @proj: a #GfsDomainProjection.
 *
 * Updates the mesh for projection @proj.
 */
void gfs_domain_projection_reshape (GfsDomainProjection * proj)
{
  g_return_if_fail (proj != NULL);

  ProjData p = { proj };
  GfsDomain * domain = GFS_DOMAIN (proj);
  p.maxlevel = gfs_temporary_variable (domain);
  gfs_domain_cell_traverse (domain, FTT_PRE_ORDER, FTT_TRAVERSE_ALL, -1,
			    (FttCellTraverseFunc) project_refine, &p);
  gts_container_foreach (GTS_CONTAINER (domain), (GtsFunc) project_coarsen_box, &p);
  gts_object_destroy (GTS_OBJECT (p.maxlevel));
  gfs_domain_reshape (domain, gfs_domain_depth (domain));
}

/**
 * gfs_domain_projection_new:
 * @domain: a #GfsDomain.
 * @c: the component aligned with the projection direction.
 *
 * Returns: a new #GfsDomainProjection, projection of @domain along @c.
 */
GfsDomainProjection * gfs_domain_projection_new (GfsDomain * domain,
						 FttComponent c)
{
  g_return_val_if_fail (domain != NULL, NULL);
  g_return_val_if_fail (c < FTT_DIMENSION, NULL);

  /* clone domain */
  char * buf;
  size_t len;
  FILE * f = open_memstream (&buf, &len);
  if (f == NULL)
    g_error ("gfs_domain_projection_new(): could not open_memstream:\n%s", strerror (errno));
  gint depth = domain->max_depth_write;
  domain->max_depth_write = -2; /* no variables, no cells */
  GtsObjectClass * klass = GTS_OBJECT (domain)->klass;
  GTS_OBJECT (domain)->klass = GTS_OBJECT_CLASS (gfs_domain_projection_class ());
  gts_graph_write (GTS_GRAPH (domain), f);
  GTS_OBJECT (domain)->klass = klass;
  domain->max_depth_write = depth;
  fclose (f);

  GtsFile * fp = gts_file_new_from_buffer (buf, len);
  GfsDomainProjection * proj = GFS_DOMAIN_PROJECTION (gfs_domain_read (fp));
  if (fp->type == GTS_ERROR)
    g_error ("gfs_domain_projection_new:\n%d:%d:%s", fp->line, fp->pos, fp->error);
  gts_file_destroy (fp);
  free (buf);
  gfs_clock_start (GFS_DOMAIN (proj)->timer);

  /* project domain */
  proj->c = c;
  proj->domain = domain;
  domain->projections = g_slist_prepend (domain->projections, proj);
  gfs_domain_projection_reshape (proj);

  return proj;
}

typedef struct {
  GfsDomainProjection * proj;
  FttCell * cell;
  GfsProjectionTraverseFunc func;
  gpointer data;
} ProjectionTraverse;

static void apply_func (FttCell * cell, ProjectionTraverse * p)
{
  (* p->func) (p->cell, cell, p->data);
}

static void traverse_overlapping (FttCell * cell, ProjectionTraverse * p)
{
  p->cell = cell;
  gfs_domain_cell_traverse_condition (p->proj->domain, FTT_PRE_ORDER, FTT_TRAVERSE_LEAFS, 
				      ftt_cell_level (cell),
				      (FttCellTraverseFunc) apply_func, p,
				      overlap, p);
}

/**
 * gfs_domain_projection_traverse:
 * @domain: a #GfsDomainProjection.
 * @order: the order in which the cells are visited - %FTT_PRE_ORDER,
 * %FTT_POST_ORDER. 
 * @flags: which types of children are to be visited.
 * @max_depth: the maximum depth of the traversal. Cells below this
 * depth will not be traversed. If @max_depth is -1 all cells in the
 * tree are visited.
 * @func: the function to call for each visited #FttCell.
 * @data: user data to pass to @func.
 *
 * For each cell of @domain defined by the traversal flags, traverses
 * all the overlapping leaf cells of its parent domain.
 */
void gfs_domain_projection_traverse (GfsDomainProjection * domain,
				     FttTraverseType order,
				     FttTraverseFlags flags,
				     gint max_depth,
				     GfsProjectionTraverseFunc func,
				     gpointer data)
{
  g_return_if_fail (domain != NULL);
  g_return_if_fail (func != NULL);

  ProjectionTraverse p = { domain, NULL, func, data };
  gfs_domain_cell_traverse (GFS_DOMAIN (domain), order, flags, max_depth,
			    (FttCellTraverseFunc) traverse_overlapping, &p);
}

/** \endobject{GfsDomainProjection} */
