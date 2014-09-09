static void traverse_face (FttCell * cell, gpointer * datum)
{
  FttDirection * d = datum[0];
  gint max_depth = *((gint *) datum[1]);
  FttFaceTraverseFunc func = (FttFaceTraverseFunc) datum[2];
  gpointer data = datum[3];
  gboolean check = *((gboolean *) datum[4]);
  gboolean boundary_faces = *((gboolean *) datum[5]);  
  FttCellFace face;
  
  face.d = *d;
  face.cell = cell;
  face.neighbor = ftt_cell_neighbor (cell, face.d);
  if (face.neighbor) {
    if (!check || (face.neighbor->flags & FTT_FLAG_TRAVERSED) == 0) {
      if (FTT_CELL_IS_LEAF (cell) && 
	  !FTT_CELL_IS_LEAF (face.neighbor) && 
	  (max_depth < 0 || ftt_cell_level (face.neighbor) < max_depth)) {
	/* coarse -> fine */
	FttCellChildren children;
	guint i, n;
	
	face.d = FTT_OPPOSITE_DIRECTION (face.d);
	n = ftt_cell_children_direction (face.neighbor, face.d, &children);
	face.neighbor = face.cell;
	for (i = 0; i < n; i++)
	  if ((face.cell = children.c[i]) && 
	      (!check || (face.cell->flags & FTT_FLAG_TRAVERSED) == 0))
	    (* func) (&face, data);
      }
      else
	(* func) (&face, data);
    }
  }
  else if (boundary_faces)
    (* func) (&face, data);
}

static void traverse_all_faces (FttCell * cell, gpointer * datum)
{
  FttDirection d;

  datum[0] = &d;
  for (d = 0; d < FTT_NEIGHBORS; d++)
    traverse_face (cell, datum);
  cell->flags |= FTT_FLAG_TRAVERSED;
}

static void traverse_all_direct_faces (FttCell * cell, gpointer * datum)
{
  FttDirection d;

  datum[0] = &d;
  for (d = 0; d < FTT_NEIGHBORS; d += 2)
    traverse_face (cell, datum);
  cell->flags |= FTT_FLAG_TRAVERSED;
}

static void traverse_face_direction (FttCell * cell, gpointer * datum)
{
  traverse_face (cell, datum);
  cell->flags |= FTT_FLAG_TRAVERSED;
}

static void traverse_face_component (FttCell * cell, gpointer * datum)
{
  FttComponent * c = datum[0];
  FttDirection d;

  datum[0] = &d;
  d = 2*(*c);
  traverse_face (cell, datum);
  d++;
  traverse_face (cell, datum);
  cell->flags |= FTT_FLAG_TRAVERSED;
  datum[0] = c;
}

static void reset_flag (FttCell * cell)
{
  cell->flags &= ~FTT_FLAG_TRAVERSED;
}
