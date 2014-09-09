/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2012 National Institute of Water and Atmospheric Research
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
 * \brief GfsFunction and various utility functions.
 */

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/times.h>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#include <math.h>
#include "config.h"
#include "solid.h"
#include "simulation.h"
#include "cartesian.h"

/*
 * get_tmp_file based on the mkstemp implementation from the GNU C library.
 * Copyright (C) 1991,92,93,94,95,96,97,98,99 Free Software Foundation, Inc.
 */
typedef gint (*GTmpFileCallback) (gchar *, gint, gint);

static gint get_tmp_file (gchar            *tmpl,
			  GTmpFileCallback  f,
			  int               flags,
			  int               mode)
{
  char *XXXXXX;
  int count, fd;
  static const char letters[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
  static const int NLETTERS = sizeof (letters) - 1;
  glong value;
  GTimeVal tv;
  static int counter = 0;

  g_return_val_if_fail (tmpl != NULL, -1);

  /* find the last occurrence of "XXXXXX" */
  XXXXXX = g_strrstr (tmpl, "XXXXXX");

  if (!XXXXXX || strncmp (XXXXXX, "XXXXXX", 6))
    {
      errno = EINVAL;
      return -1;
    }

  /* Get some more or less random data.  */
  g_get_current_time (&tv);
  value = (tv.tv_usec ^ tv.tv_sec) + getpid () + counter++;

  for (count = 0; count < 100; value += 7777, ++count)
    {
      glong v = value;

      /* Fill in the random bits.  */
      XXXXXX[0] = letters[v % NLETTERS];
      v /= NLETTERS;
      XXXXXX[1] = letters[v % NLETTERS];
      v /= NLETTERS;
      XXXXXX[2] = letters[v % NLETTERS];
      v /= NLETTERS;
      XXXXXX[3] = letters[v % NLETTERS];
      v /= NLETTERS;
      XXXXXX[4] = letters[v % NLETTERS];
      v /= NLETTERS;
      XXXXXX[5] = letters[v % NLETTERS];

      fd = f (tmpl, flags, mode);

      if (fd >= 0)
        return fd;
      else if (errno != EEXIST)
        /* Any other error will apply also to other names we might
         *  try, and there are 2^32 or so of them, so give up now.
         */
        return -1;
    }

  /* We got out of the loop because we ran out of combinations to try.  */
  errno = EEXIST;
  return -1;
}

#if !HAVE_G_MKDTEMP
/* we could use mkdtemp() but we don't trust some implementations
   (e.g. on AIX) */
static gint wrap_mkdir (gchar *tmpl,
			int    flags G_GNUC_UNUSED,
			int    mode)
{
  return mkdir (tmpl, mode);
}

gchar * g_mkdtemp (gchar * tmpl)
{
  if (get_tmp_file (tmpl, wrap_mkdir, 0, 0700) == -1)
    return NULL;
  else
    return tmpl;
}
#endif /* !HAVE_G_MKDTEMP */

static gint wrap_mkfifo (gchar *tmpl,
			 int    flags G_GNUC_UNUSED,
			 int    mode)
{
  return mkfifo (tmpl, mode);
}

/**
 * gfs_mkftemp:
 * @tmpl: template FIFO name
 *
 * Creates a temporary FIFO. See the mkfifo() documentation
 * on most UNIX-like systems.
 *
 * Returns: A pointer to @tmpl, which has been modified
 *     to hold the directory name.  In case of errors, %NULL is
 *     returned and %errno will be set.
 */
gchar * gfs_mkftemp (gchar * tmpl)
{
  if (get_tmp_file (tmpl, wrap_mkfifo, 0, 0600) == -1)
    return NULL;
  else
    return tmpl;
}

/**
 * gfs_template:
 *
 * Returns: the template for a temporary file name.
 */
gchar * gfs_template (void)
{
  gchar * tmpdir = getenv ("TMPDIR");
  return tmpdir ? g_strconcat (tmpdir, "/gfsXXXXXX", NULL) : g_strdup ("/tmp/gfsXXXXXX");
}

/**
 * gfs_object_clone:
 * @object: the #GtsObject to clone.
 * @clone: the clone of @object.
 *
 * Makes @clone the clone of @object using the write() and read()
 * methods of @object.
 */
void gfs_object_clone (GtsObject * object, GtsObject * clone)
{
  g_return_if_fail (object != NULL);
  g_return_if_fail (clone != NULL);
  g_return_if_fail (gts_object_class_is_from_class (clone->klass, object->klass));
  
  char * buf;
  size_t len;
  FILE * fp = open_memstream (&buf, &len);
  if (fp == NULL)
    g_error ("open_memstream: %s", strerror (errno));
  (* object->klass->write) (object, fp);
  fclose (fp);
  GtsFile * gfp = gts_file_new_from_buffer (buf, len);
  (* object->klass->read) (&clone, gfp);
  g_assert (gfp->type != GTS_ERROR);
  gts_file_destroy (gfp);
  free (buf);
}

/**
 * @c: a character.
 * @s: a string.
 *
 * Returns: %TRUE if @c belongs to @s, %FALSE otherwise.
 */
gboolean gfs_char_in_string (char c, const char * s)
{
  if (s == NULL)
    return FALSE;
  while (*s != '\0')
    if (*(s++) == c)
      return TRUE;
  return FALSE;
}

/**
 * gfs_file_statement:
 * @fp: a #GtsFile.
 *
 * Reads the next brackets-delimited ({...}) statemement in @fp,
 * including all comments.
 *
 * Returns: a newly allocated string containing the text of the next
 * statement in @fp, or %NULL if an error occured in which case
 * @fp->error is set.
 */
gchar * gfs_file_statement (GtsFile * fp)
{
  g_return_val_if_fail (fp != NULL, NULL);

  if (fp->type != '{') {
    gts_file_error (fp, "expecting an opening brace");
    return NULL;
  }
  GString * s = g_string_new ("");
  gchar empty[] = "", * comments = fp->comments;
  fp->comments = empty;
  guint scope = fp->scope_max;
  gint c = gts_file_getc (fp);
  while (c != EOF && fp->scope > scope) {
    g_string_append_c (s, c);
    c = gts_file_getc (fp);
  }
  fp->comments = comments;
  if (fp->scope != scope) {
    gts_file_error (fp, "parse error");
    g_string_free (s, TRUE);
    return NULL;
  }
  gchar * ret = s->str;
  g_string_free (s, FALSE);
  return ret;
}

typedef gdouble (* GfsFunctionFunc) (const FttCell * cell,
				     const FttCellFace * face,
				     GfsSimulation * sim,
				     GfsVariable ** var,
				     GfsDerivedVariable ** dvar);
typedef gdouble (* GfsFunctionDerivedFunc) (const FttCell * cell,
					    const FttCellFace * face,
					    GfsSimulation * sim,
					    gpointer data);

/**
 * Global functions.
 * \beginobject{GfsGlobal}
 */

struct _GfsGlobal {
  /*< private >*/
  GtsObject parent;

  /*< public >*/
  gchar * s;
  guint line;
  gboolean appended;
};

static void global_destroy (GtsObject * object)
{
  g_free (GFS_GLOBAL (object)->s);
  (* gfs_global_class ()->parent_class->destroy) (object);
}

static void global_write (GtsObject * object, FILE * fp)
{
  fprintf (fp, "%s {", object->klass->info.name);
  fputs (GFS_GLOBAL (object)->s, fp);
  fputs ("}\n", fp);
}

static void global_read (GtsObject ** object, GtsFile * fp)
{
  GfsGlobal * global = GFS_GLOBAL (*object);
  GtsObjectClass * klass;

  if (fp->type != GTS_STRING) {
    gts_file_error (fp, "expecting a string (GfsGlobalClass)");
    return;
  }
  klass = gfs_object_class_from_name (fp->token->str);
  if (klass == NULL) {
    gts_file_error (fp, "unknown class `%s'", fp->token->str);
    return;
  }
  if (!gts_object_class_is_from_class (klass, gfs_global_class ())) {
    gts_file_error (fp, "`%s' is not a GfsGlobal", fp->token->str);
    return;
  }
  gts_file_next_token (fp);
  global->line = fp->line;
  g_free (global->s);
  if ((global->s = gfs_file_statement (fp)))
    gts_file_next_token (fp);
}

static void gfs_global_class_init (GtsObjectClass * klass)
{
  klass->destroy = global_destroy;
  klass->read =    global_read;
  klass->write =   global_write;
}

GtsObjectClass * gfs_global_class (void)
{
  static GtsObjectClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_global_info = {
      "GfsGlobal",
      sizeof (GfsGlobal),
      sizeof (GtsObjectClass),
      (GtsObjectClassInitFunc) gfs_global_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
				  &gfs_global_info);
  }

  return klass;
}

static void gfs_global_append (GfsGlobal * g, GString * s)
{
  g_string_append_printf (s, "#line %d \"GfsGlobal\"\n", g->line);
  g_string_append (s, g->s);
  g->appended = TRUE;
  g_string_append_c (s, '\n');
}

/* GfsModule: Header */

typedef struct _GfsModule GfsModule;

struct _GfsModule {
  gchar * key;
  guint id;
  GModule * module;
  GSList * l; /* list of GfsFunctions */
};

/* GfsFunction: Header */

struct _GfsFunction {
  GtsObject parent;
  GString * expr;
  gboolean isexpr;
  GfsModule * module;
  GfsFunctionFunc f;
  gchar * sname;
  GtsSurface * s;
  GfsCartesianGrid * g;
  guint index[4];
  GfsVariable * v;
  GfsDerivedVariable * dv;
  gdouble val;
  gboolean spatial, constant, nomap;
  gdouble units;
  GfsVariable ** var;
  GfsDerivedVariable ** dvar;
};

/** \endobject{GfsGlobal} */

/* GfsModule: object */

static gchar * function_key (const GfsFunction * f)
{
  GString * s = g_string_new (f->expr->str);
  if (f->spatial)
    g_string_append (s, "spatial");
  else if (f->constant)
    g_string_append (s, "constant");
  GfsSimulation * sim = gfs_object_simulation (f);
  GSList * i = sim->globals;
  while (i) {
    g_string_append (s, GFS_GLOBAL (i->data)->s);
    i = i->next;
  }
  gchar * key = s->str;
  g_string_free (s, FALSE);
  return key;
}

static gchar * find_identifier (const gchar * s, const gchar * i)
{
  gchar * f = strstr (s, i);
  static gchar allowed[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_1234567890";

  while (f) {
    if (gfs_char_in_string (f[strlen(i)], allowed) ||
	(f > s && gfs_char_in_string (f[-1], allowed)))
      f = strstr (++f, i);
    else
      return f;
  }
  return NULL;
}

/* source code for functions pending compilation */
static GString * pending_functions = NULL;
static guint n_pending_functions = 0;

/* append source code for @f to pending compilations */
static void append_pending_function (const GfsFunction * f, guint line, guint id)
{
  GfsSimulation * sim = gfs_object_simulation (f);
  GfsDomain * domain = GFS_DOMAIN (sim);
  GSList * lv = NULL, * ldv = NULL, * i;

  if (!pending_functions) {
    pending_functions = g_string_new ("#include <stdlib.h>\n"
				      "#include <stdio.h>\n"
				      "#include <math.h>\n"
				      "#include <gfs.h>\n"
				      "#include <gerris/spatial.h>\n"
				      "#include <gerris/function.h>\n"
				      "typedef double (* Func) (const FttCell * cell,\n"
				      "                         const FttCellFace * face,\n"
				      "                         GfsSimulation * sim,\n"
				      "                         gpointer data);\n");
    g_slist_foreach (sim->globals, (GFunc) gfs_global_append, pending_functions);
  }
  else {
    i = sim->globals;
    while (i) {
      if (!GFS_GLOBAL (i->data)->appended)
	gfs_global_append (i->data, pending_functions);
      i = i->next;
    }
  }
    
  if (f->spatial)
    g_string_append_printf (pending_functions, 
			    "\ndouble f%u (double x, double y, double z, double t) {\n"
			    "  _x = x; _y = y; _z = z;\n", id);
  else if (f->constant)
    g_string_append_printf (pending_functions, "\ndouble f%u (void) {\n", id);
  else {
    g_string_append_printf (pending_functions, "char * variables%u[] = {", id);
    i = domain->variables;
    while (i) {
      if (GFS_VARIABLE (i->data)->name && 
	  find_identifier (f->expr->str, GFS_VARIABLE (i->data)->name)) {
	lv = g_slist_prepend (lv, i->data);
	g_string_append_printf (pending_functions, "\"%s\", ", GFS_VARIABLE (i->data)->name);
      }
      i = i->next;
    }
    g_string_append (pending_functions, "NULL};\n");
    lv = g_slist_reverse (lv);

    g_string_append_printf (pending_functions, "char * dvariables%u[] = {", id);
    i = domain->derived_variables;
    while (i) {
      GfsDerivedVariable * v = i->data;
      if (find_identifier (f->expr->str, v->name)) {
	ldv = g_slist_prepend (ldv, v);
	g_string_append_printf (pending_functions, "\"%s\", ", v->name);
      }
      i = i->next;
    }
    g_string_append (pending_functions, "NULL};\n");
    ldv = g_slist_reverse (ldv);

    g_string_append_printf (pending_functions,
			    "\ndouble f%u (FttCell * cell, FttCellFace * face,\n"
			    "            GfsSimulation * sim, GfsVariable ** var,\n"
			    "            GfsDerivedVariable ** dvar) {\n"
			    "  _sim = sim; _cell = cell;\n",
			    id);
    if (lv || ldv) {
      GSList * i = lv;

      while (i) {
	GfsVariable * v = i->data;
	g_string_append_printf (pending_functions, "  double %s;\n", v->name);
	i = i->next;
      }
      i = ldv;
      while (i) {
	GfsDerivedVariable * v = i->data;
	g_string_append_printf (pending_functions, "  double %s;\n", v->name);
	i = i->next;
      }
      if (lv) {
	int index = 0;
	g_string_append (pending_functions, "  if (cell) {\n");
	i = lv;
	while (i) {
	  GfsVariable * v = i->data;
	  g_string_append_printf (pending_functions,
		   "    %s = gfs_dimensional_value (var[%d], GFS_VALUE (cell, var[%d]));\n", 
				  v->name, index, index);
	  i = i->next; index++;
	}
	g_string_append (pending_functions, "  } else {\n");
	i = lv; index = 0;
	while (i) {
	  GfsVariable * v = i->data;
	  g_string_append_printf (pending_functions,
		   "    %s = gfs_dimensional_value (var[%d],\n"
		   "           gfs_face_interpolated_value_generic (face, var[%d]));\n", 
		   v->name, index, index);
	  i = i->next; index++;
	}
	g_string_append (pending_functions, "  }\n");
	g_slist_free (lv);
      }
      if (ldv) {
	i = ldv; int index = 0;
	while (i) {
	  GfsDerivedVariable * v = i->data;
	  g_string_append_printf (pending_functions,
		   "  %s = (* (Func) dvar[%d]->func) (cell, face, sim, dvar[%d]->data);\n", 
		   v->name, index, index);
	  i = i->next; index++;
	}
	g_slist_free (ldv);
      }
    }
  }
  g_string_append_printf (pending_functions, "#line %d \"GfsFunction\"\n", line);

  if (f->isexpr)
    g_string_append_printf (pending_functions, "return %s;\n}\n", f->expr->str);
  else {
    gchar * s = f->expr->str;
    guint len = strlen (s);
    g_assert (s[0] == '{' && s[len-1] == '}');
    s[len-1] = '\0';
    g_string_append_printf (pending_functions, "%s\n}\n", &s[1]);
    s[len-1] = '}';
  }
}

static GHashTable * get_function_cache (void)
{
  static GHashTable * function_cache = NULL;
  if (!function_cache)
    function_cache = g_hash_table_new (g_str_hash, g_str_equal);
  return function_cache;
}

static void link_module (GfsFunction * f)
{
  g_assert (f->module);
  GModule * module = f->module->module;
  g_assert (module);
  guint id = f->module->id;
  gchar * name = g_strdup_printf ("f%u", id);
  g_assert (g_module_symbol (module, name, (gpointer) &f->f));
  g_free (name);
  if (f->constant) {
    f->val = (* f->f) (NULL, NULL, NULL, NULL, NULL);
    f->f = NULL;
    if (f->expr) g_string_free (f->expr, TRUE);
    f->expr = NULL;
  }
  else if (!f->spatial) {
    char ** variables;
    name = g_strdup_printf ("variables%u", id);
    g_assert (g_module_symbol (module, name, (gpointer) &variables));
    g_free (name);
    char ** s = variables;
    int n = 0;
    while (*s) { n++; s++; }
    if (n > 0) {
      f->var = g_malloc (n*sizeof (GfsVariable *));
      GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (f));
      s = variables; n = 0;
      while (*s) {
	g_assert ((f->var[n] = gfs_variable_from_name (domain->variables, *s)));
	n++; s++;
      }
    }
    name = g_strdup_printf ("dvariables%u", id);
    g_assert (g_module_symbol (module, name, (gpointer) &variables));
    g_free (name);
    s = variables;
    n = 0;
    while (*s) { n++; s++; }
    if (n > 0) {
      f->dvar = g_malloc (n*sizeof (GfsDerivedVariable *));
      GfsDomain * domain = GFS_DOMAIN (gfs_object_simulation (f));
      s = variables; n = 0;
      while (*s) {
	g_assert ((f->dvar[n] = 
		   gfs_derived_variable_from_name (domain->derived_variables, *s)));
	n++; s++;
      }
    }
  }
}

static void gfs_module_ref (GfsModule * m, GfsFunction * f)
{
  f->module = m;
  m->l = g_slist_prepend (m->l, f);
  if (m->module)
    link_module (f);
}

static void gfs_module_unref (GfsModule * m, GfsFunction * f)
{
  m->l = g_slist_remove (m->l, f);
  /* modules are kept "forever" in case they come up again e.g. during
     dynamic load-balancing */
}

static gboolean lookup_function (GfsFunction * f)
{
  GHashTable * function_cache = get_function_cache ();
  gchar * key = function_key (f);
  GfsModule * module = g_hash_table_lookup (function_cache, key);
  g_free (key);
  if (module) {
    gfs_module_ref (module, f);
    return TRUE;
  }
  return FALSE;
}

static void gfs_module_new (GfsFunction * f, guint line)
{
  if (!lookup_function (f)) {
    GfsModule * m = g_malloc0 (sizeof (GfsModule));
    m->key = function_key (f);
    m->id = n_pending_functions++;
    g_hash_table_insert (get_function_cache (), m->key, m);
    append_pending_function (f, line, m->id);
    gfs_module_ref (m, f);
  }
}

static double current_time (void)
{
  GTimeVal r;
  g_get_current_time (&r);
  return r.tv_sec + 1e-6*r.tv_usec;
}

static GModule * compile (GtsFile * fp, const gchar * dirname, const gchar * finname)
{
  GModule * module = NULL;
  gfs_debug ("starting compilation");
  double start = current_time ();
  char pwd[512];
  g_assert (getcwd (pwd, 512));
  GString * build_command = g_string_new ("");
  g_string_printf (build_command,
		   "cd %s && %s/build_function "
#if FTT_2D
		   "gerris2D"
#else /* 3D */
		   "gerris3D"
#endif
		   " \"%s\""
		   , dirname, GFS_DATA_DIR, pwd);
  g_string_append (build_command, " > log 2>&1");
  gint status = system (build_command->str);
  g_string_free (build_command, TRUE);
  if (WIFSIGNALED (status) && (WTERMSIG (status) == SIGINT || WTERMSIG (status) == SIGQUIT))
    status = SIGQUIT;
  else if (status == -1 || WEXITSTATUS (status) != 0) {
    gchar * errname = g_strconcat (dirname, "/log", NULL);
    FILE * ferr = fopen (errname, "r");
    g_free (errname);
    gchar * needle;
    GString * msg = g_string_new ("");
    gint c;

    while ((c = fgetc (ferr)) != EOF)
      g_string_append_c (msg, c);
    fclose (ferr);
    while ((needle = strstr (msg->str, "GfsFunction:")))
      g_string_erase (msg, needle - msg->str, strlen ("GfsFunction:"));
    gts_file_error (fp, "error compiling expression\n%s", msg->str);
    g_string_free (msg, TRUE);
    status = SIGABRT;
  }
  else {
    gchar * mname = g_strconcat (dirname, "/module.so", NULL);
    gchar * path = g_module_build_path (GFS_MODULES_DIR, mname);
    module = g_module_open (path, 0);
    g_free (path);
    if (module == NULL)
      module = g_module_open (mname, 0);
    if (module == NULL)
      gts_file_error (fp, "cannot load module: %s", g_module_error ());
    g_free (mname);
  }
#if 1
  gchar * cleanup = g_strconcat ("rm -r -f ", dirname, NULL);
  gint status1 = system (cleanup);
  g_free (cleanup);
  if (status1 == -1 || WEXITSTATUS (status1) != 0)
    g_warning ("error when cleaning up %s", dirname);
#else
  g_warning ("not cleaning up %s", dirname);
#endif
  gfs_debug ("compilation completed in %g s", current_time () - start);
  return module;
}

static void update_module (gchar * key, GfsModule * m, GModule * module)
{
  if (m->module == NULL) {
    m->module = module;
    g_slist_foreach (m->l, (GFunc) link_module, module);
  }
}

/**
 * gfs_pending_functions_compilation:
 * @fp: a #GtsFile.
 *
 * Compiles and links pending #GfsFunction definitions.
 *
 * Compilation errors are reported in @fp.
 */
void gfs_pending_functions_compilation (GtsFile * fp)
{
  g_return_if_fail (fp != NULL);

  if (pending_functions && fp->type != GTS_ERROR) {
    gchar * dirname = gfs_template ();
    if (g_mkdtemp (dirname) == NULL) {
      gts_file_error (fp, "cannot create temporary directory\n%s", strerror (errno));
      g_free (dirname);
      return;
    }
    gchar * finname = g_strdup_printf ("%s/function.c", dirname);
    FILE * fin = fopen (finname, "w");
    fputs (pending_functions->str, fin);
    fclose (fin);
    GModule * module = compile (fp, dirname, finname);
    if (module)
      g_hash_table_foreach (get_function_cache (), (GHFunc) update_module, module);
    /* note that if there is an error in some pending functions
       (i.e. fp->type == GTS_ERROR) something needs to be done to fix
       this (e.g. abort the whole thing). */
    g_string_free (pending_functions, TRUE);
    pending_functions = NULL;
    n_pending_functions = 0;
    g_free (dirname);
    g_free (finname);
  }
}

/**
 * Numerical constants and expressions.
 * \beginobject{GfsFunction}
 */

static GtsSurface * read_surface (gchar * name, GtsFile * fp)
{
  FILE * fptr = fopen (name, "r");
  GtsSurface * s;
  GtsFile * fp1;

  if (fptr == NULL) {
    gts_file_error (fp, "cannot open file `%s'", name);
    return NULL;
  }
  fp1 = gts_file_new (fptr);
  s = gts_surface_new (gts_surface_class (), gts_face_class (), 
		       gts_edge_class (), gts_vertex_class ());
  if (gts_surface_read (s, fp1)) {
    gts_file_error (fp, "%s:%d:%d: %s", name, fp1->line, fp1->pos, fp1->error);
    gts_object_destroy (GTS_OBJECT (s));
    s = NULL;
  }
  gts_file_destroy (fp1);
  fclose (fptr);
  return s;
}

#define INDEX_T 6

static gboolean fit_index_dimension (GfsCartesianGrid * grid, guint * val, GtsFile * fp)
{
  guint i, j;
  gchar * list[INDEX_T + 2] = {"x", "y", "z", "rx", "ry", "rz", "t", NULL};

  if (grid->N > 4) {
    gts_file_error (fp, "Cartesian grids can only use four dimensions or less");
    return FALSE;
  }

  for (i = 0; i < grid->N; i++) {
    gchar ** name = list; j = 0;
    while (*name != NULL && strcmp (*name, grid->name[i])) {
      name++; j++;
    }
    if (!name) {
      gts_file_error (fp, "unknown Cartesian grid index `%s'", grid->name[i]);
      return FALSE;
    }
    val[i] = j;
  }
  return TRUE;
}

static gboolean cgd_is_spatial (GfsFunction * f)
{
  guint i;
  for (i = 0; i < f->g->N; i++)
    if (f->index[i] < INDEX_T)
      return TRUE;
  return FALSE;
}

static gdouble interpolated_cgd (const GfsFunction * f, const FttVector * r)
{
  gdouble vector[4];
  gdouble val;
  guint i;

  FttVector p = *r;
  gboolean mapped = FALSE;
  for (i = 0; i < f->g->N; i++)
    if (f->index[i] < 3) { /* x,y,z */
      if (!mapped) {
	gfs_simulation_map_inverse (gfs_object_simulation (f), &p);
	mapped = TRUE;
      }
      vector[i] = (&p.x)[f->index[i]];
    }
    else if (f->index[i] < 6) /* rx,ry,rz */
      vector[i] = (&r->x)[f->index[i] - 3];
    else { /* t */
      g_assert (f->index[i] == INDEX_T);
      vector[i] = gfs_object_simulation (f)->time.t;
    }

  if (!gfs_cartesian_grid_interpolate (f->g, vector, &val))
    return 0.;
  return val;
}

/**
 * gfs_function_expression:
 * @fp: a #GtsFile.
 * @is_expression: a pointer to a boolean or %NULL.
 *
 * Reads the expression (in which case @is_expression is set to %TRUE)
 * or function from @fp.
 *
 * Returns: a newly allocated GString containing the result or %NULL
 * in case of error.
 */
GString * gfs_function_expression (GtsFile * fp, gboolean * is_expression)
{
  GString * expr = NULL;

  g_return_val_if_fail (fp != NULL, NULL);

  if (is_expression)
    *is_expression = TRUE;
  if (fp->type == '{') {
    gchar * s = gfs_file_statement (fp);
    if (fp->type == GTS_ERROR)
      return NULL;
    expr = g_string_new ("{");
    g_string_append (expr, s);
    g_free (s);
    g_string_append_c (expr, '}');
    if (is_expression)
      *is_expression = FALSE;
    return expr;
  }
  else {
    static gchar spaces[] = " \t\f\r";
    static gchar operators[] = "+-*/%<>=&^|?:!";
    gboolean is_constant = (fp->type == GTS_INT || fp->type == GTS_FLOAT);
    gint c, scope = 0;
    gchar * s;
    gchar empty[] = "", * comments = fp->comments;

    fp->comments = empty;
    expr = g_string_new (fp->token->str);
    s = expr->str;
    while (*s != '\0') {
      if (*s == '(') scope++;
      else if (*s == ')') scope--;
      s++;
    }
    if (fp->next_token != '\0')
      c = fp->next_token;
    else {
      if (fp->type != '(')
	c = ' ';
      else
	c = gts_file_getc (fp);
    }
    if (strlen (expr->str) == 1 && gfs_char_in_string (expr->str[0], operators))
      while (c != EOF && gfs_char_in_string (c, spaces)) {
	g_string_append_c (expr, c);
	c = gts_file_getc (fp);
      }
    while (c != EOF) {
      if (gfs_char_in_string (c, "{}\n")) {
	fp->next_token = c;
	break;
      }
      else if (scope > 0) {
	while (c != EOF && scope > 0) {
	  if (c == '(') scope++;
	  else if (c == ')') scope--;
	  g_string_append_c (expr, c);
	  c = gts_file_getc (fp);
	}
      }
      else if (gfs_char_in_string (c, spaces)) {
	while (c != EOF && gfs_char_in_string (c, spaces)) {
	  g_string_append_c (expr, c);
	  c = gts_file_getc (fp);
	}
	if (c == '(') {
	  if (is_constant) {
	    fp->next_token = c;
	    break;
	  }
	  scope++;
	  g_string_append_c (expr, c);
	  c = gts_file_getc (fp);
	}
	else {
	  if (!gfs_char_in_string (c, operators)) {
	    fp->next_token = c;
	    break;
	  }
	  is_constant = FALSE;
	  g_string_append_c (expr, c);
	  c = gts_file_getc (fp);
	  while (c != EOF && gfs_char_in_string (c, spaces)) {
	    g_string_append_c (expr, c);
	    c = gts_file_getc (fp);
	  }
	}
      }
      else if (gfs_char_in_string (c, operators)) {
	is_constant = FALSE;
	g_string_append_c (expr, c);
	c = gts_file_getc (fp);
	while (c != EOF && gfs_char_in_string (c, spaces)) {
	  g_string_append_c (expr, c);
	  c = gts_file_getc (fp);
	}
      }
      else {
	is_constant = FALSE;
	if (c == '(') scope++;
	else if (c == ')') scope--;
	if (scope < 0) {
	  fp->next_token = c;
	  break;
	}
	g_string_append_c (expr, c);
	c = gts_file_getc (fp);
      }
    }
    g_strchomp (expr->str);
    fp->comments = comments;
    return expr;
  }
}

static void function_read (GtsObject ** o, GtsFile * fp)
{
  GfsFunction * f = GFS_FUNCTION (*o);
  GfsSimulation * sim;
  GfsDomain * domain;

  if (GTS_OBJECT_CLASS (gfs_function_class ())->parent_class->read)
    (* GTS_OBJECT_CLASS (gfs_function_class ())->parent_class->read) (o, fp);
  if (fp->type == GTS_ERROR)
    return;

  sim = gfs_object_simulation (*o);
  domain = GFS_DOMAIN (sim);
  if (fp->type != GTS_INT && fp->type != GTS_FLOAT && fp->type != GTS_STRING &&
      fp->type != '(' && fp->type != '{') {
    gts_file_error (fp, "expecting an expression (val)");
    return;
  }

  if (fp->type == GTS_STRING && !f->spatial && !f->constant && strlen (fp->token->str) > 3) {
    if (!strcmp (&(fp->token->str[strlen (fp->token->str) - 4]), ".gts")) {
      if ((f->s = read_surface (fp->token->str, fp)) == NULL)
	return;
      f->sname = g_strdup (fp->token->str);
      gts_file_next_token (fp);
      return;
    }
    else if (!strcmp (&(fp->token->str[strlen (fp->token->str) - 4]), ".cgd")) {
      if ((f->g = gfs_cartesian_grid_read (fp->token->str, fp)) == NULL)
	return;
      if (!fit_index_dimension (f->g, f->index, fp))
	return;
      f->sname = g_strdup (fp->token->str);
      gts_file_next_token (fp);
      return;
    }   
  }

  if ((f->expr = gfs_function_expression (fp, &f->isexpr)) == NULL)
    return;

  if (f->isexpr) {
    if (fp->type == GTS_INT || fp->type == GTS_FLOAT) {
      if (!strcmp (fp->token->str, f->expr->str)) {
	f->val = atof (fp->token->str);
	f->constant = TRUE;
	gts_file_next_token (fp);
	return;
      }
    }
    else if (fp->type == GTS_STRING && !f->spatial && !f->constant) {
      if ((f->v = gfs_variable_from_name (domain->variables, f->expr->str)) ||
	  (f->dv = gfs_derived_variable_from_name (domain->derived_variables, f->expr->str))) {
	gts_file_next_token (fp);
	return;
      }
    }
  }

  gfs_module_new (f, fp->line);

  if (fp->type == GTS_ERROR)
    return;
  gts_file_next_token (fp);
}

static void function_write (GtsObject * o, FILE * fp)
{
  GfsFunction * f = GFS_FUNCTION (o);

  if (GTS_OBJECT_CLASS (gfs_function_class ())->parent_class->write)
    (* GTS_OBJECT_CLASS (gfs_function_class ())->parent_class->write) (o, fp);

  if (f->expr)
    fprintf (fp, " %s", f->expr->str);
  else if (f->v)
    fprintf (fp, " %s", f->v->name);
  else if (f->s || f->g)
    fprintf (fp, " %s", f->sname);
  else
    fprintf (fp, " %g", f->val);
}

static void function_destroy (GtsObject * object)
{
  GfsFunction * f = GFS_FUNCTION (object);

  if (f->module)
    gfs_module_unref (f->module, f);
  if (f->expr) g_string_free (f->expr, TRUE);
  if (f->s) {
    gts_object_destroy (GTS_OBJECT (f->s));
    g_free (f->sname);
  }
  if (f->g) {
    gts_object_destroy (GTS_OBJECT (f->g));
    g_free (f->sname);
  }
  g_free (f->var);
  g_free (f->dvar);

  (* GTS_OBJECT_CLASS (gfs_function_class ())->parent_class->destroy) 
    (object);
}

static void gfs_function_class_init (GfsFunctionClass * klass)
{
  GTS_OBJECT_CLASS (klass)->read = function_read;
  GTS_OBJECT_CLASS (klass)->write = function_write;
  GTS_OBJECT_CLASS (klass)->destroy = function_destroy;
}

GfsFunctionClass * gfs_function_class (void)
{
  static GfsFunctionClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_function_info = {
      "GfsFunction",
      sizeof (GfsFunction),
      sizeof (GfsFunctionClass),
      (GtsObjectClassInitFunc) gfs_function_class_init,
      (GtsObjectInitFunc) NULL,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gts_object_class ()),
				  &gfs_function_info);
  }

  return klass;
}

/**
 * gfs_function_new:
 * @klass: a #GfsFunctionClass.
 * @val: a value.
 *
 * Returns: a new #GfsFunction with constant value @val.
 */
GfsFunction * gfs_function_new (GfsFunctionClass * klass, 
				gdouble val)
{
  GfsFunction * object;

  object = GFS_FUNCTION (gts_object_new (GTS_OBJECT_CLASS (klass)));
  object->val = val;

  return object;
}

GfsFunction * gfs_function_new_from_variable (GfsFunctionClass * klass, 
					      GfsVariable * v)
{
  GfsFunction * object;

  g_return_val_if_fail (v != NULL, NULL);

  object = GFS_FUNCTION (gts_object_new (GTS_OBJECT_CLASS (klass)));
  object->v = v;

  return object;
}

/**
 * gfs_function_set_units:
 * @f: a #GfsFunction.
 * @units: the units of @f.
 *
 * Sets the units of @f.
 */
void gfs_function_set_units (GfsFunction * f, 
			     gdouble units)
{
  g_return_if_fail (f != NULL);
  f->units = units;
}

static gdouble interpolated_value (GfsFunction * f, FttVector * p)
{
  GtsPoint q;
  GtsFace * t;

  gfs_simulation_map_inverse (gfs_object_simulation (f), p);
  q.x = p->x; q.y = p->y;
  t = gts_point_locate (&q, f->s, NULL);
  if (t == NULL)
    return 0.;
  gts_triangle_interpolate_height (GTS_TRIANGLE (t), &q);
  return q.z;
}

/**
 * gfs_function_description:
 * @f: a #GfsFunction.
 * @truncate: whether to truncate long descriptions.
 *
 * Returns: a newly allocated string describing @f.
 */
gchar * gfs_function_description (GfsFunction * f,
				  gboolean truncate)
{
  gchar * s;

  g_return_val_if_fail (f != NULL, NULL);

  if (f->s)
    s = g_strdup (f->sname);
  else if (f->v)
    s = g_strdup (f->v->name);
  else if (f->expr) {
    s = g_strdup (f->expr->str);    
    if (truncate) {
      gchar * c = s;
      guint n = 0;
      
      while (*c != '\0' && !isspace (*c))
	c++;
      while (*c != '\0' && n < 3) {
	*c = '.';
	c++; n++;
      }
      *c = '\0';
    }
  }
  else
    s = g_strdup_printf ("%g", f->val);
  return s;
}

static gdouble adimensional_value (GfsFunction * f, gdouble v)
{
  gdouble L;
  if (v == GFS_NODATA || f->units == 0. || 
      (L = gfs_object_simulation (f)->physical_params.L) == 1.)
    return v;
  return v*pow (L, - f->units);
}

/**
 * gfs_function_value:
 * @f: a #GfsFunction.
 * @cell: a #FttCell or %NULL.
 *
 * Returns: the value of function @f in @cell.
 */
gdouble gfs_function_value (GfsFunction * f, FttCell * cell)
{
  g_return_val_if_fail (f != NULL, 0.);
  g_assert (!pending_functions);

  gdouble dimensional;
  if (f->s) {
    FttVector p;
    gfs_cell_cm (cell, &p);
    dimensional = interpolated_value (f, &p);
  }
  else if (f->g) {
    FttVector p = {0.,0.,0.};
    if (cgd_is_spatial (f))
      gfs_cell_cm (cell, &p);
    dimensional = interpolated_cgd (f, &p);
  }
  else if (f->v)
    dimensional = gfs_dimensional_value (f->v, GFS_VALUE (cell, f->v));
  else if (f->dv)
    dimensional = (* (GfsFunctionDerivedFunc) f->dv->func) (cell, NULL,
							    gfs_object_simulation (f),
							    f->dv->data);
  else if (f->f)
    dimensional = (* f->f) (cell, NULL, gfs_object_simulation (f), f->var, f->dvar);
  else
    dimensional = f->val;
  return adimensional_value (f, dimensional);
}

/**
 * gfs_function_face_value:
 * @f: a #GfsFunction.
 * @fa: a #FttCellFace.
 *
 * Returns: the value of function @f at the center of face @fa.
 */
gdouble gfs_function_face_value (GfsFunction * f, FttCellFace * fa)
{
  g_return_val_if_fail (f != NULL, 0.);
  g_return_val_if_fail (fa != NULL, 0.);
  g_assert (!pending_functions);

  gdouble dimensional;
  if (f->s) {
    FttVector p;
    ftt_face_pos (fa, &p);
    dimensional = interpolated_value (f, &p);
  }
  else if (f->g) {
    FttVector p = {0.,0.,0.};
    if (cgd_is_spatial (f))
      ftt_face_pos (fa, &p);
    dimensional = interpolated_cgd (f, &p);
  }
  else if (f->v)
    dimensional = gfs_dimensional_value (f->v, gfs_face_interpolated_value_generic (fa, f->v));
  else if (f->dv)
    dimensional = (* (GfsFunctionDerivedFunc) f->dv->func) (NULL, fa,
							    gfs_object_simulation (f),
							    f->dv->data);
  else if (f->f)
    dimensional = (* f->f) (NULL, fa, gfs_object_simulation (f), f->var, f->dvar);
  else
    dimensional = f->val;
  return adimensional_value (f, dimensional);
}

/**
 * gfs_function_set_constant_value:
 * @f: a #GfsFunction.
 * @val: the value.
 *
 * Sets the value of the constant function @f to @val.
 */
void gfs_function_set_constant_value (GfsFunction * f, gdouble val)
{
  g_return_if_fail (f != NULL);
  g_return_if_fail (!f->f && !f->s && !f->v && !f->dv);

  f->val = val;
  f->constant = TRUE;
}

/**
 * gfs_function_get_constant_value:
 * @f: a #GfsFunction.
 *
 * Returns: the value of function @f if @f is constant, G_MAXDOUBLE
 * otherwise.
 */
gdouble gfs_function_get_constant_value (GfsFunction * f)
{
  g_return_val_if_fail (f != NULL, G_MAXDOUBLE);
  g_assert (!pending_functions);

  if (f->f || f->s || f->v || f->dv)
    return G_MAXDOUBLE;
  else
    return adimensional_value (f, f->val);
}

/**
 * gfs_function_is_constant:
 * @f: a #GfsFunction.
 *
 * Returns: %TRUE if @f is a constant, %FALSE otherwise.
 */
gboolean gfs_function_is_constant (const GfsFunction * f)
{
  g_return_val_if_fail (f != NULL, FALSE);

  return f->constant;
}

/**
 * gfs_function_get_variable:
 * @f: a #GfsFunction.
 *
 * Returns: the variable containing the value of @f if @f is a simple
 * variable, NULL otherwise.
 */
GfsVariable * gfs_function_get_variable (GfsFunction * f)
{
  g_return_val_if_fail (f != NULL, NULL);

  return f->v;
}

/**
 * gfs_function_read:
 * @f: a #GfsFunction.
 * @domain: a #GfsDomain.
 * @fp: a #GtsFile.
 *
 * Calls the read() method of @f.
 */
void gfs_function_read (GfsFunction * f, gpointer domain, GtsFile * fp)
{
  GtsObject * o = (GtsObject *) f;

  g_return_if_fail (f != NULL);
  g_return_if_fail (domain != NULL);
  g_return_if_fail (fp != NULL);

  GTS_OBJECT (f)->reserved = domain;
  (* GTS_OBJECT (f)->klass->read) (&o, fp);
}

/**
 * gfs_function_write:
 * @f: a #GfsFunction.
 * @fp: a file pointer.
 *
 * Calls the write() method of @f.
 */
void gfs_function_write (GfsFunction * f, FILE * fp)
{
  g_return_if_fail (f != NULL);
  g_return_if_fail (fp != NULL);

  (* GTS_OBJECT (f)->klass->write) (GTS_OBJECT (f), fp);
}

/** \endobject{GfsFunction} */

/**
 * Functions of (x,y,z,t).
 * \beginobject{GfsFunctionSpatial}
 */

static void gfs_function_spatial_init (GfsFunction * f)
{
  f->spatial = TRUE;
}

GfsFunctionClass * gfs_function_spatial_class (void)
{
  static GfsFunctionClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_function_info = {
      "GfsFunctionSpatial",
      sizeof (GfsFunction),
      sizeof (GfsFunctionClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_function_spatial_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_function_class ()),
				  &gfs_function_info);
  }

  return klass;
}

/**
 * gfs_function_spatial_value:
 * @f: a #GfsFunction.
 * @p: a #FttVector.
 *
 * Returns: the value of function @f at location @p.
 */
gdouble gfs_function_spatial_value (GfsFunction * f, const FttVector * p)
{
  g_return_val_if_fail (f != NULL, 0.);
  g_return_val_if_fail (GFS_IS_FUNCTION_SPATIAL (f), 0.);
  g_return_val_if_fail (p != NULL, 0.);
  g_assert (!pending_functions);

  gdouble dimensional;  
  if (f->f) {
    GfsSimulation * sim = gfs_object_simulation (f);
    FttVector q = *p;
    if (!f->nomap)
      gfs_simulation_map_inverse (sim, &q);
    dimensional = (* (GfsFunctionSpatialFunc) f->f) (q.x, q.y, q.z, sim->time.t);
  }
  else
    dimensional = f->val;
  return adimensional_value (f, dimensional);
}

GfsFunction * gfs_function_spatial_new (GfsFunctionClass * klass, 
					GfsFunctionSpatialFunc func)
{
  GfsFunction * object;

  g_return_val_if_fail (func != NULL, NULL);

  object = GFS_FUNCTION (gts_object_new (GTS_OBJECT_CLASS (klass)));
  object->f = (GfsFunctionFunc) func;

  return object;
}

/** \endobject{GfsFunctionSpatial} */

/**
 * Abstract class for coordinate transformations.
 * \beginobject{GfsFunctionMap}
 */

static void gfs_function_map_init (GfsFunction * f)
{
  f->nomap = TRUE;
}

GfsFunctionClass * gfs_function_map_class (void)
{
  static GfsFunctionClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_function_info = {
      "GfsFunctionMap",
      sizeof (GfsFunction),
      sizeof (GfsFunctionClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_function_map_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_function_spatial_class ()),
				  &gfs_function_info);
  }

  return klass;
}

/** \endobject{GfsFunctionMap} */

/**
 * Symbolic constants.
 * \beginobject{GfsFunctionConstant}
 */

static void gfs_function_constant_init (GfsFunction * f)
{
  f->constant = TRUE;
}

GfsFunctionClass * gfs_function_constant_class (void)
{
  static GfsFunctionClass * klass = NULL;

  if (klass == NULL) {
    GtsObjectClassInfo gfs_function_info = {
      "GfsFunctionConstant",
      sizeof (GfsFunction),
      sizeof (GfsFunctionClass),
      (GtsObjectClassInitFunc) NULL,
      (GtsObjectInitFunc) gfs_function_constant_init,
      (GtsArgSetFunc) NULL,
      (GtsArgGetFunc) NULL
    };
    klass = gts_object_class_new (GTS_OBJECT_CLASS (gfs_function_class ()),
				  &gfs_function_info);
  }

  return klass;
}

/**
 * gfs_read_constant:
 * @fp: a #GtsFile.
 * @domain: a #GfsDomain.
 *
 * Reads a constant value from @fp.
 *
 * Returns: the value of the constant or G_MAXDOUBLE if an error
 * occured.
 */
gdouble gfs_read_constant (GtsFile * fp, gpointer domain)
{
  g_return_val_if_fail (fp != NULL, G_MAXDOUBLE);
  g_return_val_if_fail (domain != NULL, G_MAXDOUBLE);

  GfsFunction * f = gfs_function_new (gfs_function_constant_class (), 0.);
  gfs_function_read (f, domain, fp);
  gfs_pending_functions_compilation (fp);
  if (fp->type == GTS_ERROR)
    return G_MAXDOUBLE;
  gdouble val = gfs_function_get_constant_value (f);
  gts_object_destroy (GTS_OBJECT (f));
  if (val == G_MAXDOUBLE)
    gts_file_error (fp, "expecting a constant");
  return val;
}

/** \endobject{GfsFunctionConstant} */

/**
 * gfs_object_class_from_name:
 * @name: the name of the class.
 *
 * Looks for a class called @name. If none is found append the "Gfs"
 * prefix and look again.
 *
 * Returns: the class or %NULL if none was found.
 */
GtsObjectClass * gfs_object_class_from_name (const gchar * name)
{
  GtsObjectClass * klass;

  g_return_val_if_fail (name != NULL, NULL);

  if ((klass = gts_object_class_from_name (name)))
    return klass;
  /* for backward parameter file compatibility */
  if (!strcmp (name, "GtsSurfaceFile"))
    return GTS_OBJECT_CLASS (gfs_solid_class ());
  gchar * ename = g_strconcat ("Gfs", name, NULL);
  klass = gts_object_class_from_name (ename);
  g_free (ename);
  return klass;
}

static void eigsrt (gdouble d[FTT_DIMENSION], gdouble v[FTT_DIMENSION][FTT_DIMENSION])
{
  gint k, j, i;
  gdouble p;

  for (i = 0; i < FTT_DIMENSION - 1; i++) {
    p = d[k = i];

    for (j = i + 1; j < FTT_DIMENSION; j++)
      if (d[j] >= p) 
	p = d[k = j];
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (j = 0; j < FTT_DIMENSION; j++) {
	p = v[j][i];
	v[j][i] = v[j][k];
	v[j][k] = p;
      }
    }
  }
}

#define ROTATE(a,i,j,k,l) {g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);}

/**
 * gfs_eigenvalues:
 * @a: a symmetric matrix.
 * @d: a vector.
 * @v: another matrix.
 *
 * Fills @d (resp. @v) with the eigenvalues (resp. eigenvectors) of
 * matrix @a.
 */
void gfs_eigenvalues (gdouble a[FTT_DIMENSION][FTT_DIMENSION],
		      gdouble d[FTT_DIMENSION],
		      gdouble v[FTT_DIMENSION][FTT_DIMENSION])
{
  gint j, iq, ip, i;
  gdouble tresh, theta, tau, t, sm, s, h, g, c, b[FTT_DIMENSION], z[FTT_DIMENSION];

  for (ip = 0; ip < FTT_DIMENSION; ip++) {
    for (iq = 0; iq < FTT_DIMENSION; iq++)
      v[ip][iq] = 0.0;
    v[ip][ip] = 1.0;
  }

  for (ip = 0; ip < FTT_DIMENSION; ip++) {
    b[ip] = d[ip] = a[ip][ip];
    z[ip] = 0.0;
  }

  for (i = 1; i <= 50; i++) {
    sm = 0.0;
    for (ip = 0; ip < FTT_DIMENSION - 1; ip++) {
      for (iq = ip + 1; iq < FTT_DIMENSION; iq++)
	sm += fabs (a[ip][iq]);
    }
    if (sm == 0.0) {
      eigsrt (d, v);
      return;
    }
    if (i < 4)
      tresh = 0.2*sm/(FTT_DIMENSION*FTT_DIMENSION);
    else
      tresh = 0.0;
    for (ip = 0; ip < FTT_DIMENSION - 1; ip++) {
      for (iq = ip + 1; iq < FTT_DIMENSION; iq++) {
	g = 100.0*fabs (a[ip][iq]);
	if (i > 4 && fabs(d[ip]) + g == fabs(d[ip]) && fabs(d[iq]) + g == fabs(d[iq]))
	  a[ip][iq] = 0.0;
	else if (fabs (a[ip][iq]) > tresh) {
	  h = d[iq] - d[ip];
	  if (fabs(h) + g == fabs(h))
	    t = a[ip][iq]/h;
	  else {
	    theta = 0.5*h/a[ip][iq];
	    t = 1.0/(fabs (theta) + sqrt (1.0 + theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c = 1.0/sqrt (1 + t*t);
	  s = t*c;
	  tau = s/(1.0 + c);
	  h = t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq] = 0.0;
	  for (j = 0; j <= ip - 1; j++)
	    ROTATE (a, j, ip, j, iq);
	  for (j = ip + 1; j <= iq - 1; j++)
	    ROTATE (a, ip, j, j, iq);
	  for (j = iq + 1; j < FTT_DIMENSION; j++)
	    ROTATE(a, ip, j, iq, j);
	  for (j = 0; j < FTT_DIMENSION; j++)
	    ROTATE(v, j, ip, j, iq);
	}
      }
    }
    for (ip = 0; ip < FTT_DIMENSION; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }
  /* Too many iterations */
  for (i = 0; i < FTT_DIMENSION; i++) {
    for (j = 0; j < FTT_DIMENSION; j++)
      fprintf (stderr, "%10.3g ", a[i][j]);
    fprintf (stderr, "\n");
  }
  g_assert_not_reached ();
}

/**
 * gfs_matrix_inverse:
 * @m: a square matrix.
 * @n: size of the matrix.
 * @pivmin: the minimum value of the pivoting coefficient.
 *
 * Replaces @m with its inverse.
 *
 * Returns: 0. if the inversion encounters a pivot coefficient smaller
 * than or equal to @pivmin (i.e. @m is non-invertible), the minimum
 * absolute value of the pivoting coefficient otherwise.
 */
gdouble gfs_matrix_inverse (gdouble ** m, guint n, gdouble pivmin)
{
  gint * indxc, * indxr, * ipiv;
  gint i, icol = 0, irow = 0, j, k, l, ll;
  gdouble big, dum, pivinv, temp, minpiv = G_MAXDOUBLE;

  g_return_val_if_fail (m != NULL, 0.);
  g_return_val_if_fail (pivmin >= 0., 0.);

  indxc = g_malloc (sizeof (gint)*n);
  indxr = g_malloc (sizeof (gint)*n);
  ipiv = g_malloc (sizeof (gint)*n);
  
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

  for (j = 0; j < n; j++)
    ipiv[j] = -1;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
	for (k = 0; k < n; k++) {
	  if (ipiv[k] == -1) {
	    if (fabs (m[j][k]) >= big) {
	      big = fabs (m[j][k]);
	      irow = j;
	      icol = k;
	    }
	  }
	}
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++) 
	SWAP (m[irow][l], m[icol][l]);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin) {
      g_free (indxc);
      g_free (indxr);
      g_free (ipiv);
      return 0.;
    }
    if (fabs (m[icol][icol]) < minpiv)
      minpiv = fabs (m[icol][icol]);
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
	dum = m[ll][icol];
	m[ll][icol] = 0.0;
	for (l = 0; l < n; l++)
	  m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
	SWAP (m[k][indxr[l]], m[k][indxc[l]]);
  }
  g_free (indxc);
  g_free (indxr);
  g_free (ipiv);
  return minpiv;
}

/**
 * gfs_matrix_new:
 * @n: the size of the matrix.
 * @p: the size of the matrix.
 * @size: the size of the matrix elements.
 *
 * The matrix elements are initialised to zero.
 *
 * Returns: a newly allocated matrix.
 */
gpointer gfs_matrix_new (guint n, guint p, guint size)
{
  guint i;
  gpointer * m;
  gchar * a;
  
  g_return_val_if_fail (n > 0, NULL);
  g_return_val_if_fail (p > 0, NULL);
  g_return_val_if_fail (size > 0, NULL);

  m = g_malloc (n*sizeof (gpointer));
  a = g_malloc0 (n*p*size);
  for (i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return m;
}

/**
 * gfs_matrix_free:
 * @m: a matrix allocated with gfs_matrix_new().
 *
 * Frees the memory occupied by @m.
 */
void gfs_matrix_free (gpointer m)
{
  g_return_if_fail (m != NULL);

  g_free (((gpointer *) m)[0]);
  g_free (m);
}

/**
 * gfs_clock_new:
 *
 * Returns: a new #GfsClock.
 */
GfsClock * gfs_clock_new (void)
{
  GfsClock * t = g_malloc (sizeof (GfsClock));

  t->start = -1;
  t->started = FALSE;
  return t;
}

/**
 * gfs_clock_start:
 * @t: a #GfsClock.
 *
 * Starts clock @t.
 */
void gfs_clock_start (GfsClock * t)
{
  struct tms tm;

  g_return_if_fail (t != NULL);
  g_return_if_fail (!t->started);

  if (times (&tm) == (clock_t) -1)
    g_warning ("cannot read clock");
  t->start = tm.tms_utime;
  t->started = TRUE;
}

/**
 * gfs_clock_stop:
 * @t: a #GfsClock.
 *
 * Stops clock @t.
 */
void gfs_clock_stop (GfsClock * t)
{
  struct tms tm;

  g_return_if_fail (t != NULL);
  g_return_if_fail (t->started);

  if (times (&tm) == (clock_t) -1)
    g_warning ("cannot read clock");
  t->stop = tm.tms_utime;
  t->started = FALSE;
}

/**
 * gfs_clock_elapsed:
 * @t: a #GfsClock.
 *
 * Returns: the time elapsed in seconds since @t was started.
 */
gdouble gfs_clock_elapsed (GfsClock * t)
{
  g_return_val_if_fail (t != NULL, 0.);
  g_return_val_if_fail (t->start >= 0, 0.);

  if (t->started == FALSE)
    return (t->stop - t->start)/(gdouble) sysconf (_SC_CLK_TCK);
  else {
    struct tms tm;
    if (times (&tm) == (clock_t) -1)
      g_warning ("cannot read clock");
    return (tm.tms_utime - t->start)/(gdouble) sysconf (_SC_CLK_TCK);
  }
}

/**
 * gfs_clock_destroy:
 * @t: a #GfsClock.
 *
 * Destroys the clock, freeing the memory allocated for it.
 */
void gfs_clock_destroy (GfsClock * t)
{
  g_return_if_fail (t != NULL);

  g_free (t);
}

/**
 * gfs_union_open:
 * @fp: a file pointer.
 * @rank: the rank of the current parallel process.
 * @file: a #GfsUnionFile.
 *
 * Opens a "parallel" file which serialises multiple parallel (write)
 * accesses to the file pointed to by @fp.
 *
 * This file must be closed with gfs_union_close().
 *
 * Returns: a "parallel" file pointer associated with @fp.
 */
FILE * gfs_union_open (FILE * fp, int rank, GfsUnionFile * file)
{
  g_return_val_if_fail (fp != NULL, NULL);
  g_return_val_if_fail (file != NULL, NULL);

  if (rank <= 0) /* master */
    return fp;
  else { /* slaves */
#ifdef HAVE_MPI
    MPI_Status status;
    int pe;
    MPI_Recv (&pe, 1, MPI_INT, 0, rank, MPI_COMM_WORLD, &status);
    g_assert (rank == pe);
#endif /* HAVE_MPI */
    file->fp = open_memstream (&file->buf, &file->len);
    if (file->fp == NULL)
      g_error ("gfs_union_open(): could not open_memstream:\n%s", strerror (errno));
    return file->fp;
  }
}

/**
 * gfs_union_close:
 * @fp: a file pointer.
 * @rank: the rank of the current parallel process.
 * @file: a #GfsUnionFile returned by a call to gfs_union_open().
 *
 * Closes a "parallel" file previously opened using gfs_union_open().
 */
void gfs_union_close (FILE * fp, int rank, GfsUnionFile * file)
{
  g_return_if_fail (fp != NULL);
  g_return_if_fail (file != NULL);

  if (rank == 0) { /* master */
#ifdef HAVE_MPI
    int pe, npe;
    MPI_Comm_size (MPI_COMM_WORLD, &npe);
    for (pe = 1; pe < npe; pe++) {
      MPI_Send (&pe, 1, MPI_INT, pe, pe, MPI_COMM_WORLD);
      MPI_Status status;
      long length;
      MPI_Recv (&length, 1, MPI_LONG, pe, pe, MPI_COMM_WORLD, &status);
      /*      fprintf (stderr, "receiving %ld bytes from PE %d\n", length, pe); */
      if (length > 0) {
	char * buf = g_malloc (length);
	MPI_Recv (buf, length, MPI_BYTE, pe, pe + 1, MPI_COMM_WORLD, &status);
	int rcvcount;
	MPI_Get_count (&status, MPI_BYTE, &rcvcount);
	fwrite (buf, 1, rcvcount, fp);
	g_free (buf);
      }
    }
#endif /* HAVE_MPI */
  }
  else if (rank > 0) { /* slaves */
    fclose (file->fp);
    long length = file->len;
#ifdef HAVE_MPI
    MPI_Send (&length, 1, MPI_LONG, 0, rank, MPI_COMM_WORLD);
    if (length > 0)
      MPI_Send (file->buf, length, MPI_BYTE, 0, rank + 1, MPI_COMM_WORLD);
#endif /* HAVE_MPI */
    if (length > 0)
      g_free (file->buf);
  }
}

static GfsFormat * format_new (const gchar * s, 
			       guint len, 
			       GfsFormatType t)
{
  GfsFormat * f = g_malloc (sizeof (GfsFormat));
  
  f->s = g_strndup (s, len);
  f->t = t;
  
  return f;
}

static void format_destroy (GfsFormat * f)
{
  g_free (f->s);
  g_free (f);
}

/**
 * gfs_format_new:
 * @format: a string.
 * @fp: a #GtsFile or %NULL.
 * @dynamic: set to %TRUE if the format is time-dependent.
 * @parallel: set to %TRUE if the format is PID-dependent.
 *
 * If @fp is not %NULL and and error occurs, an error message is set
 * in @fp.
 *
 * Returns: a list of #GfsFormat composing @format.
 */
GSList * gfs_format_new (const gchar * format,
			 GtsFile * fp,
			 gboolean * dynamic,
			 gboolean * parallel)
{
  g_return_val_if_fail (format != NULL, NULL);

  GSList * formats = NULL;
  const gchar * c, * start;
  guint len;

  start = c = format;
  while (*c != '\0') {
    if (*c == '%') {
      const gchar * startf = c, * prev = c;
	
      len = startf - start;
      if (len > 0)
	formats = g_slist_prepend (formats, format_new (start, len, GFS_NONE_FORMAT));
	
      len = 1;
      c++;
      while (*c != '\0' && !gfs_char_in_string (*c, "diouxXeEfFgGaAcsCSpn%")) {
	prev = c;
	c++;
	len++;
      }
      len++;
      if (*c == '%')
	formats = g_slist_prepend (formats, format_new ("%", 1, GFS_NONE_FORMAT));
      else if (gfs_char_in_string (*c, "diouxXc")) {
	if (*prev == 'l') {
	  formats = g_slist_prepend (formats, format_new (startf, len, GFS_ITER_FORMAT));
	  if (dynamic)
	    *dynamic = TRUE;
	}
	else {
	  formats = g_slist_prepend (formats, format_new (startf, len, GFS_PID_FORMAT));
	  if (parallel)
	    *parallel = TRUE;
	}
      }
      else if (gfs_char_in_string (*c, "eEfFgGaA")) {
	formats = g_slist_prepend (formats, format_new (startf, len, GFS_TIME_FORMAT));
	if (dynamic)
	  *dynamic = TRUE;
      }
      else {
	if (fp)
	  gts_file_error (fp, 
			  "unknown conversion specifier `%c' of format `%s'",
			  *c, format);
	return NULL;
      }
      start = c;
      start++;
    }
    c++;
  }
  len = c - start;
  if (len > 0)
    formats = g_slist_prepend (formats, format_new (start, len, GFS_NONE_FORMAT));
  formats = g_slist_reverse (formats);
  
  return formats;
}

/**
 * gfs_format_destroy:
 * @f: a list of #GfsFormat.
 *
 * Frees all memory allocated for @f.
 */
void gfs_format_destroy (GSList * f)
{
  g_slist_foreach (f, (GFunc) format_destroy, NULL);
  g_slist_free (f);
}

/**
 * gfs_format_string:
 * @format: a GSList of #GfsFormat.
 * @pid: the PID.
 * @niter: number of iterations done in the simulation.
 * @time: simulation time.
 *
 * Returns: a newly-allocated string (file name) built from the
 * GfsFormat contained in GSList.  It typically includes informations
 * on the PID, the time and/or the number of iterations done in the
 * simulation.
 */
gchar * gfs_format_string (GSList * format, 
			   gint pid, 
			   guint niter,
			   gdouble time)
{
  gchar * s = g_strdup ("");

  while (format) {
    GfsFormat * f = format->data;
    gchar * s1, * s2 = NULL;

    switch (f->t) {
    case GFS_NONE_FORMAT:
      s2 = g_strconcat (s, f->s, NULL);
      break;
    case GFS_PID_FORMAT:
      s1 = g_strdup_printf (f->s, pid);
      s2 = g_strconcat (s, s1, NULL);
      g_free (s1);
      break;
    case GFS_ITER_FORMAT:
      s1 = g_strdup_printf (f->s, niter);
      s2 = g_strconcat (s, s1, NULL);
      g_free (s1);
      break;
    case GFS_TIME_FORMAT:
      s1 = g_strdup_printf (f->s, time);
      s2 = g_strconcat (s, s1, NULL);
      g_free (s1);
      break;
    default:
      g_assert_not_reached ();
    }
    g_free (s);
    s = s2;

    format = format->next;
  }

  return s;
}

/**
 * gfs_format_time_value:
 * @format: a list of #GfsFormat.
 * @string: a string formatted according to @format.
 *
 * Returns: the value of the time or iteration contained in @s or
 * %G_MAXDOUBLE if @s is not formatted according to @format.
 *
 */
gdouble gfs_format_time_value (GSList * format, const gchar * string)
{
  gdouble val = G_MAXDOUBLE, tmp;

  g_return_val_if_fail (string != NULL, val);

  gchar * copy = g_strdup (string), * s = copy;
  while (format) {
    GfsFormat * f = format->data;
    gchar * c, c1;

    switch (f->t) {
    case GFS_NONE_FORMAT:
      c = f->s;
      while (*c != '\0' && *c == *s) { c++; s++; }
      if (*c != '\0') {
	g_free (copy);
	return val;
      }
      break;
    case GFS_ITER_FORMAT:
      c = s;
      while (gfs_char_in_string (*s, "0123456789")) s++;
      c1 = *s; *s = '\0'; tmp = atoi (c); *s = c1;
      if (val != G_MAXDOUBLE && tmp != val) {
	g_free (copy);
	return G_MAXDOUBLE;
      }
      val = tmp;
      break;
    case GFS_TIME_FORMAT:
      c = s;
      while (gfs_char_in_string (*s, "0123456789eE-+.")) s++;
      c1 = *s; *s = '\0'; tmp = atof (c); *s = c1;
      if (val != G_MAXDOUBLE && tmp != val) {
	g_free (copy);
	return G_MAXDOUBLE;
      }
      val = tmp;
      break;
    default:
      g_assert_not_reached ();
    }
    format = format->next;
  }
  g_free (copy);
  return val;
}

/**
 * gfs_cell_message:
 * @cell: a #FttCell.
 * @format: a string format.
 * ...: arguments for format.
 *
 * Logs a message preceded by the pointer, position and level of @cell.
 */
void gfs_cell_message (const FttCell * cell, 
		       const gchar * format,
		       ...)
{
  g_return_if_fail (cell != NULL);
  g_return_if_fail (format != NULL);

  FttVector p;
  ftt_cell_pos (cell, &p);
  gchar * s = g_strdup_printf ("%p:(%g,%g,%g):%d", cell, p.x, p.y, p.z, ftt_cell_level (cell));
  va_list ap;
  va_start (ap, format);
  gchar * s1 = g_strdup_vprintf (format, ap);
  g_message ("%s\n%s", s, s1);
  g_free (s);
  g_free (s1);
}

static gboolean GfsDebug = FALSE;

/**
 * gfs_debug:
 * @format: a string format.
 * ...: arguments for format.
 *
 * Logs a debugging message (only when gfs_debug_enabled() is set to
 * %TRUE).
 */
void gfs_debug (const gchar * format,
		...)
{
  if (GfsDebug) {
    g_return_if_fail (format != NULL);
    va_list ap;
    va_start (ap, format);    
    g_logv (G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG, format, ap);
  }  
}

/**
 * gfs_debug_enabled:
 * @enabled: whether to enable debug message logging.
 * 
 * Enables or disables debug message logging (see also gfs_debug()).
 */
void gfs_debug_enabled (gboolean enabled)
{
  GfsDebug = enabled;
}

/**
 * gfs_read_vector:
 * @fp: a #GtsFile.
 * @component: an array of vector component names.
 *
 * Reads a vector from @fp and fills @component. Note that the
 * component strings need to be freed after use.
 *
 * Returns: %TRUE if all vector components have been read, %FALSE
 * otherwise, in which case the error is reported in @fp.
 */
gboolean gfs_read_vector (GtsFile * fp, gchar * component[FTT_DIMENSION])
{
  g_return_val_if_fail (fp != NULL, FALSE);
  g_return_val_if_fail (component != NULL, FALSE);

  if (fp->type != '(') {
    gts_file_error (fp, "expecting an opening bracket '('");
    return FALSE;
  }
  gint c = gts_file_getc (fp), parlevel = 0, n = 0;
  GString * s = g_string_new ("");
  while ((parlevel > 0 || c != ')') && c != EOF) {
    if (c == ',' && parlevel == 0) {
      /* end of component */
      if (n == FTT_DIMENSION) {
	gts_file_error (fp, "too many vector components");
	gint i;
	for (i = 0; i < n; i++)
	  g_free (component[i]);
	g_string_free (s, TRUE);
	return FALSE;
      }
      component[n++] = s->str;
      g_string_free (s, FALSE);
      s = g_string_new ("");
    }
    else {
      /* expression */
      g_string_append_c (s, c);
      if (c == '(' || c == '[' || c == '{')
	parlevel++;
      else if (c == ')' || c == ']' || c == '}')
	parlevel--;
    }
    c = gts_file_getc (fp);
  }
  if (c != ')') {
    gts_file_error (fp, "parse error (missing closing bracket ')'?)"); 
    gint i;
    for (i = 0; i < n; i++)
      g_free (component[i]);
    g_string_free (s, TRUE);
    return FALSE;
  }
  if (n != FTT_DIMENSION - 1) {
    gts_file_error (fp, "not enough vector components");
    gint i;
    for (i = 0; i < n; i++)
      g_free (component[i]);
    g_string_free (s, TRUE);
    return FALSE;
  }
  component[n++] = s->str;
  g_string_free (s, FALSE);
  gts_file_next_token (fp);
  return TRUE;
}

/**
 * gfs_read_variable_vector:
 * @fp: a #GtsFile.
 * @vector: an array of vector component variables.
 * @domain: the #GfsDomain to which @vector belongs.
 *
 * Reads a vector from @fp and fills @vector.
 *
 * Returns: %TRUE if all vector components have been read, %FALSE
 * otherwise, in which case the error is reported in @fp.
 */
gboolean gfs_read_variable_vector (GtsFile * fp, 
				   GfsVariable * vector[FTT_DIMENSION], 
				   GfsDomain * domain)
{
  g_return_val_if_fail (fp != NULL, FALSE);
  g_return_val_if_fail (vector != NULL, FALSE);
  g_return_val_if_fail (domain != NULL, FALSE);

  gchar * component[FTT_DIMENSION];
  if (!gfs_read_vector (fp, component))
    return FALSE;

  gint i;
  gboolean status = TRUE;
  for (i = 0; i < FTT_DIMENSION && status; i++) {
    if (!(vector[i] = gfs_variable_from_name (domain->variables, component[i]))) {
      gts_file_error (fp, "unknown variable '%s'", component[i]);
      status = FALSE;
    }
    else if (vector[i]->component != i) {
      gts_file_error (fp, "variable '%s' is not the correct vector component", component[i]);
      status = FALSE;
    }
  }
  for (i = 0; i < FTT_DIMENSION; i++)
    g_free (component[i]);
  return status;
}

/**
 * gfs_read_function_vector:
 * @fp: a #GtsFile.
 * @vector: an array of vector component variables.
 * @function: an array of vector component #GfsFunction corresponding
 * to @vector.
 * @sim: the #GfsSimulation to which @vector belongs.
 *
 * Reads a vector from @fp and fills @vector.
 *
 * Returns: %TRUE if all vector components have been read, %FALSE
 * otherwise, in which case the error is reported in @fp.
 */
gboolean gfs_read_function_vector (GtsFile * fp, 
				   GfsVariable * vector[FTT_DIMENSION], 
				   GfsFunction * function[FTT_DIMENSION],
				   gpointer sim)
{
  g_return_val_if_fail (fp != NULL, FALSE);
  g_return_val_if_fail (vector != NULL, FALSE);
  g_return_val_if_fail (function != NULL, FALSE);
  g_return_val_if_fail (sim != NULL, FALSE);

  gchar * component[FTT_DIMENSION];
  if (!gfs_read_vector (fp, component))
    return FALSE;

  gboolean status = TRUE;
  gint i;
  for (i = 0; i < FTT_DIMENSION && status; i++) {
    function[i] = gfs_function_new (gfs_function_class (), 0.);
    gfs_function_set_units (function[i], vector[i]->units);
    GtsFile * fp1 = gts_file_new_from_string (component[i]);
    gfs_function_read (function[i], sim, fp1);
    if (fp1->type == GTS_ERROR) {
      gts_file_error (fp, "%s", fp1->error);
      gint j;
      for (j = 0; j <= i; j++) {
	gts_object_destroy (GTS_OBJECT (function[i]));
	function[i] = NULL;
      }	
      status = FALSE;
    }
    gts_file_destroy (fp1);
  }
  for (i = 0; i < FTT_DIMENSION; i++)
    g_free (component[i]);
  return status;
}
