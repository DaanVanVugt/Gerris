/* Gerris - The GNU Flow Solver
 * Copyright (C) 2010-2012 National Institute of Water and Atmospheric Research
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
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>

#include "kdt.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define VERSION 20120405 /* the file format version */

/* same as the system tmpfile() but uses the working directory (rather
   than /tmp) */
FILE * kdt_tmpfile (void)
{
  char name[] = "kdtXXXXXX";
  int fd = mkstemp (name);
  if (fd == -1) {
    perror ("kdt_tmpfile");
    exit (1);
  }
  FILE * fp = fdopen (fd, "r+w");
  assert (unlink (name) == 0);
  if (fp == NULL) {
    perror ("kdt_tmpfile");
    exit (1);
  }
  return fp;
}

/* refcounted buffer */

typedef struct {
  KdtPoint * p;
  int ref;
} Buffer;

static Buffer * buffer_new (long len)
{
  Buffer * b = malloc (sizeof (Buffer));
  b->p = malloc (len*sizeof (KdtPoint));
  b->ref = 1;
  return b;
}

static Buffer * buffer_ref (Buffer * b)
{
  b->ref++;
  return b;
}

static void buffer_unref (Buffer * b)
{
  b->ref--;
  if (b->ref == 0) {
    free (b->p);
    free (b);
  }
}

/* KdtHeap */

void kdt_heap_resize (KdtHeap * h, long len)
{
  assert (h->len < 0 || len < h->len);
  if (h->len == h->buflen) {
    h->buflen = len;
    h->end = h->buflen;
  }
  else if (len <= h->buflen) {
    h->buflen = len;
    kdt_heap_rewind (h);
    assert (h->end == len);
  }
  h->len = len;
}

static long heap_read (KdtHeap * h, long len)
{
  if (ftell (h->fp) != h->current)
    assert (fseek (h->fp, h->current, SEEK_SET) == 0);
  if (h->len > 0) {
    long maxlen = h->start + h->len - h->current/sizeof (KdtPoint);
    if (len > maxlen)
      len = maxlen;
  }
  long n = 0;
  if (len > 0) {
    n = fread (h->p, sizeof (KdtPoint), len, h->fp);
    h->current = ftell (h->fp);
  }
  return n;
}

static void heap_write (KdtHeap * h, long len)
{
  if (ftell (h->fp) != h->current)
    assert (fseek (h->fp, h->current, SEEK_SET) == 0);
  if (fwrite (h->p, sizeof (KdtPoint), len, h->fp) != len) {
    perror ("heap_write");
    exit (1);
  }
  h->current = ftell (h->fp);
}

void kdt_heap_create (KdtHeap * h, FILE * fp, long start, long len, long buflen)
{
  h->fp = fp;
  h->start = start;
  if (len > 0 && len < buflen)
    buflen = len;
  h->len = len;
  h->buflen = buflen;
  h->i = 0;
  h->buf = buffer_new (buflen);
  h->p = ((Buffer *) h->buf)->p;
  h->current = start*sizeof (KdtPoint);
  if (fp != NULL) {
    assert (fseek (fp, start*sizeof (KdtPoint), SEEK_SET) == 0);
    assert (ftell (fp) == h->current);
    h->end = heap_read (h, buflen);
    if (buflen == len)
      assert (h->end == len);
  }
  else
    h->end = 0;
}

void kdt_heap_rewind (KdtHeap * h)
{
  if (h->len == h->buflen) {
    h->i = 0;
    assert (h->end == h->buflen);
  }
  else {
    assert (fseek (h->fp, h->start*sizeof (KdtPoint), SEEK_SET) == 0);
    h->current = ftell (h->fp);
    h->end = heap_read (h, h->buflen);
    h->i = 0;
  }
}

int kdt_heap_get (KdtHeap * h, KdtPoint * p)
{
  if (h->len == h->buflen && h->i >= h->len)
    return 0;
  if (h->i < h->end) {
    *p = h->p[h->i++];
    return 1;
  }
  if (h->end < h->buflen)
    return 0;
  h->end = heap_read (h, h->buflen);
  h->i = 0;
  return kdt_heap_get (h, p);
}

void kdt_heap_split (KdtHeap * h1, long len1, KdtHeap * h2)
{
  assert (len1 < h1->len);
  if (h1->len == h1->buflen) {
    h2->fp = NULL;
    h2->start = 0;
    h2->len = h1->len - len1;
    h2->buflen = h2->len;
    h2->i = 0;
    h2->p = &h1->p[len1];
    h2->buf = buffer_ref (h1->buf);
    h2->end = h2->len;
    kdt_heap_resize (h1, len1);
  }
  else {
    kdt_heap_create (h2, h1->fp, h1->start + len1, h1->len - len1, h1->buflen);

    KdtHeap h;
    kdt_heap_create (&h, NULL, 0, len1, h1->buflen);
    if (len1 > h1->buflen)
      h.fp = kdt_tmpfile ();
    else
      h.end = h.len;

    kdt_heap_rewind (h1);
    long i;
    for (i = 0; i < len1; i++) {
      KdtPoint p;
      assert (kdt_heap_get (h1, &p));
      kdt_heap_put (&h, &p);
    }
    kdt_heap_flush (&h);
    h1->fp = NULL;
    kdt_heap_free (h1);
    *h1 = h;
  }
}

void kdt_heap_put (KdtHeap * h, KdtPoint * p)
{
  if (h->i == h->buflen) {
    heap_write (h, h->buflen);
    h->i = 0;
  }
  h->p[h->i++] = *p;
}

void kdt_heap_flush (KdtHeap * h)
{
  if (h->i > 0 && h->fp != NULL)
    heap_write (h, h->i);
}

void kdt_heap_free (KdtHeap * h)
{
  buffer_unref (h->buf);
  if (h->fp != NULL)
    assert (fclose (h->fp) == 0);
}

/* sort */

static int put (KdtHeap * h, KdtPoint * p, KdtHeap * merged)
{
  kdt_heap_put (merged, p);
  return kdt_heap_get (h, p);
}

static void kdt_write (KdtHeap * h, FILE * fp)
{
  kdt_heap_rewind (h);
  long i = 0;
  KdtPoint p;
  while (kdt_heap_get (h, &p)) {
    fprintf (fp, "%ld %g %g\n", i, p.x, p.y);
    i++;
  }
}

static void merge (KdtHeap * h1, KdtHeap * h2, 
		   int (*compar) (const void *, const void *),
		   long buflen)
{
  KdtHeap hm;
  assert (h1->len + h2->len > buflen);
  kdt_heap_create (&hm, NULL, h2->start - h1->len, h1->len + h2->len, buflen);
  hm.fp = h2->fp;
  KdtPoint p1, p2;
  kdt_heap_rewind (h1);
  int r1 = kdt_heap_get (h1, &p1);
  kdt_heap_rewind (h2);
  int r2 = kdt_heap_get (h2, &p2);
  while (r1 && r2) {
    if ((* compar) (&p2, &p1))
      r1 = put (h1, &p1, &hm);
    else
      r2 = put (h2, &p2, &hm);
  }
  while (r1)
    r1 = put (h1, &p1, &hm);
  while (r2)
    r2 = put (h2, &p2, &hm);
  kdt_heap_free (h1);
  h2->fp = NULL;
  kdt_heap_free (h2);
  kdt_heap_flush (&hm);
  *h1 = hm;
}

#if TIMING
static double elapsed (const struct timeval * start, const struct timeval * end)
{
  return (double) (end->tv_usec - start->tv_usec) + 
    1e6*(double) (end->tv_sec - start->tv_sec);
}
#endif

static void kdt_heap_sort (KdtHeap * h,
			   int  (*compar)   (const void *, const void *),
			   void (*progress) (void *), void * data)
{
#if TIMING
  struct timeval start;
  gettimeofday (&start, NULL);
#endif
  if (h->len == h->buflen) {
    qsort (h->p, h->len, sizeof (KdtPoint), compar);
    if (progress)
      (* progress) (data);
  }
  else {
    KdtHeap h2;
    long buflen = h->buflen;
    kdt_heap_split (h, h->len/2, &h2);
    kdt_heap_sort (h, compar, progress, data);
    kdt_heap_sort (&h2, compar, progress, data);
    merge (h, &h2, compar, buflen);
  }
#if TIMING
  struct timeval end;
  gettimeofday (&end, NULL);
  fprintf (stderr, "kdt_heap_sort %ld %g\n", len, elapsed (&start, &end));
#endif
}

/* number of qsort() calls for a given kdt_heap_sort() */
static int kdt_heap_sort_cost (long len, long buflen)
{
  int m = 1;
  while (len > buflen) {
    m *= 2;
    len /= 2;
  }
  return m;
}

/* Kdt */

typedef struct {
  KdtRect bound1, bound2;
  long len1;
  int n1;
} Node;

#define SYSTEM_32_BITS (!defined (__LP64__) && !defined (__64BIT__) && \
                        !defined (_LP64) && !(__WORDSIZE == 64))

typedef struct {
  KdtRect bound;
  long len;
  PADDING_32_BITS
  long np;
  int version;
} Header;

struct _Kdt {
  Header h;
  FILE * nodes, * sums, * leaves;
  KdtPoint * buffer;
  /* progress stuff */
  void (* progress) (float complete, void * data);
  void * data;
  int i, m;
};

static int check_32_bits (const Kdt * kdt)
{
#if SYSTEM_32_BITS
  long maxlen = (1 << 31)/sizeof (KdtPoint);
  if (kdt->h.len > maxlen) {
    fprintf (stderr, "kdt: 32-bits systems are limited to %ld data points\n", 
	     maxlen);
    return 1;
  }
#endif
  return 0;
}

static void sizes (const Node * n, long * nodes, long * sums, long * leaves)
{
  *nodes = n->n1*sizeof (Node);
  *sums  = n->n1*sizeof (KdtSumCore);
  *leaves = n->len1*sizeof (KdtPoint);
}

static void relative (const KdtRect rect, double * o, double * h)
{
  o[0] = (((double)rect[0].l) + ((double)rect[0].h))/2.;
  o[1] = (((double)rect[1].l) + ((double)rect[1].h))/2.;
  *h = ((double)rect[0].h) - ((double)rect[0].l);
  if (((double)rect[1].h) - ((double)rect[1].l) > *h)
    *h = ((double)rect[1].h) - ((double)rect[1].l);
}

static void sum_add_point (const KdtRect parent, KdtSumCore * sum, 
			   const KdtPoint * a, double w)
{
  double p[3], o[2], h;

  relative (parent, o, &h);

  p[0] = (a->x - o[0])/h; p[1] = (a->y - o[1])/h; p[2] = a->z;
#if AVG_TERRAIN
  sum->H0 += p[2];
  sum->n++;
  if (p[2] < sum->Hmin) sum->Hmin = p[2];
  if (p[2] > sum->Hmax) sum->Hmax = p[2];
#else
  sum->m01 += w*p[0];
  sum->m02 += w*p[1];
  sum->m03 += w*p[0]*p[1];
  sum->m11 += w*p[0]*p[0];
  sum->m13 += w*p[0]*p[0]*p[1];
  sum->m22 += w*p[1]*p[1];
  sum->m23 += w*p[0]*p[1]*p[1];
  sum->m33 += w*p[0]*p[0]*p[1]*p[1];
  sum->m44 += w*p[0]*p[0]*p[0];
  sum->m55 += w*p[1]*p[1]*p[1];
  sum->m66 += w*p[0]*p[0]*p[0]*p[0];
  sum->m77 += w*p[1]*p[1]*p[1]*p[1];
  sum->m67 += w*p[0]*p[0]*p[0]*p[1];
  sum->m76 += w*p[1]*p[1]*p[1]*p[0];
  sum->H0 += w*p[2];
  sum->H1 += w*p[0]*p[2];
  sum->H2 += w*p[1]*p[2];
  sum->H3 += w*p[0]*p[1]*p[2];
  sum->H4 += w*p[2]*p[2];
  sum->H5 += w*p[0]*p[0]*p[2];
  sum->H6 += w*p[1]*p[1]*p[2];
  sum->n ++;
  if (p[2] < sum->Hmin) sum->Hmin = p[2];
  if (p[2] > sum->Hmax) sum->Hmax = p[2];
#endif
}

static float area (const KdtRect rect)
{
  return (rect[0].h - rect[0].l)*(rect[1].h - rect[1].l);
}

static float length (const KdtRect rect)
{
  float w = rect[0].h - rect[0].l, h = rect[1].h - rect[1].l;
  return w > h ? w : h;
}

static void kdt_rect_write (const KdtRect rect, FILE * fp)
{
  fprintf (fp, "%f %f\n%f %f\n%f %f\n%f %f\n%f %f\n\n",
	   rect[0].l, rect[1].l, rect[0].h, rect[1].l,
	   rect[0].h, rect[1].h, rect[0].l, rect[1].h,
	   rect[0].l, rect[1].l);	   
}

static int sort_x (const void * p1, const void * p2)
{
  return ((KdtPoint *) p1)->x > ((KdtPoint *) p2)->x ? 1 : -1;
}

static int sort_y (const void * p1, const void * p2)
{
  return ((KdtPoint *) p1)->y > ((KdtPoint *) p2)->y ? 1 : -1;
}

static void kdt_sum_core_init (KdtSumCore * s)
{
  memset (s, 0, sizeof (KdtSumCore));
  s->Hmax = - 1e30;
  s->Hmin =   1e30;
}

#define GAP 0.2
#define MINLEN 6

static int update_sum (const KdtRect rect, KdtSumCore * n, KdtHeap * h, 
		       int index)
{
  kdt_sum_core_init (n);
  long i, imax = 0;
  double min, x, smax = 0;
  KdtPoint p;
  kdt_heap_rewind (h);
  assert (kdt_heap_get (h, &p));
  sum_add_point (rect, n, &p, 1.);
  min = x = (&p.x)[index];
  for (i = 1; i < h->len; i++) {
    assert (kdt_heap_get (h, &p));
    sum_add_point (rect, n, &p, 1.);
    double s = (&p.x)[index] - x;
    if (s > smax && i > MINLEN && i < h->len - MINLEN) {
      smax = s;
      imax = i;
    }
    x = (&p.x)[index];
  }
  return smax/(x - min) > GAP ? imax : -1;
}

static long update_bounds (KdtRect rect, KdtHeap * h)
{
  long len = 0;
  rect[0].h = rect[1].h = -1e30;
  rect[0].l = rect[1].l =  1e30;
  kdt_heap_rewind (h);
  KdtPoint p;
  while (kdt_heap_get (h, &p)) {
    if (p.x > rect[0].h) rect[0].h = p.x;
    if (p.x < rect[0].l) rect[0].l = p.x;
    if (p.y > rect[1].h) rect[1].h = p.y;
    if (p.y < rect[1].l) rect[1].l = p.y;
    len++;
  }
  return len;
}

static void progress (void * data)
{
  Kdt * kdt = data;
  if (kdt->progress && kdt->m > 0)
    (* kdt->progress) (++kdt->i/(float) kdt->m, kdt->data);
}

static void fwrite_check (const void * ptr, size_t size, size_t nmemb, 
			  FILE * stream)
{
  if (fwrite (ptr, size, nmemb, stream) != nmemb) {
    perror ("kdt_write");
    exit (1);
  }
}

static void union_bound (KdtRect b, const KdtRect b1, const KdtRect b2)
{
  b[0].l = MIN (b1[0].l, b2[0].l);
  b[1].l = MIN (b1[1].l, b2[1].l);
  b[0].h = MAX (b1[0].h, b2[0].h);
  b[1].h = MAX (b1[1].h, b2[1].h);
}

static int split (KdtHeap * h1, KdtRect bound, int index, Kdt * kdt, 
		  float * coverage)
{
#if TIMING
  struct timeval start;
  gettimeofday (&start, NULL);
#endif
  int ns = 0;
  if (h1->len > kdt->h.np) {
    //    fprintf (stderr, " splitting: %ld      \r", len);
    int nindex = (bound[0].h - bound[0].l < bound[1].h - bound[1].l);
    if (index != nindex) {
      kdt_heap_sort (h1, nindex ? sort_y : sort_x, progress, kdt);
      index = nindex;
    }
    else
      /* update cost estimate */
      kdt->m -= kdt_heap_sort_cost (h1->len, h1->buflen);
    KdtSumCore s;
    int imax = update_sum (bound, &s, h1, index);
    long spos = ftell (kdt->sums);
    fwrite_check (&s, sizeof (KdtSumCore), 1, kdt->sums);
    Node n;
    n.len1 = imax > 0 ? imax : h1->len/2;
#if TIMING
    struct timeval s1;
    gettimeofday (&s1, NULL);
#endif
    KdtHeap h2;
    kdt_heap_split (h1, n.len1, &h2);
#if TIMING
    struct timeval end;
    gettimeofday (&end, NULL);
    fprintf (stderr, "half %ld %f\n", len, elapsed (&s1, &end));
#endif
    update_bounds (n.bound1, h1);
    update_bounds (n.bound2, &h2);
    long pos = ftell (kdt->nodes);
    fwrite_check (&n, sizeof (Node), 1, kdt->nodes);
    float coverage1;
    n.n1 = split (h1, n.bound1, index, kdt, &coverage1);
    float coverage2;
    int n2 = split (&h2, n.bound2, index, kdt, &coverage2);
    ns = n.n1 + n2 + 1;

    /* update bound */
    union_bound (bound, n.bound1, n.bound2);

    /* update sums */
    double a = area (bound);
    if (a > 0.)
      s.coverage = (coverage1*area (n.bound1) + coverage2*area (n.bound2))/a;
    else
      s.coverage = 1.;
    assert (fseek (kdt->sums, spos + ((long) &s.coverage - (long) &s), 
		   SEEK_SET) == 0);
    fwrite_check (&s.coverage, sizeof (float), 1, kdt->sums);
    assert (fseek (kdt->sums, 0, SEEK_END) == 0);
    *coverage = s.coverage;

    /* update node */
    assert (fseek (kdt->nodes, pos, SEEK_SET) == 0);
    fwrite_check (&n, sizeof (Node), 1, kdt->nodes);
    assert (fseek (kdt->nodes, 0, SEEK_END) == 0);
  }
  else {
    assert (h1->len > 0);
    /* half the average distance between samples */
    double delta = length (bound)/sqrt (h1->len)/2.;
    bound[0].l -= delta;
    bound[1].l -= delta;
    bound[0].h += delta;
    bound[1].h += delta;
#if DEBUG
    kdt_rect_write (bound, stderr);
#endif
    assert (h1->len <= h1->buflen);
    fwrite_check (h1->p, sizeof (KdtPoint), h1->len, kdt->leaves);
    kdt_heap_free (h1);
    *coverage = 1.;
  }
#if TIMING
  struct timeval end;
  gettimeofday (&end, NULL);
  fprintf (stderr, "splitfile %ld %f\n", len, elapsed (&start, &end));
#endif
  return ns;
}

Kdt * kdt_new (void)
{
  Kdt * kdt = calloc (1, sizeof (Kdt));
  return kdt;
}

static FILE * open_ext (const char * name, const char * ext, const char * mode)
{
  int len = strlen (name), len1 = strlen (ext);
  char * fname = malloc (sizeof(char)*(len + len1 + 1));
  strcpy (fname, name);
  strcpy (&fname[len], ext);
  FILE * fp = fopen (fname, mode);
  free (fname);
  return fp;
}

static int kdt_init (Kdt * kdt, const char * name, int npmax, long len)
{
  kdt->nodes  = open_ext (name, ".kdt", "w");
  if (!kdt->nodes)
    return -1;
  
  kdt->sums   = open_ext (name, ".sum", "w");
  if (!kdt->sums)
    return -1;

  kdt->leaves = open_ext (name, ".pts", "w");
  if (!kdt->leaves)
    return -1;

  kdt->h.version = VERSION;
  kdt->h.len = len;
  kdt->h.np = npmax;
  kdt->h.bound[0].l = kdt->h.bound[1].l =  1e30;
  kdt->h.bound[0].h = kdt->h.bound[1].h = -1e30;

  if (check_32_bits (kdt))
    return -1;
  
  return 0;
}

int kdt_create (Kdt * kdt, const char * name, int blksize,
		KdtHeap * h,
		void (* progress) (float complete, void * data),
		void * data)
{
  KdtRect bound;
  long len = update_bounds (bound, h);
  kdt_heap_resize (h, len);

  int npmax = blksize/sizeof (KdtPoint);
  if (kdt_init (kdt, name, npmax, len))
    return -1;
  kdt->h.bound[0] = bound[0];
  kdt->h.bound[1] = bound[1];

  fwrite_check (&kdt->h, sizeof (Header), 1, kdt->nodes);

  /* cost estimate kdt->m (number of qsort() calls) based on balanced
     binary tree */
  kdt->m = kdt->i = 0;
  int m2 = 1;
  while (len > kdt->h.np) {
    kdt->m += kdt_heap_sort_cost (len, h->buflen)*m2;
    len /= 2;
    m2 *= 2;
  }
  kdt->progress = progress;
  kdt->data = data;
  float coverage;
  split (h, kdt->h.bound, -1, kdt, &coverage);
  /* write updated header (bounds have been updated) */
  rewind (kdt->nodes);
  fwrite_check (&kdt->h, sizeof (Header), 1, kdt->nodes);

  return 0;
}

int kdt_open (Kdt * kdt, const char * name)
{
  kdt->nodes  = open_ext (name, ".kdt", "r");
  if (!kdt->nodes)
    return -1;
  
  kdt->sums   = open_ext (name, ".sum", "r");
  if (!kdt->sums)
    return -1;

  kdt->leaves = open_ext (name, ".pts", "r");
  if (!kdt->leaves)
    return -1;

  if (fread (&kdt->h, sizeof (Header), 1, kdt->nodes) != 1)
    return -1;

  if (kdt->h.version != VERSION) {
    fprintf (stderr, 
	     "kdt: incompatible version number. Use:\n"
	     "%% kdt2kdt -v %s\n"
	     "to convert to the new format.\n", 
	     name);
    return -1;
  }

  kdt->buffer = malloc (sizeof (KdtPoint)*kdt->h.np);

  if (check_32_bits (kdt))
    return -1;

  return 0;
}

void kdt_destroy (Kdt * kdt)
{
  if (kdt->nodes)
    fclose (kdt->nodes);
  if (kdt->sums)
    fclose (kdt->sums);
  if (kdt->leaves)
    fclose (kdt->leaves);
  if (kdt->buffer)
    free (kdt->buffer);
  free (kdt);
}

int kdt_intersects (const KdtRect rect, const KdtRect query)
{
  return (rect[0].l <= query[0].h && rect[1].l <= query[1].h &&
	  rect[0].h >= query[0].l && rect[1].h >= query[1].l);
}

int kdt_includes (const KdtRect rect, const KdtRect query)
{
  return (rect[0].h <= query[0].h && rect[1].h <= query[1].h &&
	  rect[0].l >= query[0].l && rect[1].l >= query[1].l);
}

static long query (const Kdt * kdt, const KdtRect rect, long len)
{
  if (len > kdt->h.np) {
    Node node;
    if (fread (&node, sizeof (Node), 1, kdt->nodes) != 1)
      return -1;
    long pos = ftell (kdt->nodes), lpos = ftell (kdt->leaves);
    if (pos < 0 || lpos < 0)
      return -1;
    long n = 0;
    if (kdt_intersects (node.bound1, rect)) {
#if DEBUG
      kdt_rect_write (node.bound1, stderr);
#endif
      long n1 = query (kdt, rect, node.len1);
      if (n1 < 0)
	return -1;
      n += n1;
    }
    if (kdt_intersects (node.bound2, rect)) {
#if DEBUG
      kdt_rect_write (node.bound2, stderr);
#endif
      long snodes, ssums, sleaves;
      sizes (&node, &snodes, &ssums, &sleaves);
      if (fseek (kdt->nodes, pos + snodes, SEEK_SET))
	return -1;
      if (fseek (kdt->leaves, lpos + sleaves, SEEK_SET))
	return -1;
      long n1 = query (kdt, rect, len - node.len1);
      if (n1 < 0)
	return -1;
      n += n1;
    }
    return n;
  }
  else if (len > 0) {
    if (fread (kdt->buffer, sizeof (KdtPoint), len, kdt->leaves) != len)
      return -1;
    int i, n = 0;
    for (i = 0; i < len; i++)
      if (kdt->buffer[i].x >= rect[0].l && kdt->buffer[i].x <= rect[0].h && 
	  kdt->buffer[i].y >= rect[1].l && kdt->buffer[i].y <= rect[1].h) {
       	printf ("%.8f %.8f %f\n", 
		kdt->buffer[i].x, kdt->buffer[i].y, kdt->buffer[i].z);
	n++;
      }
    return n;
  }
  return 0;
}

long kdt_query (const Kdt * kdt, const KdtRect rect)
{
  rewind (kdt->nodes);
  rewind (kdt->leaves);
  Header h;
  if (fread (&h, sizeof (Header), 1, kdt->nodes) != 1)
    return -1;
  if (!kdt_intersects (rect, h.bound))
    return 0;
  return query (kdt, rect, h.len);
}

static void intersection (const KdtRect rect1, const KdtRect rect2, 
			  KdtRect inter)
{
  inter[0].l = MAX (rect1[0].l, rect2[0].l);
  inter[0].h = MIN (rect1[0].h, rect2[0].h);
  inter[1].l = MAX (rect1[1].l, rect2[1].l);
  inter[1].h = MIN (rect1[1].h, rect2[1].h);
  assert (inter[0].h >= inter[0].l && inter[1].h >= inter[1].l);  
}

static float intersection_area (const KdtRect rect1, const KdtRect rect2)
{
  KdtRect inter;
  intersection (rect1, rect2, inter);
  return area (inter);
}

static void sum_add_sum (const KdtRect parent, KdtSum * sum, 
			 const KdtRect rect, const KdtSumCore * a)
{
  /* weigh by effective area per sample */
  double w = intersection_area (rect, parent)*a->coverage/a->n;
  if (w == 0.)
    return;

  double op[2], oa[2], hp, ha;

  relative (parent, op, &hp);
  relative (rect, oa, &ha);

  double oap0 = oa[0] - op[0], oap1 = oa[1] - op[1];
  double an = a->n;
  double ha2 = ha*ha, hp2 = hp*hp;
  sum->m01 += w*(an*oap0 + a->m01*ha)/hp;
  sum->m02 += w*(an*oap1 + a->m02*ha)/hp;
  sum->m03 += w*(oap0*(an*oap1 + a->m02*ha) + ha*(a->m01*oap1 + a->m03*ha))/hp2;
  double m11 = (oap0*(an*oap0 + 2.*a->m01*ha) + a->m11*ha2)/hp2;
  sum->m11 += w*m11;
  double m13 = ha*(oap0*(a->m02*oap0 + 2.*a->m03*ha) + a->m13*ha2)/hp2;
  sum->m13 += w*(oap1*m11 + m13)/hp;
  double m22 = (oap1*(an*oap1 + 2.*a->m02*ha) + a->m22*ha2)/hp2;
  sum->m22 += w*m22;
  sum->m23 += w*(oap0*m22 + ha*(oap1*(oap1*a->m01 + 2.*a->m03*ha) + 
				a->m23*ha2)/hp2)/hp;
  sum->m33 += w*(oap1*(oap1*m11 + 2.*m13) + 
		 ha2*(oap0*(oap0*a->m22 + 2.*a->m23*ha) + ha2*a->m33)/hp2)/hp2;
  double ha3 = ha2*ha, hp3 = hp2*hp;
  sum->m44 += w*(oap0*(oap0*(oap0*an + 3.*ha*a->m01) + 3.*ha2*a->m11) + 
		 ha3*a->m44)/hp3;
  sum->m55 += w*(oap1*(oap1*(oap1*an + 3.*ha*a->m02) + 3.*ha2*a->m22) + 
		 ha3*a->m55)/hp3;
  double ha4 = ha3*ha, hp4 = hp3*hp;
  sum->m66 += w*(oap0*(oap0*(oap0*(oap0*an + 4.*ha*a->m01) + 6.*ha2*a->m11) 
		       + 4.*ha3*a->m44) + ha4*a->m66)/hp4;
  sum->m77 += w*(oap1*(oap1*(oap1*(oap1*an + 4.*ha*a->m02) + 6.*ha2*a->m22)
		       + 4.*ha3*a->m55) + ha4*a->m77)/hp4;
  sum->m67 += w*(oap1*(oap0*(oap0*(oap0*an + 3.*ha*a->m01) + 3.*ha2*a->m11) + 
		       ha3*a->m44)
		 + oap0*(oap0*(ha*a->m02*oap0 + 3.*ha2*a->m03) + 3.*ha3*a->m13) 
		 + ha4*a->m67)/hp4;
  sum->m76 += w*(oap0*(oap1*(oap1*(oap1*an + 3.*ha*a->m02) + 3.*ha2*a->m22) + 
		       ha3*a->m55)
		 + oap1*(oap1*(ha*a->m01*oap1 + 3.*ha2*a->m03) + 3.*ha3*a->m23)
		 + ha4*a->m76)/hp4;
  sum->H0 += w*a->H0;
  sum->H1 += w*(a->H0*oap0 + a->H1*ha)/hp;
  sum->H2 += w*(a->H0*oap1 + a->H2*ha)/hp;
  sum->H3 += w*(ha*(ha*a->H3 + oap0*a->H2 + oap1*a->H1) + oap0*oap1*a->H0)/hp2;
  sum->H4 += w*a->H4;
  sum->H5 += w*(oap0*(2.*ha*a->H1 + oap0*a->H0) + ha2*a->H5)/hp2;
  sum->H6 += w*(oap1*(2.*ha*a->H2 + oap1*a->H0) + ha2*a->H6)/hp2;
  sum->n += w*a->n*a->n/area (rect);
  float aparent = area (parent);
  if (aparent > 0.)
    sum->coverage += w*a->n/aparent;
  else
    sum->coverage = 1.;
  sum->w += w*a->n;
  if (a->Hmin < sum->Hmin) sum->Hmin = a->Hmin;
  if (a->Hmax > sum->Hmax) sum->Hmax = a->Hmax;
}

typedef struct {
  long np, sp, lp;
} FilePointers;

static long query_sum (const Kdt * kdt,
		       KdtCheck includes, KdtCheck intersects, void * data, 
		       KdtRect bound, long len,
		       FilePointers * f,
		       const KdtRect query, KdtSum * sum)
{
  if (len > kdt->h.np) {
    if (length (bound) <= length (query) || (* includes) (bound, data)) {
      KdtSumCore s;
      if (fseek (kdt->sums, f->sp, SEEK_SET))
	return -1;
      if (fread (&s, sizeof (KdtSumCore), 1, kdt->sums) != 1)
	return -1;
#if DEBUG
      fprintf (stderr, "read 1 sum %ld\n", sizeof (KdtSumCore));
#endif
      f->sp += sizeof (KdtSumCore);
      sum_add_sum (query, sum, bound, &s);
      return len;
    }
    f->sp += sizeof (KdtSumCore);

    Node node;
    if (fseek (kdt->nodes, f->np, SEEK_SET))
      return -1;
    if (fread (&node, sizeof (Node), 1, kdt->nodes) != 1)
      return -1;
    f->np += sizeof (Node);
#if DEBUG
    fprintf (stderr, "read 1 node %ld\n", sizeof (Node));
#endif

    long pos = f->np, lpos = f->lp, spos = f->sp;
    long n = 0;

    if ((* intersects) (node.bound1, data)) {
      long n1 = query_sum (kdt, includes, intersects, data, 
			   node.bound1, node.len1, f, query, sum);
      if (n1 < 0)
	return -1;
      n += n1;
    }

    if ((* intersects) (node.bound2, data)) {
      long snodes, ssums, sleaves;
      sizes (&node, &snodes, &ssums, &sleaves);
      f->np = pos + snodes;
      f->sp = spos + ssums;
      f->lp = lpos + sleaves;
      long n1 = query_sum (kdt, includes, intersects, data, 
			   node.bound2, len - node.len1, f, query, sum);
      if (n1 < 0)
	return -1;
      n += n1;
    }
    return n;
  }
  else {
    float h = length (bound)/sqrt (len); /* average distance between samples */
    if (h <= length (query)) {
#if DEBUG
      fprintf (stderr, "# area: %f %f\n", area (query), area (bound)/len);
      kdt_rect_write (bound, stderr);
#endif
      if (fseek (kdt->leaves, f->lp, SEEK_SET))
	return -1;
      if (fread (kdt->buffer, sizeof (KdtPoint), len, kdt->leaves) != len)
	return -1;
#if DEBUG
      fprintf (stderr, "read %ld leaves %ld\n", len, sizeof (KdtPoint));
#endif
      KdtPoint * a = kdt->buffer;
      int i, n = 0;
      for (i = 0; i < len; i++, a++) {
	KdtRect boundp;
	boundp[0].l = a->x - h/2.;
	boundp[0].h = a->x + h/2.;
	boundp[1].l = a->y - h/2.;
	boundp[1].h = a->y + h/2.;
	if ((* intersects) (boundp, data)) {
	  double w = intersection_area (boundp, query);
	  sum_add_point (query, (KdtSumCore *) sum, a, w);
	  sum->w += w;
	  sum->coverage += w/area (query);
	  n++;
	}
      }
      return n;
    }
  }
  return 0;
}

long kdt_query_sum (const Kdt * kdt, 
		    KdtCheck includes, KdtCheck intersects, void * data, 
		    const KdtRect query, KdtSum * sum)
{
  rewind (kdt->nodes);
  rewind (kdt->leaves);
  Header h;
  if (fread (&h, sizeof (Header), 1, kdt->nodes) != 1)
    return -1;
  FilePointers f;
  f.np = sizeof (Header);
  f.sp = f.lp = 0;
  if (!(* intersects) (h.bound, data))
    return 0;
  return query_sum (kdt, includes, intersects, data, h.bound, h.len, &f, 
		    query, sum);
}

void kdt_sum_init (KdtSum * s)
{
  kdt_sum_core_init ((KdtSumCore *) s);
  s->w = 0.;
}
