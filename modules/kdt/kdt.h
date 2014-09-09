/* Gerris - The GNU Flow Solver
 * Copyright (C) 2010 National Institute of Water and Atmospheric Research
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

/* padding on 32 bits systems (to match automatic 64 bits padding) */

#if defined (__LP64__) || defined (__64BIT__) || defined (_LP64) || (__WORDSIZE == 64)
  #define PADDING_32_BITS
#else
  #define PADDING_32_BITS int padding;
#endif

FILE * kdt_tmpfile (void);

typedef struct {
  double x, y, z;
} KdtPoint;

typedef struct {
  float l, h;
} KdtInterval;

typedef KdtInterval KdtRect[2];

typedef struct {
  KdtPoint * p;
  long start, len, i, end, buflen, current;
  FILE * fp;
  void * buf;
} KdtHeap;

void kdt_heap_create (KdtHeap * h, FILE * fp, long start, long len, long buflen);
void kdt_heap_resize (KdtHeap * h, long len);
int  kdt_heap_get    (KdtHeap * h, KdtPoint * p);
void kdt_heap_split  (KdtHeap * h1, long len1, KdtHeap * h2);
void kdt_heap_put    (KdtHeap * h, KdtPoint * p);
void kdt_heap_flush  (KdtHeap * h);
void kdt_heap_free   (KdtHeap * h);
void kdt_heap_rewind (KdtHeap * h);

typedef struct { /* needs to be identical to RSurfaceSum in rsurface.h */
#if AVG_TERRAIN
  double H0;
  float Hmin, Hmax;
  int n;
#else
  double m01, m02, m03;
  double m11, m13;
  double m22, m23, m33;
  double m44, m55, m66, m77;
  double m67, m76;
  double H0, H1, H2, H3, H4;
  double H5, H6;
  float Hmin, Hmax;
  float coverage;
  int n;
#endif
} KdtSumCore;

typedef struct {
#if AVG_TERRAIN
  double H0;
  float Hmin, Hmax;
  int n;
#else
  double m01, m02, m03;
  double m11, m13;
  double m22, m23, m33;
  double m44, m55, m66, m77;
  double m67, m76;
  double H0, H1, H2, H3, H4;
  double H5, H6;
  float Hmin, Hmax;
  float coverage;
  int n;
#endif
  double w;
} KdtSum;

typedef struct _Kdt Kdt;
typedef int (* KdtCheck) (const KdtRect rect, void * data);

int kdt_intersects  (const KdtRect rect, const KdtRect query);
int kdt_includes    (const KdtRect rect, const KdtRect query);

Kdt * kdt_new       (void);
int  kdt_create     (Kdt * kdt, 
		     const char * name, 
		     int blksize,
		     KdtHeap * h,
		     void (* progress) (float complete, void * data),
		     void * data);
int  kdt_open       (Kdt * kdt, const char * name);
void kdt_destroy    (Kdt * kdt);
long kdt_query      (const Kdt * kdt, const KdtRect rect);
long kdt_query_sum  (const Kdt * kdt,
		     KdtCheck includes, KdtCheck intersects, void * data,
		     const KdtRect query, KdtSum * sum);
void kdt_sum_init   (KdtSum * s);
