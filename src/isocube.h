/* Gerris - The GNU Flow Solver
 * Copyright (C) 2001-2004 National Institute of Water and Atmospheric
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

/* isocube adapted from GTS (see gts/src/iso.c and gts/doc/isocube.fig) */
static guint edge1[12][2] = {
  {0, 4}, {1, 5}, {3, 7}, {2, 6},
  {0, 2}, {1, 3}, {5, 7}, {4, 6},
  {0, 1}, {4, 5}, {6, 7}, {2, 3}
};
static FttVector vertex[8] = {
  {0.,0.,0.},{0.,0.,1.},{0.,1.,0.},{0.,1.,1.},
  {1.,0.,0.},{1.,0.,1.},{1.,1.,0.},{1.,1.,1.}
};
static guint face[6][4][2] = {
  {{7,0},{10,0},{6,1},{9,1}}, /* right */
  {{4,0},{11,0},{5,1},{8,1}}, /* left */
  {{3,0},{10,0},{2,1},{11,1}},/* top */
  {{0,0},{9,0},{1,1},{8,1}},  /* bottom */
  {{1,0},{6,0},{2,1},{5,1}},  /* front */
  {{0,0},{7,0},{3,1},{4,1}}   /* back */
};
static guint face_v[6][4] = {
  {4,6,7,5},/* right */
  {0,2,3,1},/* left */
  {2,6,7,3},/* top */
  {0,4,5,1},/* bottom */
  {1,5,7,3},/* front */
  {0,4,6,2} /* back */
};
static FttVector edge[12][2] = {
  {{0.,0.,0.},{1.,0.,0.}},{{0.,0.,1.},{1.,0.,1.}},{{0.,1.,1.},{1.,1.,1.}},{{0.,1.,0.},{1.,1.,0.}},
  {{0.,0.,0.},{0.,1.,0.}},{{0.,0.,1.},{0.,1.,1.}},{{1.,0.,1.},{1.,1.,1.}},{{1.,0.,0.},{1.,1.,0.}},
  {{0.,0.,0.},{0.,0.,1.}},{{1.,0.,0.},{1.,0.,1.}},{{1.,1.,0.},{1.,1.,1.}},{{0.,1.,0.},{0.,1.,1.}}
};
static guint connectv[12][2][4] = {
  {{4, 5, 1, 0}, {0, 2, 6, 4}}, /* 0 */
  {{5, 7, 3, 1}, {1, 0, 4, 5}}, /* 1 */
  {{7, 6, 2, 3}, {3, 1, 5, 7}}, /* 2 */
  {{6, 4, 0, 2}, {2, 3, 7, 6}}, /* 3 */
  {{2, 6, 4, 0}, {0, 1, 3, 2}}, /* 4 */
  {{3, 2, 0, 1}, {1, 5, 7, 3}}, /* 5 */
  {{7, 3, 1, 5}, {5, 4, 6, 7}}, /* 6 */
  {{6, 7, 5, 4}, {4, 0, 2, 6}}, /* 7 */
  {{1, 3, 2, 0}, {0, 4, 5, 1}}, /* 8 */
  {{5, 1, 0, 4}, {4, 6, 7, 5}}, /* 9 */
  {{7, 5, 4, 6}, {6, 2, 3, 7}}, /* 10 */
  {{3, 7, 6, 2}, {2, 0, 1, 3}}  /* 11 */
};
static FttVector cvertex[8] = {
  {0., 0., 0.}, {0., 0., 1.}, {0., 1., 0.}, {0., 1., 1.},
  {1., 0., 0.}, {1., 0., 1.}, {1., 1., 0.}, {1., 1., 1.}
};

#if FTT_2D
# define NCORNERS 4
static FttDirection corner[4][2] = {
  {FTT_LEFT, FTT_BOTTOM}, {FTT_RIGHT, FTT_BOTTOM}, {FTT_RIGHT, FTT_TOP}, {FTT_LEFT, FTT_TOP}
};
#else /* 3D */
# define NCORNERS 8
/* first index is the edge number, second index is the edge orientation 
   (0 or 1), third index are the edges which this edge may connect to
   in order and the corresponding face direction */
static guint connect[12][2][4] = {
  {{9, 1, 8, FTT_BOTTOM}, {4, 3, 7, FTT_BACK}},   /* 0 */
  {{6, 2, 5, FTT_FRONT},  {8, 0, 9, FTT_BOTTOM}}, /* 1 */
  {{10, 3, 11, FTT_TOP},  {5, 1, 6, FTT_FRONT}},  /* 2 */
  {{7, 0, 4, FTT_BACK},   {11, 2, 10, FTT_TOP}},  /* 3 */
  {{3, 7, 0, FTT_BACK},   {8, 5, 11, FTT_LEFT}},  /* 4 */
  {{11, 4, 8, FTT_LEFT},  {1, 6, 2, FTT_FRONT}},  /* 5 */
  {{2, 5, 1, FTT_FRONT},  {9, 7, 10, FTT_RIGHT}}, /* 6 */
  {{10, 6, 9, FTT_RIGHT}, {0, 4, 3, FTT_BACK}},   /* 7 */
  {{5, 11, 4, FTT_LEFT},  {0, 9, 1, FTT_BOTTOM}}, /* 8 */
  {{1, 8, 0, FTT_BOTTOM}, {7, 10, 6, FTT_RIGHT}}, /* 9 */
  {{6, 9, 7, FTT_RIGHT},  {3, 11, 2, FTT_TOP}},   /* 10 */
  {{2, 10, 3, FTT_TOP},   {4, 8, 5, FTT_LEFT}}    /* 11 */
};
static FttDirection corner[8][3] = {
  {FTT_LEFT, FTT_BOTTOM, FTT_BACK},
  {FTT_LEFT, FTT_BOTTOM, FTT_FRONT},
  {FTT_LEFT, FTT_TOP, FTT_BACK},
  {FTT_LEFT, FTT_TOP, FTT_FRONT},
  {FTT_RIGHT, FTT_BOTTOM, FTT_BACK},
  {FTT_RIGHT, FTT_BOTTOM, FTT_FRONT},
  {FTT_RIGHT, FTT_TOP, FTT_BACK},
  {FTT_RIGHT, FTT_TOP, FTT_FRONT}
};
#endif /* 3D */

#define SLIGHTLY_LARGER 1.001
