/* padding on 32 bits systems (to match automatic 64 bits padding) */

typedef struct _RSurface RSurface;

typedef struct { /* needs to be identical to typdirinfo in RStarTree.h */
  double m01, m02, m03;
  double m11, m13;
  double m22, m23, m33;
  double m44, m55, m66, m77;
  double m67, m76;
  double H0, H1, H2, H3, H4;
  double H5, H6;
  float Hmin, Hmax;
  int n;
  PADDING_32_BITS
} RSurfaceSum;

typedef struct {
  float l, h;
} RSurfaceInterval;

typedef RSurfaceInterval  RSurfaceRect[2];

typedef int (* RSurfaceCheck) (RSurfaceRect rect, void * data, int depth);

RSurface * r_surface_open     (const char * fname, const char * mode, int size);
int        r_surface_close    (RSurface * rt);
typedef int (* RSurfaceQuery) (double p[3], void * user_data);
void       r_surface_query_region (RSurface * rt, 
				   double min[2], double max[2],
				   RSurfaceQuery q, void * user_data);
