/***********************************************************************
 *                   GNU Lesser General Public License
 *
 * This file is part of the GFDL FRE NetCDF tools package (FRE-NCTools).
 *
 * FRE-NCtools is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * FRE-NCtools is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FRE-NCTools.  If not, see
 * <http://www.gnu.org/licenses/>.
 **********************************************************************/

#ifndef _GENERAL_UTILS_H
#define _GENERAL_UTILS_H 1
#endif

#define min(a,b) (a<b ? a:b)
#define max(a,b) (a>b ? a:b)

struct Node_acc{
  double x, y, z, u, u_clip;
  int intersect; /* indicate if this point is an intersection, 0 = no, 1= yes, 2=both intersect and vertices */
  int inbound;      /* -1 uninitialized, 0 coincident, 1 outbound, 2 inbound */
  int initialized; /* = 0 means empty list */
  int isInside;   /* = 1 means one point is inside the other polygon, 0 is not, -1 undecided. */
  int subj_index; /* the index of subject point that an intersection follow. */
  int clip_index; /* the index of clip point that an intersection follow */
  struct Node_acc *Next_acc;
};

#pragma omp declare target
int nearest_index_acc(double value, const double *array, int ia);
#pragma omp end declare target

#pragma omp declare target
int lon_fix_acc(double *x, double *y, int n_in, double tlon);
#pragma omp end declare target

#pragma omp declare target
double minval_double_acc(int size, const double *data);
#pragma omp end declare target

#pragma omp declare target
double maxval_double_acc(int size, const double *data);
#pragma omp end declare target

#pragma omp declare target
double avgval_double_acc(int size, const double *data);
#pragma omp end declare target

#pragma omp declare target
void latlon2xyz_acc(int size, const double *lon, const double *lat, double *x, double *y, double *z);
#pragma omp end declare target

#pragma omp declare target
void xyz2latlon_acc(int size, const double *x, const double *y, const double *z, double *lon, double *lat);
#pragma omp end declare target

#pragma omp declare target
double poly_area_main_acc(const double x[], const double y[], int n);
#pragma omp end declare target

#pragma omp declare target
double poly_area_acc(const double lon[], const double lat[], int n);
#pragma omp end declare target

#pragma omp declare target
int delete_vtx_acc(double x[], double y[], int n, int n_del);
#pragma omp end declare target

#pragma omp declare target
int insert_vtx_acc(double x[], double y[], int n, int n_ins, double lon_in, double lat_in);
#pragma omp end declare target

#pragma omp declare target
int fix_lon_acc(double lon[], double lat[], int n, double tlon);
#pragma omp end declare target

#pragma omp declare target
double spherical_angle_acc(const double *v1, const double *v2, const double *v3);
#pragma omp end declare target

#pragma omp declare target
void vect_cross_acc(const double *p1, const double *p2, double *e );
#pragma omp end declare target

#pragma omp declare target
double dot_acc(const double *p1, const double *p2);
#pragma omp end declare target

#pragma omp declare target
double metric_acc(const double *p) ;
#pragma omp end declare target

#pragma omp declare target
int intersect_tri_with_line_acc(const double *plane, const double *l1, const double *l2, double *p, double *t);
#pragma omp end declare target

#pragma omp declare target
void mult_acc(double m[], double v[], double out_v[]);
#pragma omp end declare target

#pragma omp declare target
int invert_matrix_3x3_acc(double m[], double m_inv[]);
#pragma omp end declare target

#pragma omp declare target
double great_circle_area_acc(int n, const double *x, const double *y, const double *z);
#pragma omp end declare target

#pragma omp declare target
int insidePolygon_acc(struct Node_acc *node, struct Node_acc *list );
#pragma omp end declare target

#pragma omp declare target
int inside_a_polygon_acc( double *lon1, double *lat1, int *npts, double *lon2, double *lat2);
#pragma omp end declare target

#pragma omp declare target
void rewindList_acc(void);
#pragma omp end declare target

#pragma omp declare target
struct Node_acc *getNext_acc();
#pragma omp end declare target

#pragma omp declare target
void initNode_acc(struct Node_acc *node);
#pragma omp end declare target

#pragma omp declare target
void addEnd_acc(struct Node_acc *list, double x, double y, double z, int intersect, double u, int inbound, int inside);
#pragma omp end declare target

#pragma omp declare target
int addIntersect_acc(struct Node_acc *list, double x, double y, double z, int intersect, double u1, double u2,
                int inbound, int is1, int ie1, int is2, int ie2);
#pragma omp end declare target

#pragma omp declare target
void insertIntersect_acc(struct Node_acc *list, double x, double y, double z, double u1, double u2, int inbound,
                     double x2, double y2, double z2);
#pragma omp end declare target

#pragma omp declare target
int length_acc(struct Node_acc *list);
#pragma omp end declare target

#pragma omp declare target
int samePoint_acc(double x1, double y1, double z1, double x2, double y2, double z2);
#pragma omp end declare target

#pragma omp declare target
int sameNode_acc(struct Node_acc node1, struct Node_acc node2);
#pragma omp end declare target

#pragma omp declare target
void addNode_acc(struct Node_acc *list, struct Node_acc nodeIn);
#pragma omp end declare target

#pragma omp declare target
struct Node_acc *getNode_acc(struct Node_acc *list, struct Node_acc inNode_acc);
#pragma omp end declare target

#pragma omp declare target
struct Node_acc *getNextNode_acc(struct Node_acc *list);
#pragma omp end declare target

#pragma omp declare target
void copyNode_acc(struct Node_acc *node_out, struct Node_acc node_in);
#pragma omp end declare target

#pragma omp declare target
void printNode_acc(struct Node_acc *list, char *str);
#pragma omp end declare target

#pragma omp declare target
int intersectInList_acc(struct Node_acc *list, double x, double y, double z);
#pragma omp end declare target

#pragma omp declare target
void insertAfter_acc(struct Node_acc *list, double x, double y, double z, int intersect, double u, int inbound,
                 double x2, double y2, double z2);
#pragma omp end declare target

#pragma omp declare target
double gridArea_acc(struct Node_acc *grid);
#pragma omp end declare target

#pragma omp declare target
int isIntersect_acc(struct Node_acc node);
#pragma omp end declare target

#pragma omp declare target
int getInbound_acc( struct Node_acc node );
#pragma omp end declare target

#pragma omp declare target
struct Node_acc *getLast_acc(struct Node_acc *list);
#pragma omp end declare target

#pragma omp declare target
int getFirstInbound_acc( struct Node_acc *list, struct Node_acc *nodeOut);
#pragma omp end declare target

#pragma omp declare target
void getCoordinate_acc(struct Node_acc node, double *x, double *y, double *z);
#pragma omp end declare target

#pragma omp declare target
void getCoordinates_acc(struct Node_acc *node, double *p);
#pragma omp end declare target

#pragma omp declare target
void setCoordinate_acc(struct Node_acc *node, double x, double y, double z);
#pragma omp end declare target

#pragma omp declare target
void setInbound_acc(struct Node_acc *interList, struct Node_acc *list);
#pragma omp end declare target

#pragma omp declare target
int isInside_acc(struct Node_acc *node);
#pragma omp end declare target

#pragma omp declare target
void set_rotate_poly_true_acc(void);
#pragma omp end declare target

#pragma omp declare target
int is_near_pole_acc(const double y[], int n);
#pragma omp end declare target

#pragma omp declare target
int crosses_pole_acc(const double x[], int n);
#pragma omp end declare target

#pragma omp declare target
void rotate_point_acc( double rv[], double rmat [][3]);
#pragma omp end declare target

#pragma omp declare target
void rotate_poly_acc(const double x[], const double y[], const int n, double xr[], double yr[]);
#pragma omp end declare target

#pragma omp declare target
void set_the_rotation_matrix_acc();
#pragma omp end declare target

#pragma omp declare target
void pimod_acc(double x[],int nn);
#pragma omp end declare target

#pragma omp declare target
int inside_edge_acc(double x0, double y0, double x1, double y1, double x, double y);
#pragma omp end declare target

#pragma omp declare target
int line_intersect_2D_3D_acc(double *a1, double *a2, double *q1, double *q2, double *q3,
                         double *intersect, double *u_a, double *u_q, int *inbound);
#pragma omp end declare target
