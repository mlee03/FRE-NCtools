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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "globals.h"
#include "mosaic_util.h"
#include "create_xgrid.h"
#include "create_xgrid_util.h"
#include "create_xgrid_acc.h"
#include "constant.h"
#include "openacc.h"

#define AREA_RATIO_THRESH (1.e-6)
#define MASK_THRESH       (0.5)
#define EPSLN8            (1.e-8)
#define EPSLN30           (1.0e-30)
#define EPSLN10           (1.0e-10)
#define MAX_V 8

/*******************************************************************************
  prepare_create_xgrid_2dx2d_order2_acc
*******************************************************************************/
int prepare_create_xgrid_2dx2d_acc(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out,
                                   const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                                   Minmaxavg_lists *out_minmaxavg, const double *mask_in,
                                   int *counts_per_ij1, int *ij2_start, int *ij2_end)
{

#define MAX_V 8
  int nx1p, nx2p, ny1p, ny2p, n1, n2;
  size_t approx_nxgrid;

  nx1p = nlon_in + 1;
  nx2p = nlon_out + 1;
  ny1p = nlat_in + 1;
  ny2p = nlat_out + 1;
  n1 = nlon_in*nlat_in;
  n2 = nlon_out*nlat_out;

  approx_nxgrid = 0;

#pragma acc data present(lon_out[0:nx2p*ny2p], lat_out[0:nx2p*ny2p], lon_in[0:nx1p*ny1p], lat_in[0:nx1p*ny1p], \
                         out_minmaxavg[0:1], out_minmaxavg->lon_min[0:n2], out_minmaxavg->lon_max[0:n2], \
                         out_minmaxavg->lat_min[0:n2], out_minmaxavg->lat_max[0:n2], out_minmaxavg->n[0:n2], \
                         out_minmaxavg->lon[0:n2],out_minmaxavg->lat[0:n2], out_minmaxavg->lon_avg[0:n2], \
                         counts_per_ij1[0:n1], ij2_start[0:n1], ij2_end[0:n1], mask_in[0:n1]) \
                   copy(approx_nxgrid)
#pragma acc parallel
{
#pragma acc loop independent reduction(+:approx_nxgrid)
  for( int ij1=0 ; ij1 < nlon_in*nlat_in ; ij1++) {

    int i1, j1;
    int icount=0, ij2_max=0 , ij2_min=nlon_out*nlat_out+1;

    i1 = ij1%nlon_in;
    j1 = ij1/nlon_in;

    counts_per_ij1[ij1]=0;

    if( mask_in[ij1] > MASK_THRESH ) {

      int n0, n1, n2, n3, n1_in;
      double lat_in_min, lat_in_max, lon_in_min, lon_in_max, lon_in_avg;
      double x1_in[MV], y1_in[MV];

      n0 = j1*nx1p+i1;       n1 = j1*nx1p+i1+1;
      n2 = (j1+1)*nx1p+i1+1; n3 = (j1+1)*nx1p+i1;
      x1_in[0] = lon_in[n0]; y1_in[0] = lat_in[n0];
      x1_in[1] = lon_in[n1]; y1_in[1] = lat_in[n1];
      x1_in[2] = lon_in[n2]; y1_in[2] = lat_in[n2];
      x1_in[3] = lon_in[n3]; y1_in[3] = lat_in[n3];
      lat_in_min = minval_double(4, y1_in);
      lat_in_max = maxval_double(4, y1_in);
      n1_in = fix_lon(x1_in, y1_in, 4, M_PI);
      lon_in_min = minval_double(n1_in, x1_in);
      lon_in_max = maxval_double(n1_in, x1_in);
      lon_in_avg = avgval_double(n1_in, x1_in);

#pragma acc loop independent reduction(+:approx_nxgrid) reduction(+:icount) reduction(min:ij2_min) reduction(max:ij2_max)
      for(int ij2=0; ij2<nlon_out*nlat_out; ij2++) {

        int i2, j2, l;
        double dx, lon_out_min, lon_out_max;
        double x2_in[MAX_V], y2_in[MAX_V],  x_out[MV], y_out[MV];;

        i2 = ij2%nlon_out;
        j2 = ij2/nlon_out;

        if(out_minmaxavg->lat_min[ij2] >= lat_in_max || out_minmaxavg->lat_max[ij2] <= lat_in_min ) continue;

        /* adjust x2_in according to lon_in_avg*/
        lon_out_min = out_minmaxavg->lon_min[ij2];
        lon_out_max = out_minmaxavg->lon_max[ij2];
        dx = out_minmaxavg->lon_avg[ij2] - lon_in_avg;

        if(dx < -M_PI ) {
          lon_out_min += TPI;
          lon_out_max += TPI;
        }

        else if (dx >  M_PI) {
          lon_out_min -= TPI;
          lon_out_max -= TPI;
        }

        /* x2_in should in the same range as x1_in after lon_fix, so no need to consider cyclic condition */
        if(lon_out_min >= lon_in_max || lon_out_max <= lon_in_min ) continue;


        //Note, the check for AREA_RATIO_THRESH has been removed
        //Thus, the computed value of approx_nxgrid will be equal to or greater than nxgrid
        approx_nxgrid++;
        icount++;
        ij2_min = min(ij2_min, ij2);
        ij2_max = max(ij2_max, ij2);

      } //ij2

      counts_per_ij1[ij1] = icount;
      ij2_start[ij1] = ij2_min ;
      ij2_end[ij1] = ij2_max;

    } // mask
  } //ij1
} //kernel

  return approx_nxgrid;

}

/*******************************************************************************
  create_xgrid_2dx2d_order1_acc
*******************************************************************************/
int create_xgrid_2dx2d_order1_acc(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out,
                                  const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                                  Minmaxavg_lists *out_minmaxavg, const double *mask_in, const int approx_nxgrid,
                                  const int *counts_per_ij1, const int *ij2_start, const int *ij2_end,
                                  Interp_config_mini *interp_mini, const int jstart, const int m)
{

#define MAX_V 8
  int nx1p, nx2p, ny1p, ny2p, n1, n2;
  double *area_in, *area_out;

  int *i_in, *j_in, *i_out, *j_out ;
  double *xgrid_area;

  int nxgrid;
  int *new_count;

  nx1p = nlon_in + 1;
  nx2p = nlon_out + 1;
  ny1p = nlat_in + 1;
  ny2p = nlat_out + 1;
  n1 = nlon_in*nlat_in;
  n2 = nlon_out*nlat_out;

  //Temporarily holds information about exchange grid cells because approx_nxgrid >= nxgrid
  i_in = (int *)malloc(approx_nxgrid*sizeof(int));
  j_in = (int *)malloc(approx_nxgrid*sizeof(int));
  i_out = (int *)malloc(approx_nxgrid*sizeof(int));
  j_out = (int *)malloc(approx_nxgrid*sizeof(int));
  xgrid_area = (double *)malloc(approx_nxgrid*sizeof(double));
  new_count=(int *)malloc(n1*sizeof(int));
  area_in  = (double *)malloc(n1*sizeof(double));
  area_out = (double *)malloc(n2*sizeof(double));

#pragma acc enter data create(i_in[0:approx_nxgrid], j_in[0:approx_nxgrid], \
                              i_out[0:approx_nxgrid], j_out[0:approx_nxgrid], \
                              xgrid_area[0:approx_nxgrid],new_count[0:n1], \
                              area_in[0:n1], area_out[0:n2])

  get_grid_area_acc(nlon_in, nlat_in, lon_in, lat_in, area_in);
  get_grid_area_acc(nlon_out, nlat_out, lon_out, lat_out, area_out);

  nxgrid = 0;

#pragma acc data present(lon_out[0:nx2p*ny2p], lat_out[0:nx2p*ny2p], lon_in[0:nx1p*ny1p], lat_in[0:nx1p*ny1p], \
                         out_minmaxavg[0:1], out_minmaxavg->lon_max[0:n2], out_minmaxavg->lon_min[0:n2], \
                         out_minmaxavg->lat_max[0:n2], out_minmaxavg->lat_min[0:n2], out_minmaxavg->n[0:n2], \
                         out_minmaxavg->lon[0:n2], out_minmaxavg->lat[0:n2], out_minmaxavg->lon_avg[0:n2], counts_per_ij1[0:n1], \
                         ij2_start[0:n1], ij2_end[0:n1], new_count[0:n1], mask_in[0:n1], area_in[0:n1], area_out[0:n2], \
                         xgrid_area[0:approx_nxgrid], j_in[0:approx_nxgrid], j_out[0:approx_nxgrid], \
                         i_out[0:approx_nxgrid], i_in[0:approx_nxgrid]) \
                  copyin(approx_nxgrid) \
                  copy(nxgrid)
#pragma acc parallel
{
#pragma acc loop independent reduction(+:nxgrid)
  for( int ij1=0 ; ij1<nlon_in*nlat_in ; ij1++) {
    if( mask_in[ij1] > MASK_THRESH ) {

      int n0, n1, n2, n3, n1_in;
      int i1, j1;
      int ij1_start, ixgrid;
      double lat_in_min, lat_in_max, lon_in_min, lon_in_max, lon_in_avg;
      double x1_in[MV], y1_in[MV];
      int ij1_jstart;

      i1 = ij1%nlon_in;
      j1 = ij1/nlon_in;

      n0 = j1*nx1p+i1;       n1 = j1*nx1p+i1+1;
      n2 = (j1+1)*nx1p+i1+1; n3 = (j1+1)*nx1p+i1;
      x1_in[0] = lon_in[n0]; y1_in[0] = lat_in[n0];
      x1_in[1] = lon_in[n1]; y1_in[1] = lat_in[n1];
      x1_in[2] = lon_in[n2]; y1_in[2] = lat_in[n2];
      x1_in[3] = lon_in[n3]; y1_in[3] = lat_in[n3];
      lat_in_min = minval_double(4, y1_in);
      lat_in_max = maxval_double(4, y1_in);

      n1_in = fix_lon(x1_in, y1_in, 4, M_PI);
      lon_in_min = minval_double(n1_in, x1_in);
      lon_in_max = maxval_double(n1_in, x1_in);
      lon_in_avg = avgval_double(n1_in, x1_in);

      ixgrid=0;
      // ij1_start, the total number of exchange grid cells computed for input cell ij1
      // is an approximation.
      ij1_start=0;
      if(ij1>0) {
#pragma acc loop seq
        for(int i=0 ; i<ij1 ; i++) ij1_start+=counts_per_ij1[i];
      }

#pragma acc loop seq reduction(+:nxgrid)
      for(int ij2=ij2_start[ij1]; ij2<=ij2_end[ij1]; ij2++) {

        int n_out, i2, j2, n2_in, l;
        double xarea, dx, lon_out_min, lon_out_max;
        double x2_in[MAX_V], y2_in[MAX_V],  x_out[MV], y_out[MV];;

        if(out_minmaxavg->lat_min[ij2] >= lat_in_max || out_minmaxavg->lat_max[ij2] <= lat_in_min ) continue;

        i2 = ij2%nlon_out;
        j2 = ij2/nlon_out;

        /* adjust x2_in according to lon_in_avg*/
        n2_in = out_minmaxavg->n[ij2];
#pragma acc loop seq
        for(l=0; l<n2_in; l++) {
          x2_in[l] = out_minmaxavg->lon[ij2*MAX_V+l];
          y2_in[l] = out_minmaxavg->lat[ij2*MAX_V+l];
        }
        lon_out_min = out_minmaxavg->lon_min[ij2];
        lon_out_max = out_minmaxavg->lon_max[ij2];
        dx = out_minmaxavg->lon_avg[ij2] - lon_in_avg;
        if(dx < -M_PI ) {
          lon_out_min += TPI;
          lon_out_max += TPI;
#pragma acc loop seq
          for (l=0; l<n2_in; l++) x2_in[l] += TPI;
        }
        else if (dx >  M_PI) {
          lon_out_min -= TPI;
          lon_out_max -= TPI;
#pragma acc loop seq
          for (l=0; l<n2_in; l++) x2_in[l] -= TPI;
        }

        /* x2_in should in the same range as x1_in after lon_fix, so no need to consider cyclic condition */
        if(lon_out_min >= lon_in_max || lon_out_max <= lon_in_min ) continue;

        n_out = 1;
        if (  (n_out = clip_2dx2d( x1_in, y1_in, n1_in, x2_in, y2_in, n2_in, x_out, y_out )) > 0) {
          double min_area;
          xarea = poly_area (x_out, y_out, n_out ) * mask_in[ij1];
          min_area = min(area_in[ij1], area_out[ij2]);
          if( xarea/min_area > AREA_RATIO_THRESH ) {
            xgrid_area[ij1_start+ixgrid] = xarea;
            i_in[ij1_start+ixgrid] = i1;
            j_in[ij1_start+ixgrid] = j1+jstart;
            i_out[ij1_start+ixgrid] = i2;
            j_out[ij1_start+ixgrid] = j2;
            ixgrid++;
            nxgrid++;
          } //if
        } //if
      } //ij2
      new_count[ij1]=ixgrid;
      ij1_jstart=(j1+jstart)*nlon_in+i1;
    } //mask
  } //ij1
} //kernel
//nxgrid is copied out

#pragma acc exit data delete(area_in, area_out)
 free(area_in); free(area_out);

 interp_mini->nxgrid = 0;
 if(nxgrid>0) {
   interp_mini->nxgrid = nxgrid;
   interp_mini->i_in = (int *)malloc(nxgrid*sizeof(int));
   interp_mini->j_in = (int *)malloc(nxgrid*sizeof(int));
   interp_mini->i_out = (int *)malloc(nxgrid*sizeof(int));
   interp_mini->j_out = (int *)malloc(nxgrid*sizeof(int));
   interp_mini->t_in = (int *)malloc(nxgrid*sizeof(int));
   interp_mini->area = (double *)malloc(nxgrid*sizeof(double));
#pragma acc enter data copyin(interp_mini[0:1])                         \
                       create(interp_mini->i_in[0:nxgrid], interp_mini->j_in[0:nxgrid],      \
                              interp_mini->i_out[0:nxgrid], interp_mini->j_out[0:nxgrid], interp_mini->t_in[0:nxgrid], \
                              interp_mini->area[0:nxgrid])

#pragma acc data present(counts_per_ij1[0:n1], new_count[0:n1], xgrid_area[0:nxgrid],\
                         i_in[0:nxgrid], j_in[0:nxgrid], j_out[0:nxgrid], i_out[0:nxgrid])
#pragma acc parallel loop
   for(int ij1=0 ; ij1<nlon_in*nlat_in ; ij1++){
     int ij1_start=0, ij1_start2=0;
#pragma acc loop
     for(int i=0 ; i<ij1 ; i++) {
       ij1_start+=counts_per_ij1[i];
       ij1_start2+=new_count[i];
     }
#pragma acc loop
     for(int i=0; i<new_count[ij1] ; i++){
       interp_mini->i_in[i+ij1_start2] = i_in[i+ij1_start];
       interp_mini->j_in[i+ij1_start2] = j_in[i+ij1_start];
       interp_mini->i_out[i+ij1_start2] = i_out[i+ij1_start];
       interp_mini->j_out[i+ij1_start2] = j_out[i+ij1_start];
       interp_mini->t_in[i+ij1_start2] = m;
       interp_mini->area[i+ij1_start2] = xgrid_area[i+ij1_start];
     }
   }
 }

#pragma acc exit data delete( i_in[0:approx_nxgrid], j_in[0:approx_nxgrid],\
                              i_out[0:approx_nxgrid], j_out[0:approx_nxgrid], \
                              xgrid_area[0:approx_nxgrid], new_count[0:n1])

 free(i_in); free(j_in); free(i_out); free(j_out);
 free(xgrid_area); free(new_count);

 return nxgrid;

};/* get_xgrid_2Dx2D_order1 */


/*******************************************************************************
  create_xgrid_2dx2d_order2 OPENACC version
*******************************************************************************/
int create_xgrid_2dx2d_order2_acc(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out,
                                  const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                                  Minmaxavg_lists *out_minmaxavg, const double *mask_in, const int approx_nxgrid,
                                  const int *counts_per_ij1, const int *ij2_start, const int *ij2_end,
                                  Interp_config_mini *interp_mini,  CellStruct *cell_in, const int jstart, const int m)
{

#define MAX_V 8
  int nx1p, nx2p, ny1p, ny2p, n1, n2;
  double *area_in, *area_out;

  int *i_in, *j_in, *i_out, *j_out ;
  double *xgrid_area, *xgrid_clon, *xgrid_clat;

  int nxgrid;
  int *new_count;

  nx1p = nlon_in + 1;
  nx2p = nlon_out + 1;
  ny1p = nlat_in + 1;
  ny2p = nlat_out + 1;
  n1 = nlon_in*nlat_in;
  n2 = nlon_out*nlat_out;

  //Temporarily holds information about exchange grid cells because approx_nxgrid >= nxgrid
  i_in = (int *)malloc(approx_nxgrid*sizeof(int));
  j_in = (int *)malloc(approx_nxgrid*sizeof(int));
  i_out = (int *)malloc(approx_nxgrid*sizeof(int));
  j_out = (int *)malloc(approx_nxgrid*sizeof(int));
  xgrid_area = (double *)malloc(approx_nxgrid*sizeof(double));
  xgrid_clon = (double *)malloc(approx_nxgrid*sizeof(double));
  xgrid_clat = (double *)malloc(approx_nxgrid*sizeof(double));
  new_count=(int *)malloc(n1*sizeof(int));
  area_in  = (double *)malloc(n1*sizeof(double));
  area_out = (double *)malloc(n2*sizeof(double));

#pragma acc enter data create(i_in[0:approx_nxgrid], j_in[0:approx_nxgrid], \
                              i_out[0:approx_nxgrid], j_out[0:approx_nxgrid], \
                              xgrid_area[0:approx_nxgrid], xgrid_clon[0:approx_nxgrid], \
                              xgrid_clat[0:approx_nxgrid], new_count[0:n1], area_in[0:n1], area_out[0:n2])

  get_grid_area_acc(nlon_in, nlat_in, lon_in, lat_in, area_in);
  get_grid_area_acc(nlon_out, nlat_out, lon_out, lat_out, area_out);

  nxgrid = 0;

#pragma acc data present(lon_out[0:nx2p*ny2p], lat_out[0:nx2p*ny2p], lon_in[0:nx1p*ny1p], lat_in[0:nx1p*ny1p],\
                         cell_in[0:1], out_minmaxavg[0:1],              \
                         cell_in->clat[0:n1], cell_in->clon[0:n1], cell_in->clat[0:n1], \
                         out_minmaxavg->lon_max[0:n2], out_minmaxavg->lon_min[0:n2], out_minmaxavg->lat_max[0:n2],\
                         out_minmaxavg->lat_min[0:n2], out_minmaxavg->n[0:n2], out_minmaxavg->lon[0:n2], \
                         out_minmaxavg->lat[0:n2], out_minmaxavg->lon_avg[0:n2], counts_per_ij1[0:n1],\
                         ij2_start[0:n1], ij2_end[0:n1], new_count[0:n1], mask_in[0:n1], area_in[0:n1], area_out[0:n2], \
                         xgrid_area[0:approx_nxgrid], xgrid_clon[0:approx_nxgrid], xgrid_clat[0:approx_nxgrid], \
                         j_in[0:approx_nxgrid], j_out[0:approx_nxgrid], i_out[0:approx_nxgrid], i_in[0:approx_nxgrid])\
            copyin(approx_nxgrid) \
            copy(nxgrid)
#pragma acc parallel
{
#pragma acc loop independent reduction(+:nxgrid)
  for( int ij1=0 ; ij1<nlon_in*nlat_in ; ij1++) {
    if( mask_in[ij1] > MASK_THRESH ) {

      int n0, n1, n2, n3, n1_in;
      int i1, j1;
      int ij1_start, ixgrid;
      double lat_in_min, lat_in_max, lon_in_min, lon_in_max, lon_in_avg;
      double x1_in[MV], y1_in[MV];
      double cell_in_area, cell_in_clon, cell_in_clat;
      int ij1_jstart;

      i1 = ij1%nlon_in;
      j1 = ij1/nlon_in;

      n0 = j1*nx1p+i1;       n1 = j1*nx1p+i1+1;
      n2 = (j1+1)*nx1p+i1+1; n3 = (j1+1)*nx1p+i1;
      x1_in[0] = lon_in[n0]; y1_in[0] = lat_in[n0];
      x1_in[1] = lon_in[n1]; y1_in[1] = lat_in[n1];
      x1_in[2] = lon_in[n2]; y1_in[2] = lat_in[n2];
      x1_in[3] = lon_in[n3]; y1_in[3] = lat_in[n3];
      lat_in_min = minval_double(4, y1_in);
      lat_in_max = maxval_double(4, y1_in);

      n1_in = fix_lon(x1_in, y1_in, 4, M_PI);
      lon_in_min = minval_double(n1_in, x1_in);
      lon_in_max = maxval_double(n1_in, x1_in);
      lon_in_avg = avgval_double(n1_in, x1_in);

      cell_in_area=0.;
      cell_in_clon=0.;
      cell_in_clat=0.;
      ixgrid=0;
      // ij1_start, the total number of exchange grid cells computed for input cell ij1
      // is an approximation.
      ij1_start=0;
      if(ij1>0) {
#pragma acc loop seq
        for(int i=0 ; i<ij1 ; i++) ij1_start+=counts_per_ij1[i];
      }

#pragma acc loop seq reduction(+:nxgrid) reduction(+:cell_in_area) reduction(+:cell_in_clon) reduction(+:cell_in_clat)
      for(int ij2=ij2_start[ij1]; ij2<=ij2_end[ij1]; ij2++) {

        int n_out, i2, j2, n2_in, l;
        double xarea, dx, lon_out_min, lon_out_max;
        double x2_in[MAX_V], y2_in[MAX_V],  x_out[MV], y_out[MV];;

        if(out_minmaxavg->lat_min[ij2] >= lat_in_max || out_minmaxavg->lat_max[ij2] <= lat_in_min ) continue;

        i2 = ij2%nlon_out;
        j2 = ij2/nlon_out;

        /* adjust x2_in according to lon_in_avg*/
        n2_in = out_minmaxavg->n[ij2];
#pragma acc loop seq
        for(l=0; l<n2_in; l++) {
          x2_in[l] = out_minmaxavg->lon[ij2*MAX_V+l];
          y2_in[l] = out_minmaxavg->lat[ij2*MAX_V+l];
        }
        lon_out_min = out_minmaxavg->lon_min[ij2];
        lon_out_max = out_minmaxavg->lon_max[ij2];
        dx = out_minmaxavg->lon_avg[ij2] - lon_in_avg;
        if(dx < -M_PI ) {
          lon_out_min += TPI;
          lon_out_max += TPI;
#pragma acc loop seq
          for (l=0; l<n2_in; l++) x2_in[l] += TPI;
        }
        else if (dx >  M_PI) {
          lon_out_min -= TPI;
          lon_out_max -= TPI;
#pragma acc loop seq
          for (l=0; l<n2_in; l++) x2_in[l] -= TPI;
        }

        /* x2_in should in the same range as x1_in after lon_fix, so no need to consider cyclic condition */
        if(lon_out_min >= lon_in_max || lon_out_max <= lon_in_min ) continue;

        n_out = 1;
        if (  (n_out = clip_2dx2d( x1_in, y1_in, n1_in, x2_in, y2_in, n2_in, x_out, y_out )) > 0) {
          double min_area;
          xarea = poly_area (x_out, y_out, n_out ) * mask_in[ij1];
          min_area = min(area_in[ij1], area_out[ij2]);
          if( xarea/min_area > AREA_RATIO_THRESH ) {
            xgrid_area[ij1_start+ixgrid] = xarea;
            xgrid_clon[ij1_start+ixgrid] = poly_ctrlon(x_out, y_out, n_out, lon_in_avg);
            xgrid_clat[ij1_start+ixgrid] = poly_ctrlat (x_out, y_out, n_out );
            cell_in_area += xarea;
            cell_in_clon += xgrid_clon[ij1_start+ixgrid];
            cell_in_clat += xgrid_clat[ij1_start+ixgrid];
            i_in[ij1_start+ixgrid] = i1;
            j_in[ij1_start+ixgrid] = j1+jstart;
            i_out[ij1_start+ixgrid] = i2;
            j_out[ij1_start+ixgrid] = j2;
            ixgrid++;
            nxgrid++;
          } //if
        } //if
      } //ij2
      new_count[ij1]=ixgrid;
      ij1_jstart=(j1+jstart)*nlon_in+i1;
      cell_in->area[ij1_jstart] = cell_in_area;
      cell_in->clon[ij1_jstart] = cell_in_clon;
      cell_in->clat[ij1_jstart] = cell_in_clat;
    } //mask
  } //ij1
} //kernel
//nxgrid is copied out

#pragma acc exit data delete(area_in, area_out)
 free(area_in); free(area_out);

 interp_mini->nxgrid = 0;
 if(nxgrid>0) {
   interp_mini->nxgrid = nxgrid;
   interp_mini->i_in = (int *)malloc(nxgrid*sizeof(int));
   interp_mini->j_in = (int *)malloc(nxgrid*sizeof(int));
   interp_mini->i_out = (int *)malloc(nxgrid*sizeof(int));
   interp_mini->j_out = (int *)malloc(nxgrid*sizeof(int));
   interp_mini->t_in = (int *)malloc(nxgrid*sizeof(int));
   interp_mini->area = (double *)malloc(nxgrid*sizeof(double));
   interp_mini->di_in = (double *)malloc(nxgrid*sizeof(double));
   interp_mini->dj_in = (double *)malloc(nxgrid*sizeof(double));
#pragma acc enter data copyin(interp_mini[0:1])                           \
                       create(interp_mini->i_in[0:nxgrid], interp_mini->j_in[0:nxgrid],      \
                              interp_mini->i_out[0:nxgrid], interp_mini->j_out[0:nxgrid], interp_mini->t_in[0:nxgrid], \
                              interp_mini->area[0:nxgrid], interp_mini->di_in[0:nxgrid],interp_mini->dj_in[0:nxgrid])

#pragma acc data present(counts_per_ij1[0:n1], new_count[0:n1], xgrid_area[0:nxgrid],\
                         xgrid_clon[0:nxgrid], xgrid_clat[0:nxgrid], i_in[0:nxgrid], \
                         j_in[0:nxgrid], j_out[0:nxgrid], i_out[0:nxgrid])
#pragma acc parallel loop
   for(int ij1=0 ; ij1<nlon_in*nlat_in ; ij1++){
     int ij1_start=0, ij1_start2=0;
#pragma acc loop
     for(int i=0 ; i<ij1 ; i++) {
       ij1_start+=counts_per_ij1[i];
       ij1_start2+=new_count[i];
     }
#pragma acc loop
     for(int i=0; i<new_count[ij1] ; i++){
       interp_mini->i_in[i+ij1_start2] = i_in[i+ij1_start];
       interp_mini->j_in[i+ij1_start2] = j_in[i+ij1_start];
       interp_mini->i_out[i+ij1_start2] = i_out[i+ij1_start];
       interp_mini->j_out[i+ij1_start2] = j_out[i+ij1_start];
       interp_mini->t_in[i+ij1_start2] = m;
       interp_mini->area[i+ij1_start2] = xgrid_area[i+ij1_start];
       interp_mini->di_in[i+ij1_start2] = xgrid_clon[i+ij1_start]/xgrid_area[i+ij1_start];
       interp_mini->dj_in[i+ij1_start2] = xgrid_clat[i+ij1_start]/xgrid_area[i+ij1_start];
     }
   }
 }

#pragma acc exit data delete( i_in[0:approx_nxgrid], j_in[0:approx_nxgrid],\
                              i_out[0:approx_nxgrid], j_out[0:approx_nxgrid], \
                              xgrid_area[0:approx_nxgrid], xgrid_clon[0:approx_nxgrid], \
                              xgrid_clat[0:approx_nxgrid], new_count[0:n1])

free(i_in); free(j_in); free(i_out); free(j_out);
free(xgrid_area); free(xgrid_clon); free(xgrid_clat);
free(new_count);

 return nxgrid;

};/* get_xgrid_2Dx2D_order2 */
