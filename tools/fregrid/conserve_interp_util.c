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
#include <string.h>
#include <netcdf.h>
#include <math.h>
#include "openacc.h"
#include "constant.h"
#include "globals.h"
#include "mpp.h"
#include "mpp_io.h"
#include "read_mosaic.h"
#include "mosaic_util.h"
#include "create_xgrid_util.h"
#include "conserve_interp_util.h"

#define  AREA_RATIO (1.e-3)

/*******************************************************************************
  void read_remap_file( int ntiles_in, int ntiles_out, Grid_config *grid_out,
                        Interp_config *interp, unsigned int opcode)
  Reads in the weight/remap file if provided
*******************************************************************************/
void read_remap_file( int ntiles_in, int ntiles_out, Grid_config *grid_out,
                      Interp_config *interp, unsigned int opcode){

  int *i_in=NULL, *j_in=NULL, *i_out=NULL, *j_out=NULL;
  double *xgrid_area=NULL, *tmp_area=NULL, *xgrid_clon=NULL, *xgrid_clat=NULL;
  int n, i;
  int zero=0;
  size_t nxgrid;
  double garea;

  Interp_config *pinterp;
  Interp_config_mini *pinterp_mini;
  int isc, jsc, iec, jec;

  garea = 4*M_PI*RADIUS*RADIUS;

  for(n=0; n<ntiles_out; n++) {
    if( interp[n].file_exist ) { /* reading from file */
      int *t_in, *ind;
      int fid, vid;

      nxgrid     = read_mosaic_xgrid_size(interp[n].remap_file);
      malloc_xgrid_arrays(nxgrid, &i_in, &j_in, &i_out, &j_out, &xgrid_area, &xgrid_clon, &xgrid_clat);
      t_in       = (int    *)malloc(nxgrid*sizeof(int   ));
      ind        = (int    *)malloc(nxgrid*sizeof(int   ));
      if(opcode & CONSERVE_ORDER1)
        read_mosaic_xgrid_order1(interp[n].remap_file, i_in, j_in, i_out, j_out, xgrid_area);
      else
        read_mosaic_xgrid_order2(interp[n].remap_file, i_in, j_in, i_out, j_out, xgrid_area, xgrid_clon, xgrid_clat);

      /*--- rescale the xgrid area */
#pragma acc enter data copyin(xgrid_area[0:nxgrid])
#pragma acc parallel loop present(xgrid_area[0:nxgrid]) copyin(garea,nxgrid)
      for(i=0; i<nxgrid; i++) xgrid_area[i] *= garea;

      fid = mpp_open(interp[n].remap_file, MPP_READ);
      vid = mpp_get_varid(fid, "tile1");
      mpp_get_var_value(fid, vid, t_in);
      mpp_close(fid);
      /*distribute the exchange grid on each pe according to target grid index*/
      interp[n].nxgrid = 0;
      isc = grid_out[n].isc;
      jsc = grid_out[n].jsc;
      jec = grid_out[n].jec;
      iec = grid_out[n].iec;
      for(i=0; i<nxgrid; i++) {
        if( i_out[i] <= iec && i_out[i] >= isc &&
            j_out[i] <= jec && j_out[i] >= jsc )
          ind[interp[n].nxgrid++] = i;
      }

#ifdef _OPENACC

      pinterp = interp+n;

      //OpenACC does not support array reduction yet
      for( int i=0 ; i<nxgrid ; i++) pinterp->intile[ t_in[i] ].nxgrid++;

      for(int m=0 ; m<ntiles_in ; m++) {
        nxgrid = pinterp->intile[m].nxgrid;
        pinterp->intile[m].i_in  = (int *)malloc(nxgrid*sizeof(int));
        pinterp->intile[m].j_in  = (int *)malloc(nxgrid*sizeof(int));
        pinterp->intile[m].i_out = (int *)malloc(nxgrid*sizeof(int));
        pinterp->intile[m].j_out = (int *)malloc(nxgrid*sizeof(int));
        pinterp->intile[m].area  = (double *)malloc(nxgrid*sizeof(double));
        pinterp->intile[m].di_in = (double *)malloc(nxgrid*sizeof(double));
        pinterp->intile[m].dj_in = (double *)malloc(nxgrid*sizeof(double));
      }

      for(int i=0 ; i<nxgrid ; i++) {
        int ii, iii;
        ii = t_in[i];
        iii = ind[i];
        pinterp->intile[ii].t_in[i] = t_in[iii] - 1 ;
        pinterp->intile[ii].i_in[i] = i_in[iii];
        pinterp->intile[ii].j_in[i] = j_in[iii];
        pinterp->intile[ii].i_out[i] = i_out[iii] - isc;
        pinterp->intile[ii].j_out[i] = j_out[iii] - jsc;
        if( opcode & CONSERVE_ORDER2) {
          pinterp->intile[ii].di_in[i] = xgrid_clon[iii];
          pinterp->intile[ii].dj_in[i] = xgrid_clat[iii];
        }
      }

#pragma acc enter data copyin( pinterp )
      for( int m=0 ; m < ntiles_in ; m++ ) {
        nxgrid = pinterp->intile[m].nxgrid ;
#pragma acc enter data copyin( pinterp->intile[m] )
#pragma acc enter data copyin( pinterp->intile[m].i_in[0:nxgrid] )
#pragma acc enter data copyin( pinterp->intile[m].j_in[0:nxgrid] )
#pragma acc enter data copyin( pinterp->intile[m].i_out[0:nxgrid] )
#pragma acc enter data copyin( pinterp->intile[m].j_out[0:nxgrid] )
#pragma acc enter data copyin( pinterp->intile[m].area[0:nxgrid] )
#pragma acc enter data copyin( pinterp->intile[m].di_in[0:nxgrid] )
#pragma acc enter data copyin( pinterp->intile[m].dj_in[0:nxgrid] )
      }

#else
      interp[n].i_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
      interp[n].j_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
      interp[n].i_out  = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
      interp[n].j_out  = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
      interp[n].area   = (double *)malloc(interp[n].nxgrid*sizeof(double));
      interp[n].t_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
      if(opcode & CONSERVE_ORDER2) {
        interp[n].di_in   = (double *)malloc(interp[n].nxgrid*sizeof(double));
        interp[n].dj_in   = (double *)malloc(interp[n].nxgrid*sizeof(double));
      }

      for(i=0; i< nxgrid; i++) {
        interp[n].i_in [i] = i_in [ind[i]];
        interp[n].j_in [i] = j_in [ind[i]];
        interp[n].t_in [i] = t_in [ind[i]] - 1;
        interp[n].i_out[i] = i_out[ind[i]] - isc;
        interp[n].j_out[i] = j_out[ind[i]] - jsc;
        interp[n].area [i] = xgrid_area[ind[i]];
      }
      if(opcode & CONSERVE_ORDER2) {
        for(i=0; i< nxgrid; i++) {
          interp[n].di_in[i] = xgrid_clon[ind[i]];
          interp[n].dj_in[i] = xgrid_clat[ind[i]];
        }
      }
#endif

      free(t_in);
      free(ind);
      malloc_xgrid_arrays(zero, &i_in, &j_in, &i_out, &j_out, &xgrid_area, &xgrid_clon, &xgrid_clat);

    }//if read from file

  } // ntiles
  if(mpp_pe() == mpp_root_pe())printf("NOTE: Finish reading index and weight for conservative interpolation from file.\n");

}//end read_regrid_weights

/*******************************************************************************
  void malloc_xgrid_arrays( int nsize, int **i_in, int **j_in, int **i_out, int **j_out,
                            double **xgrid_area, double **xgrid_clon, double **xgrid_clat )
  allocates arrays that will hold exchange grid information
*******************************************************************************/
void malloc_xgrid_arrays( int nsize, int **i_in, int **j_in, int **i_out, int **j_out,
                          double **xgrid_area, double **xgrid_clon, double **xgrid_clat )
{

  // free if malloc-ed
  if(*i_in!=NULL) {
    free(*i_in);
    *i_in=NULL;
  }
  if(*j_in!=NULL) {
    free(*j_in);
    *j_in=NULL;
  }
  if(*i_out!=NULL) {
    free(*i_out);
    *i_out=NULL;
  }
  if(*j_out!=NULL) {
    free(*j_out);
    *j_out=NULL;
  }
  if(*xgrid_area!=NULL) {
    free(*xgrid_area);
    *xgrid_area=NULL;
  }
  if(*xgrid_clon!=NULL) {
    free(*xgrid_clon);
    *xgrid_clon=NULL;
  }
  if(*xgrid_clat!=NULL) {
    free(*xgrid_clat);
    *xgrid_clat=NULL;
  }

  if(nsize>0) {
    *i_in       = (int *) malloc(nsize * sizeof(int   ));
    *j_in       = (int *) malloc(nsize * sizeof(int   ));
    *i_out      = (int *) malloc(nsize * sizeof(int   ));
    *j_out      = (int *) malloc(nsize * sizeof(int   ));
    *xgrid_area = (double *) malloc(nsize * sizeof(double));
    *xgrid_clon = (double *) malloc(nsize * sizeof(double));
    *xgrid_clat = (double *) malloc(nsize * sizeof(double));
  }

}

/*******************************************************************************
  void get_CellStruct
  Gathers exchange grid information from all ranks and
  stores information in cell_in structure.
  Cell_in stores exchange grid information corresponding to each input parent cell
*******************************************************************************/
void get_CellStruct(const int tile_in, const int nx_in, const int nxgrid, int *i_in, int *j_in,
                    double *xgrid_area, double *xgrid_clon, double *xgrid_clat,
                    CellStruct *cell_in)
{

  int g_nxgrid;
  int *g_i_in, *g_j_in;
  double *g_area, *g_clon, *g_clat;

  g_nxgrid = nxgrid;
  mpp_sum_int(1, &g_nxgrid);
  if(g_nxgrid > 0) {
    g_i_in = (int    *)malloc(g_nxgrid*sizeof(int   ));
    g_j_in = (int    *)malloc(g_nxgrid*sizeof(int   ));
    g_area = (double *)malloc(g_nxgrid*sizeof(double));
    g_clon = (double *)malloc(g_nxgrid*sizeof(double));
    g_clat = (double *)malloc(g_nxgrid*sizeof(double));
    mpp_gather_field_int   (nxgrid, i_in,       g_i_in);
    mpp_gather_field_int   (nxgrid, j_in,       g_j_in);
    mpp_gather_field_double(nxgrid, xgrid_area, g_area);
    mpp_gather_field_double(nxgrid, xgrid_clon, g_clon);
    mpp_gather_field_double(nxgrid, xgrid_clat, g_clat);
    for(int i=0; i<g_nxgrid; i++) {
      int ii;
      ii = g_j_in[i]*nx_in+g_i_in[i];
      cell_in[tile_in].area[ii] += g_area[i];
      cell_in[tile_in].clon[ii] += g_clon[i];
      cell_in[tile_in].clat[ii] += g_clat[i];
    }
    free(g_i_in);
    free(g_j_in);
    free(g_area);
    free(g_clon);
    free(g_clat);
  }

}
/*******************************************************************************
  void get_interp
  stores exchange grid information to the interp structure
********************************************************************************/
void get_interp( const int opcode, const int nxgrid, Interp_config *interp, const int m, const int n,
                 const int *i_in, const int *j_in, const int *i_out, const int *j_out,
                 const double *xgrid_clon, const double *xgrid_clat, const double *xgrid_area )
{

  int nxgrid_prev;
  int i;

  if(nxgrid > 0) {
    nxgrid_prev = interp[n].nxgrid;
    interp[n].nxgrid += nxgrid;
    if(nxgrid_prev == 0 ) {
      interp[n].i_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
      interp[n].j_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
      interp[n].i_out  = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
      interp[n].j_out  = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
      interp[n].area   = (double *)malloc(interp[n].nxgrid*sizeof(double));
      interp[n].t_in   = (int    *)malloc(interp[n].nxgrid*sizeof(int   ));
      for(i=0; i<interp[n].nxgrid; i++) {
        interp[n].t_in [i] = m;
        interp[n].i_in [i] = i_in [i];
        interp[n].j_in [i] = j_in [i];
        interp[n].i_out[i] = i_out[i];
        interp[n].j_out[i] = j_out[i];
        interp[n].area[i]  = xgrid_area[i];
      }
      if(opcode & CONSERVE_ORDER2) {
        interp[n].di_in   = (double *)malloc(interp[n].nxgrid*sizeof(double));
        interp[n].dj_in   = (double *)malloc(interp[n].nxgrid*sizeof(double));
        for(i=0; i<interp[n].nxgrid; i++) {
          interp[n].di_in [i] = xgrid_clon[i]/xgrid_area[i];
          interp[n].dj_in [i] = xgrid_clat[i]/xgrid_area[i];
        }
      }
    }
    else {
      interp[n].i_in   = (int    *)realloc(interp[n].i_in,  interp[n].nxgrid*sizeof(int   ));
      interp[n].j_in   = (int    *)realloc(interp[n].j_in,  interp[n].nxgrid*sizeof(int   ));
      interp[n].i_out  = (int    *)realloc(interp[n].i_out, interp[n].nxgrid*sizeof(int   ));
      interp[n].j_out  = (int    *)realloc(interp[n].j_out, interp[n].nxgrid*sizeof(int   ));
      interp[n].area   = (double *)realloc(interp[n].area,  interp[n].nxgrid*sizeof(double));
      interp[n].t_in   = (int    *)realloc(interp[n].t_in,  interp[n].nxgrid*sizeof(int   ));
      for(i=0; i<nxgrid; i++) {
        interp[n].t_in [nxgrid_prev+i] = m;
        interp[n].i_in [nxgrid_prev+i] = i_in [i];
        interp[n].j_in [nxgrid_prev+i] = j_in [i];
        interp[n].i_out[nxgrid_prev+i] = i_out[i];
        interp[n].j_out[nxgrid_prev+i] = j_out[i];
        interp[n].area [nxgrid_prev+i] = xgrid_area[i];
      }
      if(opcode & CONSERVE_ORDER2) {
        interp[n].di_in   = (double *)realloc(interp[n].di_in, interp[n].nxgrid*sizeof(double));
        interp[n].dj_in   = (double *)realloc(interp[n].dj_in, interp[n].nxgrid*sizeof(double));
        for(i=0; i<nxgrid; i++) {
          interp[n].di_in [i+nxgrid_prev] = xgrid_clon[i]/xgrid_area[i];
          interp[n].dj_in [i+nxgrid_prev] = xgrid_clat[i]/xgrid_area[i];
        }
      }
    }
  }  /* if(nxgrid>0) */

}
/*******************************************************************************
  void get_jstart_jend
  get the starting and ending indices of the input grid that
  overlaps with the output grid
********************************************************************************/
void get_jstart_jend( const int nx_out, const int ny_out, const int nx_in, const int ny_in,
                      const double *lat_out, const double *lat_in,
                      int *jstart, int *jend, int *ny_now )
{

  double y_min, y_max, yy ;
  int i, j;

  y_min = minval_double((nx_out+1)*(ny_out+1), lat_out);
  y_max = maxval_double((nx_out+1)*(ny_out+1), lat_out);
  *jstart = ny_in; *jend = -1;
  for(j=0; j<=ny_in; j++) for(i=0; i<=nx_in; i++) {
      yy = lat_in[j*(nx_in+1)+i];
      if( yy > y_min ) {
        if(j < *jstart ) *jstart = j;
      }
      if( yy < y_max ) {
        if(j > *jend ) *jend = j;
      }

    }
  *jstart = max(0, *jstart-1);
  *jend   = min(ny_in-1, *jend+1);
  *ny_now = *jend-*jstart+1;

}
/*******************************************************************************
  void get_interp_dij
********************************************************************************/
void get_interp_dij(const int ntiles_in, const int ntiles_out,
                    const Grid_config *grid_in, CellStruct *cell_in, Interp_config *interp)
{

  int nx_in, ny_in;
  double x1_in[50], y1_in[50], lon_in_avg, clon, clat;
  int    i, ii, j, n, n0, n1, n2, n3, n1_in, tile;

  for(n=0 ; n<ntiles_in; n++) {

    /* calcualte cell area */
    nx_in = grid_in[n].nx;
    ny_in = grid_in[n].ny;
    for(j=0; j<ny_in; j++) {
      for(i=0; i<nx_in; i++) {
        ii = j*nx_in + i;
        if(cell_in[n].area[ii] > 0) {
          if( fabs(cell_in[n].area[ii]-grid_in[n].cell_area[ii])/grid_in[n].cell_area[ii] < AREA_RATIO ) {
            cell_in[n].clon[ii] /= cell_in[n].area[ii];
            cell_in[n].clat[ii] /= cell_in[n].area[ii];
          }
          else {
            n0 = j*(nx_in+1)+i;       n1 = j*(nx_in+1)+i+1;
            n2 = (j+1)*(nx_in+1)+i+1; n3 = (j+1)*(nx_in+1)+i;
            x1_in[0] = grid_in[n].lonc[n0]; y1_in[0] = grid_in[n].latc[n0];
            x1_in[1] = grid_in[n].lonc[n1]; y1_in[1] = grid_in[n].latc[n1];
            x1_in[2] = grid_in[n].lonc[n2]; y1_in[2] = grid_in[n].latc[n2];
            x1_in[3] = grid_in[n].lonc[n3]; y1_in[3] = grid_in[n].latc[n3];
            n1_in = fix_lon(x1_in, y1_in, 4, M_PI);
            lon_in_avg = avgval_double(n1_in, x1_in);
            clon = poly_ctrlon(x1_in, y1_in, n1_in, lon_in_avg);
            clat = poly_ctrlat (x1_in, y1_in, n1_in );
            cell_in[n].clon[ii] = clon/grid_in[n].cell_area[ii];
            cell_in[n].clat[ii] = clat/grid_in[n].cell_area[ii];
          }
        }
      }
    }
  }
  for(n=0; n<ntiles_out; n++) {
    for(i=0; i<interp[n].nxgrid; i++) {
      tile = interp[n].t_in[i];
      ii   = interp[n].j_in[i] * grid_in[tile].nx + interp[n].i_in[i];
      interp[n].di_in[i] -= cell_in[tile].clon[ii];
      interp[n].dj_in[i] -= cell_in[tile].clat[ii];
    }
  }

  /* free the memory */
  for(n=0; n<ntiles_in; n++) {
    free(cell_in[n].area);
    free(cell_in[n].clon);
    free(cell_in[n].clat);
  }
  free(cell_in);


}
/*******************************************************************************
  void get_interp_dij_acc
********************************************************************************/
void get_interp_dij_acc(const int m, const int nx_in, const int ny_in, const double *grid_area,
                        const double *latc, const double *lonc, CellStruct *cell_in, Interp_config_mini *interp_mini)
{

  int nxp, nyp, n, nxgrid;
  nxp = nx_in + 1;
  nyp = ny_in + 1;
  n = nx_in * ny_in;
  nxgrid = interp_mini->nxgrid;

#pragma acc data present(cell_in->clon[0:n], cell_in->clat[0:n], cell_in->area[0:n],\
                         grid_area[0:n], latc[0:nxp*nyp], lonc[0:nxp*nyp] )
#pragma acc parallel loop collapse(2)
  for(int j=0; j<ny_in; j++) {
    for(int i=0; i<nx_in; i++) {
      double x1_in[50], y1_in[50], lon_in_avg, clon, clat;
      int    n, n0, n1, n2, n3, n1_in, ii;

      ii = j*nx_in + i;
      if(cell_in->area[ii] > 0) {
        if( fabs(cell_in->area[ii]-grid_area[ii])/grid_area[ii] < AREA_RATIO ) {
          cell_in->clon[ii] /= cell_in->area[ii];
          cell_in->clat[ii] /= cell_in->area[ii];
        }
        else {
          n0 = j*(nx_in+1)+i;       n1 = j*(nx_in+1)+i+1;
          n2 = (j+1)*(nx_in+1)+i+1; n3 = (j+1)*(nx_in+1)+i;
          x1_in[0] = lonc[n0]; y1_in[0] = latc[n0];
          x1_in[1] = lonc[n1]; y1_in[1] = latc[n1];
          x1_in[2] = lonc[n2]; y1_in[2] = latc[n2];
          x1_in[3] = lonc[n3]; y1_in[3] = latc[n3];
          n1_in = fix_lon(x1_in, y1_in, 4, M_PI);
          lon_in_avg = avgval_double(n1_in, x1_in);
          clon = poly_ctrlon(x1_in, y1_in, n1_in, lon_in_avg);
          clat = poly_ctrlat (x1_in, y1_in, n1_in );
          cell_in->clon[ii] = clon/grid_area[ii];
          cell_in->clat[ii] = clat/grid_area[ii];
        }
      }
    }
  }

#pragma acc data present(cell_in->clon[0:n], cell_in->clat[0:n],        \
                         interp_mini->j_in[0:nxgrid], interp_mini->i_in[0:nxgrid], \
                         interp_mini->di_in[0:nxgrid], interp_mini->dj_in[0:nxgrid])
#pragma acc parallel loop
  for(int i=0 ; i< nxgrid ; i++) {
    int ii;
    ii = interp_mini->j_in[i] * nx_in + interp_mini->i_in[i];
    interp_mini->di_in[i] -= cell_in->clon[ii];
    interp_mini->dj_in[i] -= cell_in->clat[ii];
  }

}

/*******************************************************************************
  void write_remap_acc
********************************************************************************/
void write_remap_acc(const int ntiles_out, const int ntiles_in, Grid_config *grid_out, Interp_config *interp, unsigned int opcode)
{

  int n, i, ii, nxgrid;
  Interp_config *pinterp ;

  for(n=0; n<ntiles_out; n++) {

    nxgrid=interp[n].nxgrid;

    if(nxgrid > 0) {
      size_t start[4], nwrite[4];
      int    fid, dim_string, dim_ncells, dim_two, dims[4];
      int    id_xgrid_area, id_tile1_dist;
      int    id_tile1_cell, id_tile2_cell, id_tile1;
      int    *gdata_int, *ldata_int;
      double *gdata_dbl;

      fid = mpp_open( interp[n].remap_file, MPP_WRITE);
      dim_string = mpp_def_dim(fid, "string", STRING);
      dim_ncells = mpp_def_dim(fid, "ncells", nxgrid);
      dim_two    = mpp_def_dim(fid, "two", 2);
      dims[0] = dim_ncells; dims[1] = dim_two;
      id_tile1      = mpp_def_var(fid, "tile1",      NC_INT, 1, &dim_ncells, 1,
                                  "standard_name", "tile_number_in_mosaic1");
      id_tile1_cell = mpp_def_var(fid, "tile1_cell", NC_INT, 2, dims, 1,
                                  "standard_name", "parent_cell_indices_in_mosaic1");
      id_tile2_cell = mpp_def_var(fid, "tile2_cell", NC_INT, 2, dims, 1,
                                  "standard_name", "parent_cell_indices_in_mosaic2");
      id_xgrid_area = mpp_def_var(fid, "xgrid_area", NC_DOUBLE, 1, &dim_ncells, 2,
                                  "standard_name", "exchange_grid_area", "units", "m2");
      if(opcode & CONSERVE_ORDER2) id_tile1_dist = mpp_def_var(fid, "tile1_distance", NC_DOUBLE, 2, dims, 1,
                                                               "standard_name", "distance_from_parent1_cell_centroid");
      mpp_end_def(fid);
      for(i=0; i<4; i++) {
        start[i] = 0 ;
        nwrite[i] = 1;
      }
      nwrite[0] = nxgrid;

      gdata_int = (int *)malloc(nxgrid*sizeof(int));

      for(int m=0 ; m<ntiles_in ; m++) {
        int inxgrid ;
        inxgrid=interp[n].intile[m].nxgrid;
#pragma acc update self( interp[n].intile[m].t_in[0:inxgrid], \
                         interp[n].intile[m].i_in[0:inxgrid], \
                         interp[n].intile[m].j_in[0:inxgrid],   \
                         interp[n].intile[m].i_out[0:inxgrid],  \
                         interp[n].intile[m].j_out[0:inxgrid],  \
                         interp[n].intile[m].area[0:inxgrid])
        if( opcode & CONSERVE_ORDER2 ) {
#pragma acc update self(interp[n].intile[m].di_in[0:inxgrid], interp[n].intile[m].dj_in[0:inxgrid])
        }
      }

      ii = 0;
      for( int m=0 ; m<ntiles_in ; m++ ) {
        pinterp = interp[n].intile+m;
        for( int i=0 ; i<pinterp->nxgrid ; i++ )gdata_int[ii++] = pinterp->t_in[i] + 1 ;
      }
      mpp_put_var_value(fid, id_tile1, gdata_int);

      ii = 0;
      for( int m=0 ; m<ntiles_in ; m++ ) {
        pinterp = interp[n].intile+m;
        for( int i=0 ; i<pinterp->nxgrid ; i++ )gdata_int[ii++] = pinterp->i_in[i] + 1 ;
      }
      mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, gdata_int);

      ii = 0;
      for( int m=0 ; m<ntiles_in ; m++ ) {
        pinterp = interp[n].intile+m;
        for( int i=0 ; i<pinterp->nxgrid ; i++ )gdata_int[ii++] = pinterp->i_out[i] + grid_out[n].isc + 1;
      }
      mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, gdata_int);

      ii = 0;
      for( int m=0 ; m<ntiles_in ; m++ ) {
        pinterp = interp[n].intile+m;
        for( int i=0 ; i<pinterp->nxgrid ; i++ )gdata_int[ii++] = pinterp->j_in[i] + 1 ;
      }
      start[1] = 1;
      mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, gdata_int);

      ii = 0;
      for( int m=0 ; m<ntiles_in ; m++ ) {
        pinterp = interp[n].intile+m;
        for( int i=0 ; i<pinterp->nxgrid ; i++ )gdata_int[ii++] = pinterp->j_out[i] + grid_out[n].jsc + 1;
      }
      mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, gdata_int);

      free(gdata_int);
      gdata_dbl = (double *)malloc(nxgrid*sizeof(double));

      ii = 0;
      for( int m=0 ; m<ntiles_in ; m++ ) {
        pinterp = interp[n].intile+m;
        for( int i=0 ; i<pinterp->nxgrid ; i++ )gdata_dbl[ii++] = pinterp->area[i] ;
      }
      mpp_put_var_value(fid, id_xgrid_area, gdata_dbl);

      if(opcode & CONSERVE_ORDER2) {
        ii = 0 ;
        start[1] = 0;
        for( int m=0 ; m<ntiles_in ; m++ ) {
          pinterp = interp[n].intile+m;
          for( int i=0 ; i<pinterp->nxgrid ; i++ )gdata_dbl[ii++] = pinterp->di_in[i];
        }
        mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, gdata_dbl);

        ii = 0 ;
        start[1] = 1;
        for( int m=0 ; m<ntiles_in ; m++ ) {
          pinterp = interp[n].intile+m;
          for( int i=0 ; i<pinterp->nxgrid ; i++ )gdata_dbl[ii++] = pinterp->dj_in[i];
        }
        mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, gdata_dbl);
      }

      free(gdata_dbl);
      mpp_close(fid);
    }
  }

}

void write_remap(const int ntiles_out, Grid_config *grid_out, Interp_config *interp, unsigned int opcode)
{

  int n, i;

  for(n=0; n<ntiles_out; n++) {
    int nxgrid;

    nxgrid = interp[n].nxgrid;
    mpp_sum_int(1, &nxgrid);
    if(nxgrid > 0) {
      size_t start[4], nwrite[4];
      int    fid, dim_string, dim_ncells, dim_two, dims[4];
      int    id_xgrid_area, id_tile1_dist;
      int    id_tile1_cell, id_tile2_cell, id_tile1;
      int    *gdata_int, *ldata_int;
      double *gdata_dbl;

      fid = mpp_open( interp[n].remap_file, MPP_WRITE);
      dim_string = mpp_def_dim(fid, "string", STRING);
      dim_ncells = mpp_def_dim(fid, "ncells", nxgrid);
      dim_two    = mpp_def_dim(fid, "two", 2);
      dims[0] = dim_ncells; dims[1] = dim_two;
      id_tile1      = mpp_def_var(fid, "tile1",      NC_INT, 1, &dim_ncells, 1,
                                  "standard_name", "tile_number_in_mosaic1");
      id_tile1_cell = mpp_def_var(fid, "tile1_cell", NC_INT, 2, dims, 1,
                                  "standard_name", "parent_cell_indices_in_mosaic1");
      id_tile2_cell = mpp_def_var(fid, "tile2_cell", NC_INT, 2, dims, 1,
                                  "standard_name", "parent_cell_indices_in_mosaic2");
      id_xgrid_area = mpp_def_var(fid, "xgrid_area", NC_DOUBLE, 1, &dim_ncells, 2,
                                  "standard_name", "exchange_grid_area", "units", "m2");
      if(opcode & CONSERVE_ORDER2) id_tile1_dist = mpp_def_var(fid, "tile1_distance", NC_DOUBLE, 2, dims, 1,
                                                               "standard_name", "distance_from_parent1_cell_centroid");
      mpp_end_def(fid);
      for(i=0; i<4; i++) {
        start[i] = 0; nwrite[i] = 1;
      }
      nwrite[0] = nxgrid;
      gdata_int = (int *)malloc(nxgrid*sizeof(int));
      if(interp[n].nxgrid>0) ldata_int = (int *)malloc(interp[n].nxgrid*sizeof(int));
      mpp_gather_field_int(interp[n].nxgrid, interp[n].t_in, gdata_int);
      for(i=0; i<nxgrid; i++) gdata_int[i]++;
      mpp_put_var_value(fid, id_tile1, gdata_int);

      mpp_gather_field_int(interp[n].nxgrid, interp[n].i_in, gdata_int);
      for(i=0; i<nxgrid; i++) gdata_int[i]++;
      mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, gdata_int);

      for(i=0; i<interp[n].nxgrid; i++) ldata_int[i] = interp[n].i_out[i] + grid_out[n].isc + 1;
      mpp_gather_field_int(interp[n].nxgrid, ldata_int, gdata_int);
      mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, gdata_int);

      mpp_gather_field_int(interp[n].nxgrid, interp[n].j_in, gdata_int);
      for(i=0; i<nxgrid; i++) gdata_int[i]++;
      start[1] = 1;
      mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, gdata_int);

      for(i=0; i<interp[n].nxgrid; i++) ldata_int[i] = interp[n].j_out[i] + grid_out[n].jsc + 1;
      mpp_gather_field_int(interp[n].nxgrid, ldata_int, gdata_int);
      mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, gdata_int);

      free(gdata_int);
      if(interp[n].nxgrid>0)free(ldata_int);

      gdata_dbl = (double *)malloc(nxgrid*sizeof(double));
      mpp_gather_field_double(interp[n].nxgrid, interp[n].area, gdata_dbl);
      mpp_put_var_value(fid, id_xgrid_area, gdata_dbl);

      if(opcode & CONSERVE_ORDER2) {
        start[1] = 0;
        mpp_gather_field_double(interp[n].nxgrid, interp[n].di_in, gdata_dbl);
        mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, gdata_dbl);
        start[1] = 1;
        mpp_gather_field_double(interp[n].nxgrid, interp[n].dj_in, gdata_dbl);
        mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, gdata_dbl);
      }

      free(gdata_dbl);
      mpp_close(fid);
    }
  }

}

void check_conserve(const int ntiles_out, Grid_config *grid_out, Interp_config *interp)
{

  int nx1, ny1, max_i, max_j, i, j, ii, n;
  double max_ratio, ratio_change;
  double *area2;

  /* sum over exchange grid to get the area of grid_in */
  nx1  = grid_out[0].nxc;
  ny1  = grid_out[0].nyc;

  area2 = (double *)malloc(nx1*ny1*sizeof(double));

  for(n=0; n<ntiles_out; n++) {
    for(i=0; i<nx1*ny1; i++) area2[i] = 0;
    for(i=0; i<interp[n].nxgrid; i++) {
      ii = interp[n].j_out[i]*nx1 + interp[n].i_out[i];
      area2[ii] +=  interp[n].area[i];
    }
    max_ratio = 0;
    max_i = 0;
    max_j = 0;
    /* comparing area1 and area2 */
    for(j=0; j<ny1; j++) for(i=0; i<nx1; i++) {
        ii = j*nx1+i;
        ratio_change = fabs(grid_out[n].cell_area[ii]-area2[ii])/grid_out[n].cell_area[ii];
        if(ratio_change > max_ratio) {
          max_ratio = ratio_change;
          max_i = i;
          max_j = j;
        }
        if( ratio_change > 1.e-4 ) {
          printf("(i,j)=(%d,%d), change = %g, area1=%g, area2=%g\n", i, j, ratio_change, grid_out[n].cell_area[ii],area2[ii]);
        }
      }
    ii = max_j*nx1+max_i;
    printf("The maximum ratio change at (%d,%d) = %g, area1=%g, area2=%g\n",
           max_i, max_j, max_ratio, grid_out[n].cell_area[ii],area2[ii]);
  }
    free(area2);
}
