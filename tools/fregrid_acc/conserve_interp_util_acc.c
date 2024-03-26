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
#include <openacc.h>
#include "constant.h"
#include "globals.h"
#include "mpp.h"
#include "mpp_io.h"
#include "read_mosaic.h"
#include "mosaic_util.h"
#include "conserve_interp_util_acc.h"

/*******************************************************************************
  void read_remap_file
  Reads in the weight/remap file if provided
*******************************************************************************/
void read_remap_file(int ntiles_in, int ntiles_out, Grid_config *grid_out,
                     Interp_config *interp, unsigned int opcode)
{

  int *i_in=NULL, *j_in=NULL, *i_out=NULL, *j_out=NULL;
  double *xgrid_area=NULL, *xgrid_clon=NULL, *xgrid_clat=NULL;

  size_t nxgrid, nxgrid_acc;
  double garea;

  int isc, jsc, iec, jec;

  acc_copyin(interp, ntiles_out);

  for(int n=0; n<ntiles_out; n++) {
    if( interp[n].file_exist ) { //reading from file
      int *t_in, *ind, *ind_acc;
      int fid, vid, nxgrid;

      nxgrid = read_mosaic_xgrid_size(interp[n].remap_file);
      t_in  = (int *)malloc(nxgrid*sizeof(int));
      ind   = (int *)malloc(nxgrid*sizeof(int));
      i_in  = (int *)malloc(nxgrid*sizeof(int));
      j_in  = (int *)malloc(nxgrid*sizeof(int));
      i_out = (int *)malloc(nxgrid*sizeof(int));
      j_out = (int *)malloc(nxgrid*sizeof(int));
      xgrid_area = (double *)malloc(nxgrid*sizeof(double));
      if( opcode & CONSERVE_ORDER2) {
        xgrid_clon = (double *)malloc(nxgrid*sizeof(double));
        xgrid_clat = (double *)malloc(nxgrid*sizeof(double));
      }

      if(opcode & CONSERVE_ORDER1)
        read_mosaic_xgrid_order1(interp[n].remap_file, i_in, j_in, i_out, j_out, xgrid_area);
      else
        read_mosaic_xgrid_order2(interp[n].remap_file, i_in, j_in, i_out, j_out, xgrid_area, xgrid_clon, xgrid_clat);

      //rescale the xgrid area
      for(int i=0; i<nxgrid; i++) xgrid_area[i] *= garea;

      fid = mpp_open(interp[n].remap_file, MPP_READ);
      vid = mpp_get_varid(fid, "tile1");
      mpp_get_var_value(fid, vid, t_in);
      mpp_close(fid);

      //distribute the exchange grid on each pe according to target grid index
      interp[n].nxgrid = 0;
      isc = grid_out[n].isc;
      jsc = grid_out[n].jsc;
      jec = grid_out[n].jec;
      iec = grid_out[n].iec;
      for(int i=0; i<nxgrid; i++) {
        if( i_out[i] <= iec && i_out[i] >= isc &&
            j_out[i] <= jec && j_out[i] >= jsc )
          ind[interp[n].nxgrid++] = i;
      }

      for(int i=0 ; i<nxgrid ; i++) t_in[i]--;

      //OpenACC does not support array reduction yet
      for(int m=0 ; m<ntiles_in ; m++) interp[n].interp_mini[m].nxgrid=0;
      for(int i=0 ; i<nxgrid ; i++) interp[n].interp_mini[ t_in[i] ].nxgrid++;

      for(int m=0 ; m<ntiles_in ; m++) {
        nxgrid_acc = interp[n].interp_mini[m].nxgrid;
        interp[n].interp_mini[m].i_in  = (int *)malloc(nxgrid_acc*sizeof(int));
        interp[n].interp_mini[m].j_in  = (int *)malloc(nxgrid_acc*sizeof(int));
        interp[n].interp_mini[m].i_out = (int *)malloc(nxgrid_acc*sizeof(int));
        interp[n].interp_mini[m].j_out = (int *)malloc(nxgrid_acc*sizeof(int));
        interp[n].interp_mini[m].area  = (double *)malloc(nxgrid_acc*sizeof(double));
        if(opcode & CONSERVE_ORDER2) {
          interp[n].interp_mini[m].di_in = (double *)malloc(nxgrid_acc*sizeof(double));
          interp[n].interp_mini[m].dj_in = (double *)malloc(nxgrid_acc*sizeof(double));
        }
      }

      ind_acc = (int *)calloc(ntiles_in, sizeof(int));
      for(int i=0 ; i<nxgrid ; i++) {
        int ii, iii, iv;
        ii = t_in[i];
        iii = ind[i];
        iv=ind_acc[ii];
        interp[n].interp_mini[ii].i_in[iv] = i_in[iii];
        interp[n].interp_mini[ii].j_in[iv] = j_in[iii];
        interp[n].interp_mini[ii].i_out[iv] = i_out[iii] - isc;
        interp[n].interp_mini[ii].j_out[iv] = j_out[iii] - jsc;
        interp[n].interp_mini[ii].area[iv] = xgrid_area[iii];
        if( opcode & CONSERVE_ORDER2) {
          interp[n].interp_mini[ii].di_in[iv] = xgrid_clon[iii];
          interp[n].interp_mini[ii].dj_in[iv] = xgrid_clat[iii];
        }
        ind_acc[ii]++;
      }

      acc_copyin(interp[n].interp_mini, ntiles_in);
      for( int m=0 ; m < ntiles_in ; m++ ) {
        nxgrid_acc = interp[n].interp_mini[m].nxgrid ;
        acc_copyin(interp[n].interp_mini[m].i_in, nxgrid_acc);
        acc_copyin(interp[n].interp_mini[m].j_in, nxgrid_acc);
        acc_copyin(interp[n].interp_mini[m].i_out, nxgrid_acc);
        acc_copyin(interp[n].interp_mini[m].j_out, nxgrid_acc);
        if( opcode & CONSERVE_ORDER2) {
          acc_copyin(interp[n].interp_mini[m].di_in, nxgrid_acc);
          acc_copyin(interp[n].interp_mini[m].dj_in, nxgrid_acc);
        }
      }

      free(t_in) ; free(ind);
      free(i_in) ; free(j_in); free(i_out); free(j_out);
      free(xgrid_area);
      if(opcode & CONSERVE_ORDER2) {
        free(xgrid_clon); free(xgrid_clat);
      }

    }//if file exists
  } //ntiles out

} /*end read_remap*/
