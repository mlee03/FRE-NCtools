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
#include <time.h>
#include "globals_acc.h"
#include "conserve_interp_acc.h"
#include "interp_utils_acc.h"
#include "create_xgrid_acc.h"
#include "create_xgrid_utils_acc.h"
#include "general_utils_acc.h"
#include "mpp.h"
#include "mpp_io.h"
#include "read_mosaic.h"

/*******************************************************************************
  void setup_conserve_interp
  Setup the interpolation weight for conservative interpolation
*******************************************************************************/
void setup_conserve_interp_acc(int ntiles_input_grid, Grid_config *input_grid, int ntiles_output_grid,
			   Grid_config *output_grid, Interp_config_acc *interp, unsigned int opcode)
{
  int nlon_input_cells, nlat_input_cells, ncells_input;
  int nlon_output_cells, nlat_output_cells, tile;

  Grid_cells_struct_config output_grid_cells;

  if( opcode & READ) {
    read_remap_file_acc(ntiles_input_grid, ntiles_output_grid, output_grid, input_grid, interp, opcode);
    copy_interp_to_device_acc(ntiles_input_grid, ntiles_output_grid, interp, opcode);
    return;
  }

  for(int n=0; n<ntiles_output_grid; n++) {

    int nlon_output_cells = output_grid[n].nxc;
    int nlat_output_cells = output_grid[n].nyc;
    int npts_output_grid = (nlon_output_cells+1)*(nlat_output_cells+1);

    interp[n].nxcells = 0;

    copy_grid_to_device_acc(npts_output_grid, output_grid[n].latc, output_grid[n].lonc);

    get_grid_cells_struct_acc( nlon_output_cells, nlat_output_cells,
                               output_grid[n].lonc, output_grid[n].latc, &output_grid_cells );

    for(int m=0; m<ntiles_input_grid; m++){

      int nlon_input_cells = input_grid[m].nx;
      int nlat_input_cells = input_grid[m].ny;
      int ncells_input = nlon_input_cells * nlat_input_cells;
      int npts_input_grid = (nlon_input_cells+1)*(nlat_input_cells+1);
      int jlat_overlap_starts=0, jlat_overlap_ends=0, nxcells=0, upbound_nxcells=0;
      int *approx_nxcells_per_ij1, *ij2_start, *ij2_end;
      double *skip_input_cells = (double *)malloc(ncells_input *sizeof(double));

      copy_grid_to_device_acc(npts_input_grid, input_grid[m].latc, input_grid[m].lonc);

      get_skip_cells_acc(nlon_input_cells*nlat_input_cells, skip_input_cells);

      //get the input grid portion (bounding indices) that overlaps with the output grid in the latitudonal direction.
      get_bounding_indices_acc(nlon_output_cells, nlat_output_cells, nlon_input_cells, nlat_input_cells,
                               output_grid[n].latc, input_grid[m].latc, &jlat_overlap_starts, &jlat_overlap_ends);

      create_upbound_nxcells_arrays_on_device_acc( ncells_input, &approx_nxcells_per_ij1, &ij2_start, &ij2_end);

      upbound_nxcells = get_upbound_nxcells_2dx2d_acc( nlon_input_cells, nlat_input_cells,
                                                       nlon_output_cells, nlat_output_cells,
                                                       jlat_overlap_starts, jlat_overlap_ends,
                                                       input_grid[m].lonc, input_grid[m].latc,
                                                       output_grid[n].lonc, output_grid[n].latc,
                                                       skip_input_cells,
                                                       &output_grid_cells,
                                                       approx_nxcells_per_ij1, ij2_start, ij2_end);

      if(opcode & GREAT_CIRCLE) {
        printf("GREAT_CIRCLE HAS NOT BEEN IMPLEMENTED YET\n");
        exit(0);
      }
      else {

        if(opcode & CONSERVE_ORDER1) {
          nxcells= create_xgrid_2dx2d_order1_acc( nlon_input_cells, nlat_input_cells,
                                                  nlon_output_cells, nlat_output_cells,
                                                  jlat_overlap_starts, jlat_overlap_ends,
                                                  input_grid[m].lonc, input_grid[m].latc,
                                                  output_grid[n].lonc, output_grid[n].latc,
                                                  upbound_nxcells,
                                                  skip_input_cells,
                                                  &output_grid_cells,
                                                  approx_nxcells_per_ij1, ij2_start, ij2_end,
                                                  interp[n].per_intile+m);
          interp[n].nxcells+=nxcells;
        }
        else if(opcode & CONSERVE_ORDER2) {
          nxcells= create_xgrid_2dx2d_order2_acc( nlon_input_cells, nlat_input_cells,
                                                  nlon_output_cells, nlat_output_cells,
                                                  jlat_overlap_starts, jlat_overlap_ends,
                                                  input_grid[m].lonc, input_grid[m].latc,
                                                  output_grid[n].lonc, output_grid[n].latc,
                                                  upbound_nxcells,
                                                  skip_input_cells,
                                                  &output_grid_cells,
                                                  approx_nxcells_per_ij1, ij2_start, ij2_end,
                                                  interp[n].per_intile+m, input_grid[m].cell_area);
          interp[n].nxcells+=nxcells;
        }
        else mpp_error("conserve_interp: interp_method should be CONSERVE_ORDER1 or CONSERVE_ORDER2");
      } //conserve_order methods

      free_upbound_nxcells_array_from_all_acc(ncells_input, approx_nxcells_per_ij1, ij2_start, ij2_end);
      free_skip_cells_on_all_acc( ncells_input, skip_input_cells);
      delete_grid_from_device_acc(npts_input_grid, input_grid[m].lonc, input_grid[m].latc);

    } //input tile

    free_grid_cell_struct_acc( nlon_output_cells*nlat_output_cells, &output_grid_cells);
    delete_grid_from_device_acc(npts_output_grid, output_grid[n].lonc, output_grid[n].latc);

  }//output tile

  if( opcode & WRITE) write_remap_file(ntiles_output_grid, ntiles_input_grid, output_grid, input_grid, interp, opcode);
  if(opcode & CHECK_CONSERVE) check_area_conservation(ntiles_output_grid, ntiles_input_grid, output_grid, interp);

  if(mpp_pe() == mpp_root_pe())printf("NOTE: done calculating index and weight for conservative interpolation\n");

}; /* setup_conserve_interp */


/*******************************************************************************
 void read_remap_file
 Reads in the weight/remap file if provided and copies the data to the device
*******************************************************************************/
void read_remap_file_acc(int ntiles_input_grid, int ntiles_output_grid, Grid_config *output_grid, Grid_config *input_grid,
                         Interp_config_acc *interp, unsigned int opcode)
{

  int *input_lon=NULL, *input_lat=NULL, *output_lon=NULL, *output_lat=NULL;
  double *xcell_area=NULL, *xcell_centroid_lon=NULL, *xcell_centroid_lat=NULL;

  int nxcells_acc;

  for(int n=0; n<ntiles_output_grid; n++) {
    if( interp[n].file_exist ) {
      int *t_in, *ind_acc;
      int fid, vid;

      int nlon_output = output_grid[n].nxc;
      int nxcells = read_mosaic_xgrid_size(interp[n].remap_file);
      interp[n].nxcells = nxcells;

      t_in  = (int *)malloc(nxcells*sizeof(int));
      input_lon = (int *)malloc(nxcells*sizeof(int));
      input_lat = (int *)malloc(nxcells*sizeof(int));
      output_lat = (int *)malloc(nxcells*sizeof(int));
      output_lon = (int *)malloc(nxcells*sizeof(int));
      xcell_area = (double *)malloc(nxcells*sizeof(double));
      if( opcode & CONSERVE_ORDER2) {
        xcell_centroid_lon = (double *)malloc(nxcells*sizeof(double));
        xcell_centroid_lat = (double *)malloc(nxcells*sizeof(double));
      }

      if(opcode & CONSERVE_ORDER1)
        read_mosaic_xgrid_order1(interp[n].remap_file, input_lon, input_lat, output_lon, output_lat, xcell_area);
      else
        read_mosaic_xgrid_order2(interp[n].remap_file, input_lon, input_lat, output_lon, output_lat, xcell_area,
                                 xcell_centroid_lon, xcell_centroid_lat);

      //rescale the xgrid area
      for(int i=0; i<interp[n].nxcells; i++) xcell_area[i] *= GAREA;

      //read in tile number of input parent cells
      fid = mpp_open(interp[n].remap_file, MPP_READ);
      vid = mpp_get_varid(fid, "tile1");
      mpp_get_var_value(fid, vid, t_in);
      mpp_close(fid);

      //tile number starts from 1, not 0, in the weight files
      for(int i=0 ; i<nxcells ; i++) t_in[i]--;

      //get number of nxcells per input tile
      for(int m=0 ; m<ntiles_input_grid ; m++) interp[n].per_intile[m].nxcells=0;
      for(int i=0 ; i<nxcells ; i++) interp[n].per_intile[ t_in[i] ].nxcells++;

      for(int m=0 ; m<ntiles_input_grid ; m++) {
        nxcells_acc = interp[n].per_intile[m].nxcells;
        interp[n].per_intile[m].input_parent_cell_indices  = (int *)malloc(nxcells_acc*sizeof(int));
        interp[n].per_intile[m].output_parent_cell_indices = (int *)malloc(nxcells_acc*sizeof(int));
        interp[n].per_intile[m].xcell_area  = (double *)malloc(nxcells_acc*sizeof(double));
        if(opcode & CONSERVE_ORDER2) {
          interp[n].per_intile[m].dcentroid_lon = (double *)malloc(nxcells_acc*sizeof(double));
          interp[n].per_intile[m].dcentroid_lat = (double *)malloc(nxcells_acc*sizeof(double));
        }
      }

      ind_acc = (int *)calloc(ntiles_input_grid, sizeof(int));
      for(int i=0 ; i<nxcells ; i++) {
        int itile, ii, nlon_input;
        itile = t_in[i];
        ii=ind_acc[itile];
        nlon_input = input_grid[itile].nxc;
        interp[n].per_intile[itile].input_parent_cell_indices[ii] = input_lat[i]*nlon_input + input_lon[i];
        interp[n].per_intile[itile].output_parent_cell_indices[ii] = output_lat[i]*nlon_output + output_lon[i];
        interp[n].per_intile[itile].xcell_area[ii] = xcell_area[i];
        if( opcode & CONSERVE_ORDER2) {
          interp[n].per_intile[itile].dcentroid_lon[ii] = xcell_centroid_lon[i];
          interp[n].per_intile[itile].dcentroid_lat[ii] = xcell_centroid_lat[i];
        }
        ind_acc[itile]++;
      }

      free(t_in) ; free(ind_acc);
      free(input_lon) ; free(input_lat); free(output_lon); free(output_lat);
      free(xcell_area); free(xcell_centroid_lon); free(xcell_centroid_lat);

    }//if file exists
  } //ntiles out

  printf("NOTE: Finish reading index and weight for conservative interpolation from file.\n");

} //end read_remap_file

/*******************************************************************************
 void write_remap_file
 write out the traditional remap file.
*******************************************************************************/
void write_remap_file(const int ntiles_out, const int ntiles_in, Grid_config *output_grid,
                      Grid_config *input_grid, Interp_config_acc *interp, unsigned int opcode)
{

  //copy and pasted from the original code with start and write

  for(int n=0 ; n<ntiles_out ; n++) {

    Interp_config_acc *p_interp = interp+n;
    int nxcells=p_interp->nxcells;
    int nlon_cells, ii;

    size_t start[4] = {0,0,0,0}, nwrite[4] = {1, 1, 1, 1};
    int *data_int;
    double *data_double;

    int fid = mpp_open( interp[n].remap_file, MPP_WRITE);
    int dim_string = mpp_def_dim(fid, "string", STRING);
    int dim_ncells = mpp_def_dim(fid, "ncells", nxcells);
    int dim_two    = mpp_def_dim(fid, "two", 2);
    int dims[4] = {dim_ncells, dim_two, 0, 0};

    int id_tile1 = mpp_def_var(fid, "tile1", NC_INT, 1, &dim_ncells, 1,
                               "standard_name", "tile_number_in_mosaic1");
    int id_tile1_cell = mpp_def_var(fid, "tile1_cell", NC_INT, 2, dims, 1,
                                    "standard_name", "parent_cell_indices_in_mosaic1");
    int id_tile2_cell = mpp_def_var(fid, "tile2_cell", NC_INT, 2, dims, 1,
                                    "standard_name", "parent_cell_indices_in_mosaic2");
    int id_xgrid_area = mpp_def_var(fid, "xgrid_area", NC_DOUBLE, 1, &dim_ncells, 2,
                                    "standard_name", "exchange_grid_area", "units", "m2");
    int id_tile1_dist = (opcode & CONSERVE_ORDER2) ?
      mpp_def_var(fid, "tile1_distance", NC_DOUBLE, 2, dims, 1,
                  "standard_name", "distance_from_parent1_cell_centroid") : 0 ;

    nwrite[0] = nxcells;
    mpp_end_def(fid);

    data_int = (int *)malloc(nxcells*sizeof(int));

    //update data on host
    for(int m=0 ; m<ntiles_in ; m++) {
      Interp_per_input_tile *p_interp_for_tile_m = p_interp->per_intile+m;
      int m_nxcells = p_interp_for_tile_m->nxcells;
#pragma acc update host( p_interp_for_tile_m->input_parent_cell_indices[:m_nxcells], \
                         p_interp_for_tile_m->output_parent_cell_indices[:m_nxcells], \
                         p_interp_for_tile_m->xcell_area[:m_nxcells])
#pragma acc update if(opcode &CONSERVE_ORDER2) host(p_interp_for_tile_m->dcentroid_lon[:m_nxcells], \
                                                    p_interp_for_tile_m->dcentroid_lat[:m_nxcells])
    }

    //input tile
    ii = 0;
    for( int m=0 ; m<ntiles_in ; m++ ) {
      Interp_per_input_tile *p_interp_for_tile_m = p_interp->per_intile+m;
      int m_nxcells = p_interp_for_tile_m->nxcells;
      for( int i=0 ; i<m_nxcells ; i++ ){
        data_int[ii] = m+1;
        ii++;
      }
    }
    mpp_put_var_value(fid, id_tile1, data_int);

    // i (x, lon) indices of input parent
    ii=0;
    for( int m=0 ; m<ntiles_in ; m++){
      Interp_per_input_tile *p_interp_for_tile_m = p_interp->per_intile+m;
      int m_nxcells = p_interp_for_tile_m->nxcells;
      nlon_cells = input_grid[m].nxc;
      for( int i=0 ; i<m_nxcells; i++) {
        data_int[ii] = p_interp_for_tile_m->input_parent_cell_indices[i]%nlon_cells+1;
        ii++;
      }
    }
    mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, data_int);

    // i (x, lon) indices of output parent
    ii=0;
    nlon_cells = output_grid[n].nxc;
    for( int m=0 ; m<ntiles_in ; m++){
      Interp_per_input_tile *p_interp_for_tile_m = p_interp->per_intile+m;
      int m_nxcells = p_interp_for_tile_m->nxcells;
      for( int i=0 ; i<m_nxcells; i++) {
        data_int[ii] = p_interp_for_tile_m->output_parent_cell_indices[i]%nlon_cells+1;
        ii++;
      }
    }
    mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, data_int);

    start[1]=1;

    // j (y, lat) indices of input parent
    ii=0;
    for( int m=0 ; m<ntiles_in ; m++){
      Interp_per_input_tile *p_interp_for_tile_m = p_interp->per_intile+m;
      int m_nxcells = p_interp_for_tile_m->nxcells;
      nlon_cells = input_grid[m].nxc;
      for( int i=0 ; i<m_nxcells; i++) {
        data_int[ii] = p_interp_for_tile_m->input_parent_cell_indices[i]/nlon_cells+1;
        ii++;
      }
    }
    mpp_put_var_value_block(fid, id_tile1_cell, start, nwrite, data_int);

    // j (y, lat) indices of output parent
    ii=0;
    nlon_cells = output_grid[n].nxc;
    for( int m=0 ; m<ntiles_in ; m++){
      Interp_per_input_tile *p_interp_for_tile_m = p_interp->per_intile+m;
      int m_nxcells = p_interp_for_tile_m->nxcells;
      for( int i=0 ; i<m_nxcells; i++) {
        data_int[ii] = p_interp_for_tile_m->output_parent_cell_indices[i]/nlon_cells+1;
        ii++;
      }
    }
    mpp_put_var_value_block(fid, id_tile2_cell, start, nwrite, data_int);

    free(data_int);
    data_double = (double *)malloc(nxcells*sizeof(double));

    // exchange cell area
    ii=0;
    for( int m=0 ; m<ntiles_in ; m++){
      Interp_per_input_tile *p_interp_for_tile_m = p_interp->per_intile+m;
      int m_nxcells = p_interp_for_tile_m->nxcells;
      for( int i=0 ; i<m_nxcells; i++) {
        data_double[ii] = p_interp_for_tile_m->xcell_area[i];
        ii++;
      }
    }
    mpp_put_var_value(fid, id_xgrid_area, data_double);

    if(opcode & CONSERVE_ORDER2) {
      ii=0; start[1] = 0 ;
      for( int m=0 ; m<ntiles_in ; m++) {
        Interp_per_input_tile *p_interp_for_tile_m = p_interp->per_intile+m;
        int m_nxcells = p_interp_for_tile_m->nxcells;
        for( int i=0 ; i<m_nxcells ; i++) {
          data_double[ii] = p_interp_for_tile_m->dcentroid_lon[i];
          ii++;
        }
      }
      mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, data_double);

      ii=0; start[1] = 1 ;
      for( int m=0 ; m<ntiles_in ; m++) {
        Interp_per_input_tile *p_interp_for_tile_m = p_interp->per_intile+m;
        int m_nxcells = p_interp_for_tile_m->nxcells;
        for( int i=0 ; i<m_nxcells ; i++) {
          data_double[ii] = p_interp_for_tile_m->dcentroid_lat[i];
          ii++;
        }
      }
      mpp_put_var_value_block(fid, id_tile1_dist, start, nwrite, data_double);
    }
    free(data_double);

    mpp_close(fid);

  }//ntiles_out

}

void check_area_conservation(const int ntiles_output_grid, const int ntiles_input_grid, Grid_config *output_grid,
                             Interp_config_acc *interp)
{

  for(int n=0; n<ntiles_output_grid; n++) {

    int nlon_cells = output_grid[n].nxc;
    int nlat_cells = output_grid[n].nyc;
    int ncells = nlon_cells*nlat_cells;

    int max_ij=0;
    double max_ratio=0.0, ratio_change=0.0, max_area=0.0;
    double *recomputed_output_area = (double *)calloc(ncells,sizeof(double));

    /* sum over exchange grid to get the area of output grid cells*/
    for(int m=0; m<ntiles_input_grid; m++) {
      Interp_per_input_tile *minterp = interp[n].per_intile+m;
      int nxcells = minterp->nxcells;
      for(int i=0; i<nxcells; i++) {
        int ii = minterp->output_parent_cell_indices[i];
        recomputed_output_area[ii] += minterp->xcell_area[i];
      }
    }

    /* compare actual area and recomputed_output_area */
    for(int ij=0 ; ij<ncells ; ij++) {
      double actual_cell_area=output_grid[n].cell_area[ij];
      ratio_change = fabs(actual_cell_area-recomputed_output_area[ij])/actual_cell_area;
      if(ratio_change > max_ratio) {
        max_ratio = ratio_change;
        max_area = recomputed_output_area[ij];
        max_ij = ij;
      }
      if( ratio_change > 1.e-4 ) {
        printf("(i,j)=(%d,%d), change = %g, area1=%g, recomputed_output_area=%g\n",
               ij%nlon_cells, ij/nlon_cells, ratio_change, output_grid[n].cell_area[ij],recomputed_output_area[ij]);
      }
    }
    printf("The maximum ratio change at (%d,%d) = %g, area1=%g, recomputed_output_area=%g\n",
           max_ij%nlon_cells, max_ij/nlon_cells, max_ratio, output_grid[n].cell_area[max_ij], max_area);

    free(recomputed_output_area);
  }//for each output tile

}

/*******************************************************************************
 void do_scalar_conserve_interp( )
 doing conservative interpolation
*******************************************************************************/
void do_scalar_conserve_interp_acc(Interp_config_acc *interp, int varid, int ntiles_input_grid, const Grid_config *input_grid,
                                          int ntiles_output_grid, const Grid_config *output_grid, const Field_config *field_in,
                                          Field_config *field_out, unsigned int opcode, int nz)
{
  int weights_exist = input_grid[0].weight_exist;;
  int cell_measures = field_in->var[varid].cell_measures;
  int cell_methods = field_in->var[varid].cell_methods;
  int target_grid = ( field_in->var[varid].use_volume ) ? 0 : (opcode & TARGET);
  int has_missing = field_in->var[varid].has_missing;
  double missing = (has_missing) ? field_in->var[varid].missing : -MAXVAL;
  double gsum_out=0.0;

  for(int n=0; n<ntiles_output_grid; n++) {

    int output_ncells = output_grid[n].nxc * output_grid[n].nyc;
    int *out_miss = (int *)malloc(output_ncells*sizeof(int));
    double *out_area = (double *)malloc(output_ncells*sizeof(double));
    double *p_fieldout_data = field_out[n].data;

#pragma acc enter data create( p_fieldout_data[:output_ncells],   \
                               out_area[:output_ncells],          \
                               out_miss[:output_ncells])

#pragma acc parallel loop present(p_fieldout_data[:output_ncells], out_area[:output_ncells], out_miss[:output_ncells])
    for(int i=0; i<output_ncells; i++) {
      p_fieldout_data[i] = 0.0;
      out_area[i] = 0.0;
      out_miss[i] = 0;
    }

    for(int m=0 ; m<ntiles_input_grid; m++) {

      int input_ncells = input_grid[m].nxc * input_grid[m].nyc;
      double *input_area_weight = (double *)malloc(input_ncells*sizeof(double));

#pragma acc enter data create(input_area_weight[:input_ncells])
      get_input_area_weight(weights_exist, cell_measures, cell_methods, field_in+m, input_grid+m, input_area_weight);

      if(opcode & CONSERVE_ORDER1)
        interp_data_order1(output_grid+n, input_grid+m, interp[n].per_intile+m, input_area_weight, field_in[m].data,
                           p_fieldout_data, out_area, out_miss, missing);
      if(opcode & CONSERVE_ORDER2)
        interp_data_order2(output_grid+n, input_grid+m, interp[n].per_intile+m, input_area_weight, field_in[m].data,
                           p_fieldout_data, out_area, out_miss,
                           field_in[m].grad_mask, field_in[m].grad_y, field_in[m].grad_x,  missing);

#pragma acc exit data delete(input_area_weight[:input_ncells])

    } //m

    if(opcode & CHECK_CONSERVE) {
#pragma acc enter data copyin(gsum_out)
#pragma acc parallel loop present(out_area[:output_ncells], p_fieldout_data[:output_ncells]) reduction(+:gsum_out)
      for(int i=0; i<output_ncells; i++) {
        if(out_area[i] > 0) gsum_out += p_fieldout_data[i];
      }
    }

    if ( cell_methods == CELL_METHODS_SUM ) {
#pragma acc parallel loop present(out_area[:output_ncells], out_miss[:output_ncells], p_fieldout_data[:output_ncells])
      for(int i=0; i<output_ncells; i++) {
        if(out_area[i] == 0) {
          p_fieldout_data[i] = 0.0;
          if(out_miss[i] == 0) p_fieldout_data[i] = missing;
        }
      }
    }
    else {
#pragma acc parallel loop present(out_area[:output_ncells], out_miss[:output_ncells], p_fieldout_data[:output_ncells])
      for(int i=0; i<output_ncells; i++) {
        if(out_area[i] > 0) {
          p_fieldout_data[i] /= out_area[i];
        }
        else {
          p_fieldout_data[i] = 0.0;
          if(out_miss[i] == 0) p_fieldout_data[i] = missing;
        }
      }

      if( (target_grid) ) {
#pragma acc parallel loop present(out_area[:output_ncells])
        for(int i=0; i<output_ncells; i++) out_area[i] = 0.0;

        for(int m=0 ; m<ntiles_input_grid; m++) {
          Interp_per_input_tile *minterp = interp[n].per_intile+m;
          int input_ncells = input_grid[m].nxc * input_grid[m].nyc;
          int ixcells = minterp->nxcells;
          double *p_gridin_area = input_grid[m].cell_area;
          double *p_fieldin_area = field_in[m].area;
#pragma acc parallel loop present(minterp->output_parent_cell_indices[:ixcells], \
                                  minterp->input_parent_cell_indices[:ixcells], \
                                  minterp->xcell_area[:ixcells],         \
                                  out_area[:output_ncells])\
                           copyin(p_fieldin_area[:input_ncells],\
                                  p_gridin_area[:input_ncells])
          for(int ix=0; ix<ixcells; ix++) {
            int ij2   = minterp->output_parent_cell_indices[ix];
            int ij1   = minterp->input_parent_cell_indices[ix];
            double area = minterp->xcell_area[ix];
            if(cell_measures ) out_area[ij2] += (area*p_fieldin_area[ij1]/p_gridin_area[ij1]);
            else out_area[ij2] += area;
          }
#pragma acc parallel loop present(p_fieldout_data[:output_ncells], out_area[:output_ncells], \
                                  output_grid[n].cell_area[:output_ncells])
          for(int i=0; i<output_ncells; i++) {
            if(p_fieldout_data[i] != missing)
              p_fieldout_data[i] *=  (out_area[i]/output_grid[n].cell_area[i]);
          }
        }
      }
    }

#pragma acc exit data copyout(p_fieldout_data[:output_ncells])
#pragma acc exit data delete(out_area[:output_ncells],\
                             out_miss[:output_ncells])

    free(out_area);
    free(out_miss);
  } // n

  return;

  /* conservation check if needed */
  if(opcode & CHECK_CONSERVE) {
    double gsum_in = 0.0;
    for(int m=0; m<ntiles_input_grid; m++) {

      int nx1  = input_grid[m].nx;
      int ny1  = input_grid[m].ny;

      if( cell_measures ) {
        for(int ij=0; ij<nx1*ny1; ij++){
          double dd = field_in[m].data[ij];
          if(dd != missing) gsum_in += dd*field_in[m].area[ij];
        }
        continue;
      }
      if ( cell_methods == CELL_METHODS_SUM ) {
        for(int ij=0; ij<nx1*ny1; ij++) {
          double dd = field_in[m].data[ij];
          if(dd != missing) gsum_in += dd;
        }
        continue;
      }
      else {
        for(int ij=0; ij<nx1*ny1; ij++) {
          double dd = field_in[m].data[ij];
          if(dd != missing) gsum_in += dd*input_grid[m].cell_area[ij];
        }
      }
    }
    mpp_sum_double(1, &gsum_out);

    if(mpp_pe() == mpp_root_pe()) printf("the flux(data*area) sum of %s: input = %g, output = %g, diff = %g. \n",
                                         field_in->var[varid].name, gsum_in, gsum_out, gsum_out-gsum_in);

  }


}; /* do_scalar_conserve_interp_order */

void get_input_area_weight(const int weights_exist, const int cell_measures, const int cell_methods,
                           const Field_config *field_in, const Grid_config *input_grid,
                           double *input_area_weight)
{
  int input_ncells = input_grid->nxc * input_grid->nyc;
  double *p_gridin_area  = input_grid->cell_area;
  double *p_fieldin_area = field_in->area;
  double *p_weight = input_grid->weight;

  if(cell_methods == CELL_METHODS_SUM) {
#pragma acc parallel loop present(input_area_weight[:input_ncells]) copyin(p_gridin_area[:input_ncells])
    for(int i=0 ; i<input_ncells ; i++) input_area_weight[i] = 1.0/p_gridin_area[i];
  }

  else if(cell_measures) {
#pragma acc parallel loop present(input_area_weight[:input_ncells]) \
                          copyin(p_gridin_area[:input_ncells],     \
                                 p_fieldin_area[:input_ncells])
    for(int i=0 ; i<input_ncells ; i++) input_area_weight[i] = p_fieldin_area[i]/p_gridin_area[i];
  }

  else {
#pragma acc parallel loop present(input_area_weight[:input_ncells])
    for(int i=0 ; i<input_ncells ; i++) input_area_weight[i]=1.0;
  }

  if(weights_exist){
#pragma acc parallel loop independent present(input_area_weight[:input_ncells]) copyin(p_weight[:input_ncells])
    for(int i=0; i<input_ncells; i++) input_area_weight[i] *= p_weight[i];
  }

}

void interp_data_order1( const Grid_config *output_grid, const Grid_config *input_grid,
                         Interp_per_input_tile *minterp, double *input_area_weight, double *fieldin_data,
                         double *fieldout_data, double *out_area, int *out_miss, double missing)
{

  int nxcells = minterp->nxcells;
  int input_ncells = input_grid->nxc * input_grid->nyc;
  int output_ncells = output_grid->nxc * output_grid->nyc;

#pragma acc data present(minterp[:1],                                    \
                         minterp->input_parent_cell_indices[:nxcells],   \
                         minterp->output_parent_cell_indices[:nxcells],  \
                         minterp->xcell_area[:nxcells],                  \
                         input_area_weight[:input_ncells],              \
                         fieldout_data[:output_ncells],                 \
                         out_area[:output_ncells],                      \
                         out_miss[:output_ncells])                      \
  copyin(fieldin_data[:input_ncells])
#pragma acc parallel loop
  for(int ix=0; ix<nxcells; ix++) {
    int ij1 = minterp->input_parent_cell_indices[ix];
    int ij2 = minterp->output_parent_cell_indices[ix];
    double area = minterp->xcell_area[ix];

    if( fieldin_data[ij1] == missing ) continue;

    area *= input_area_weight[ij1] ;
#pragma acc atomic update
    fieldout_data[ij2] += fieldin_data[ij1]*area ;
#pragma acc atomic update
    out_area[ij2] += area;
    out_miss[ij2] = 1;
  }

}

void interp_data_order2( const Grid_config *output_grid, const Grid_config *input_grid,
                         Interp_per_input_tile *minterp, double *input_area_weight, double *fieldin_data,
                         double *fieldout_data, double *out_area, int *out_miss,
                         int *grad_mask, double *grad_y, double *grad_x, double missing)
{

  int nxcells = minterp->nxcells;
  int n_halo_cells = 2;
  int input_nlon_cells = input_grid->nxc;
  int input_nlat_cells = input_grid->nyc;
  int input_data_ncells = (input_nlon_cells+n_halo_cells)*(input_nlat_cells+n_halo_cells);
  int input_ncells = input_nlon_cells * input_nlat_cells;

  int output_nlon_cells = output_grid->nxc;
  int output_ncells = output_nlon_cells * (output_grid->nyc);

#pragma acc data present( minterp[:1],                                   \
                          minterp->input_parent_cell_indices[:nxcells],  \
                          minterp->output_parent_cell_indices[:nxcells], \
                          minterp->xcell_area[:nxcells],                 \
                          input_area_weight[:input_ncells],             \
                          fieldout_data[:output_ncells],                \
                          out_area[:output_ncells],                     \
                          out_miss[:output_ncells])                     \
  copyin(fieldin_data[:input_data_ncells],                              \
         grad_mask[:input_ncells],                                      \
         grad_x[:input_ncells], grad_y[:input_ncells])
#pragma acc parallel loop
  for(int ix=0; ix<nxcells; ix++){
    int ij1 = minterp->input_parent_cell_indices[ix];
    int ij2 = minterp->output_parent_cell_indices[ix];
    double area = minterp->xcell_area[ix];
    double dx = minterp->dcentroid_lon[ix];
    double dy = minterp->dcentroid_lat[ix];

    int i1 = ij1%input_nlon_cells;
    int j1=ij1/input_nlon_cells;
    int data_pt = (j1+1)*(input_nlon_cells+n_halo_cells)+(i1+1);

    if( fieldin_data[data_pt] == missing ) continue;

    area *= input_area_weight[ij1] ;
#pragma acc atomic update
    fieldout_data[ij2] += (fieldin_data[data_pt] + (1-grad_mask[ij1])*
                           (grad_x[ij1]*dx + grad_y[ij1]*dy) )*area;
#pragma acc atomic update
    out_area[ij2] += area;
    out_miss[ij2] = 1;
  }

}
