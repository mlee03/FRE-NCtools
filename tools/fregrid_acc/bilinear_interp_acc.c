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
#include <time.h>
#include <openacc.h>
#include "globals_acc.h"
#include "general_utils_acc.h"
#include "bilinear_interp_acc.h"
#include "interp_utils_acc.h"
#include "mpp_io.h"
#include "mpp.h"

#define min(a,b) (a<b ? a:b)
#define max(a,b) (a>b ? a:b)
#define sign(a,b)(b<0 ? -fabs(a):fabs(a))

#pragma acc routine seq
int max_weight_index_acc( double *var, int nvar);
double normalize_great_circle_distance_acc(const double *v1, const double *v2);
#pragma acc routine seq
double dist2side_acc(const double *v1, const double *v2, const double *point);
void redu2x_acc(const double *varfin, const double *yfin, int nxfin, int nyfin, double *varcrs,
                int nxcrs, int nycrs, int has_missing, double missvalue);
void do_latlon_coarsening_acc(const double *var_latlon, const double *ylat, int nlon, int nlat, int nz,
			  double *var_latlon_crs, int finer_steps, int has_missing, double missvalue);
void do_c2l_interp_acc(const Interp_config_acc *interp_acc, int nlon_input_cells, int nlat_input_cells,
                       int input_ntiles, const Field_config *field_in,
                       int nx_out, int ny_out, double *data_out, int has_missing, double missing, int fill_missing);
#pragma acc routine seq
int get_index_acc(const Grid_config *input_grid, const Grid_config *output_grid, int *index,
                  int i_in, int j_in, int l_in, int i_out, int j_out);
#pragma acc routine seq
int get_closest_index_acc(const Grid_config *input_grid, const Grid_config *output_grid, int *index,
                          int i_in, int j_in, int l_in, int i_out, int j_out);

void read_remap_file(const Grid_config *output_grid, const Interp_config_acc *interp_acc);

void get_interp_index_test( const int input_ntiles, const int iter, const Grid_config *input_grid,
                            const Grid_config *output_grid, const double dlon_in, const double dlat_in,
                            const double lonbegin, const double latbegin, Interp_config_acc *interp_acc);
void get_interp_weights(const int input_ntiles, const int output_ntiles, const Grid_config *output_grid,
                        const Grid_config *input_grid, Interp_config_acc *interp_acc);
#pragma acc routine seq
void get_vector(const Grid_config *grid, const int icell, double *vector);
void write_bilinear_remap_file(const int nlon_output_cells, const int nlat_output_cells, Interp_config_acc *interp_acc);
void write_bilinear_interp_in_not_zhi_format(const int ncells, Interp_config_acc *interp_acc);

/*******************************************************************************
  void setup_bilinear_interp( )
    !------------------------------------------------------------------!
    ! calculate weights for bilinear interpolation                     !
    ! from cubed sphere to latlon grid                                 !
    !                                                                  !
    ! input:                                                           !
    ! sph_corner      cubed sphere corner location in spherical coor   !
    ! npx, npy        number of corners per tile                       !
    ! ntiles          number of tiles                                  !
    ! xlon, ylat      latlon grid coor                                 !
    ! nlon, nlat      latlon grid dimension                            !
    !                                                                  !
    ! output:                                                          !
    ! c2l_index       cubed sphere index for latlon interpolation      !
    ! c2l_weight      weights for cubsph_to_latlon interpolation       !
    ! elon_cubsph     lon unit vector for cubed sphere center          !
    ! elat_cubsph     lat unit vector for cubed sphere center          !
    ! elon_latlon     lon unit vector for latlon grid                  !
    ! elat_latlon     lat unit vector for latlon grid                  !
    !------------------------------------------------------------------!
*******************************************************************************/
void setup_bilinear_interp_acc(int input_ntiles, const Grid_config *input_grid,
                               int output_ntiles, const Grid_config *output_grid,
                               Interp_config_acc *interp_acc, unsigned int opcode, double dlon_in, double dlat_in,
                               double lonbegin, double latbegin)
{
  const int max_iter = 10;
  int nlon_output_cells = output_grid->nx_fine;
  int nlat_output_cells = output_grid->ny_fine;
  int ncells_output = nlon_output_cells * nlat_output_cells;
  int nxd = input_grid->nx + 2;
  int nyd = input_grid->ny + 2;

  // input_ntiles must be six and output_ntiles must be one
  if(input_ntiles != 6) mpp_error("Error from bilinear_interp: source mosaic should be cubic mosaic "
                      			       "and have six tiles when using bilinear option");
  if(output_ntiles != 1) mpp_error("Error from bilinear_interp: destination mosaic should be "
                            				"one tile lat-lon grid when using bilinear option");

  /*------------------------------------------------------------------!
   ! cubed sphere: cartesian coordinates of cell corners,             !
   !               cell lenghts between corners,                      !
   !               cartesian and spherical coordinates of cell centers!
   !               calculate latlon unit vector                       !
   !-----------------------------------------------------------------*/


  interp_acc->index  = (int    *)malloc(3*ncells_output*sizeof(int));
  interp_acc->weight = (double *)malloc(4*ncells_output*sizeof(double));

  // read from file
  if( (opcode & READ) && interp_acc->file_exist ) {
    read_remap_file(output_grid, interp_acc);
    return;
  }

  enter_bilinear_interp_to_device_acc(ncells_output, interp_acc);
  copy_bilinear_grid_to_device_acc(output_ntiles, ncells_output, output_grid);
  copy_bilinear_grid_to_device_acc(input_ntiles, nxd*nyd, input_grid);


  //find lower left corner on cubed sphere for given latlon location
  get_interp_index_test(input_ntiles, max_iter, input_grid, output_grid, dlon_in, dlat_in, lonbegin, latbegin, interp_acc);

  // calculate weights for interpolation
  get_interp_weights(input_ntiles, output_ntiles, output_grid, input_grid, interp_acc);


#pragma acc update host(interp_acc->index[:3*ncells_output], interp_acc->weight[:4*ncells_output])


  /* write out weight information if needed */
  if( opcode & WRITE ) write_bilinear_remap_file(nlon_output_cells, nlat_output_cells, interp_acc);
  write_bilinear_interp_in_not_zhi_format(ncells_output, interp_acc);

  delete_bilinear_grid_from_device_acc(output_ntiles, ncells_output, output_grid);
  delete_bilinear_grid_from_device_acc(input_ntiles, nxd*nyd, input_grid);

  printf("\n done calculating interp_index and interp_weight\n");

}; // setup_bilinear_interp

/*----------------------------------------------------------------------------
   void do_scalar_bilinear_interp(Mosaic_config *input, Mosaic_config *output, int varid )
   interpolate scalar data to latlon,                               !
   --------------------------------------------------------------------------*/
void do_scalar_bilinear_interp_acc(const Interp_config_acc *interp_acc, int vid, int input_ntiles,
                                   const Grid_config *input_grid, const Grid_config *output_grid,
                                   const Field_config *field_in, Field_config *field_out, int finer_step,
                                   int fill_missing)
{

  int nlon_output_cells = output_grid->nx_fine;
  int nlat_output_cells = output_grid->ny_fine;
  int nlon_input_cells  = input_grid->nx;
  int nlat_input_cells  = input_grid->ny;
  // currently we are regridding one vertical level for each call to reduce the memory usage
  double missing  = field_in[0].var[vid].missing;
  int has_missing = field_in[0].var[vid].has_missing;

  double *data_fine;  data_fine = (double *)malloc(nlon_output_cells*nlat_output_cells*sizeof(double));

#pragma acc enter data create(data_fine[:nlon_output_cells*nlat_output_cells])
#pragma acc enter data copyin(field_in[:input_ntiles])
  for(int itile=0 ; itile<input_ntiles ; itile++) {
#pragma acc enter data copyin(field_in[itile].data[:(nlon_input_cells+2)*(nlat_input_cells+2)])
  }


  do_c2l_interp_acc(interp_acc, nlon_input_cells, nlat_input_cells, input_ntiles, field_in, nlon_output_cells,
                    nlat_output_cells, data_fine, has_missing, missing, fill_missing);


#pragma acc exit data copyout(data_fine[:nlon_output_cells*nlat_output_cells])
  for(int itile=0 ; itile<input_ntiles ; itile++) {
#pragma acc exit data delete(field_in[itile].data[:(nlon_input_cells+2)*(nlat_input_cells+2)])
  }
#pragma acc exit data delete(field_in[:input_ntiles])


  do_latlon_coarsening_acc(data_fine, output_grid->latt1D_fine, nlon_output_cells, nlat_output_cells, 1, field_out->data,
                           finer_step, has_missing, missing);

  free(data_fine);

}; /* do_c2l_scalar_interp */



/*----------------------------------------------------------------------------
   void do_vector_bilinear_interp()
   interpolate vector data to latlon,                               !
   --------------------------------------------------------------------------*/
void do_vector_bilinear_interp_acc(Interp_config_acc *interp_acc, int vid, int input_ntiles,
                                   const Grid_config *input_grid, int output_ntiles,
                                   const Grid_config *output_grid, const Field_config *u_in,  const Field_config *v_in,
                                   Field_config *u_out, Field_config *v_out, int finer_step, int fill_missing)
{

  printf("BILINEAR VECTOR IS NOT AVAILABLE YET\n");

}; /* do_vector_bilinear_interp */


void do_c2l_interp_acc(const Interp_config_acc *interp_acc, int nlon_input_cells, int nlat_input_cells, int input_ntiles,
                       const Field_config *field_in, int nlon_output_cells, int nlat_output_cells,
                       double *data_out, int has_missing, double missing, int fill_missing)
{
  int ind;
  int nxd = nlon_input_cells + 2;
  int ncells_output = nlon_output_cells * nlat_output_cells;

  if (has_missing) {
    if(fill_missing) {
#pragma acc parallel loop present(interp_acc[:1], field_in[:input_ntiles], data_out[:ncells_output])
      for(int n2=0; n2<ncells_output; n2++) {
        int ic = interp_acc->index[3*n2];
        int jc = interp_acc->index[3*n2+1];
        int tile = interp_acc->index[3*n2+2];
        double d_in[4] = {field_in[tile].data[jc*nxd+ic],
                          field_in[tile].data[(jc+1)*nxd+ic],
                          field_in[tile].data[(jc+1)*nxd+ic+1],
                          field_in[tile].data[jc*nxd+ic+1]};
        if (d_in[0] == missing || d_in[1] == missing || d_in[2] == missing || d_in[3] == missing ) {
          ind = max_weight_index_acc(interp_acc->weight+4*n2, 4);
          data_out[n2] = d_in[ind];
          continue;
        }
        data_out[n2] = d_in[0]*interp_acc->weight[4*n2] + d_in[1]*interp_acc->weight[4*n2+1] +
                       d_in[2]*interp_acc->weight[4*n2+2] + d_in[3]*interp_acc->weight[4*n2+3];
      }
    } // if fill_missing

    else {
#pragma acc parallel loop present(interp_acc[:1], field_in[:input_ntiles], data_out[:ncells_output])
      for(int n2=0; n2<ncells_output; n2++) {
        int ic = interp_acc->index[3*n2];
        int jc = interp_acc->index[3*n2+1];
        int tile = interp_acc->index[3*n2+2];
        double d_in[4] = {field_in[tile].data[jc*nxd+ic],
                          field_in[tile].data[(jc+1)*nxd+ic],
                          field_in[tile].data[(jc+1)*nxd+ic+1],
                          field_in[tile].data[jc*nxd+ic+1]};
        if (d_in[0] == missing || d_in[1] == missing || d_in[2] == missing || d_in[3] == missing ) {
          data_out[n2] = missing;
          continue;
        }
        data_out[n2] = d_in[0]*interp_acc->weight[4*n2] + d_in[1]*interp_acc->weight[4*n2+1] +
                       d_in[2]*interp_acc->weight[4*n2+2] + d_in[3]*interp_acc->weight[4*n2+3];
      }
    } //if not fill_missing

  } // if missing

  else {
#pragma acc parallel loop present(interp_acc[:1], field_in[:input_ntiles], data_out[:ncells_output])
    for(int n2=0; n2<ncells_output; n2++) {
      int ic = interp_acc->index[3*n2];
      int jc = interp_acc->index[3*n2+1];
      int tile = interp_acc->index[3*n2+2];
      double d_in[4] = {field_in[tile].data[jc*nxd+ic],
                        field_in[tile].data[(jc+1)*nxd+ic],
                        field_in[tile].data[(jc+1)*nxd+ic+1],
                        field_in[tile].data[jc*nxd+ic+1]};
      data_out[n2] = d_in[0]*interp_acc->weight[4*n2] + d_in[1]*interp_acc->weight[4*n2+1] +
                     d_in[2]*interp_acc->weight[4*n2+2] + d_in[3]*interp_acc->weight[4*n2+3];
    }
  } // if not missing

}; /* do_c2l_interp_acc */


/*------------------------------------------------------------------
  void get_index_acc(ig, jg, lg)
  determine lower left corner
  ----------------------------------------------------------------*/
int get_index_acc(const Grid_config *input_grid, const Grid_config *output_grid, int *index,
	       int i_in, int j_in, int l_in, int i_out, int j_out)
{
  int    ok, n0, n1, n2, n3, n4, n5;
  int    nlon_input_cells, nlat_input_cells, nlon_output_cells, nlat_output_cells;
  double v0[3], v1[3], v2[3], v3[3], v4[3], v5[3];
  double angle_1_2_3, angle_1_2_0, angle_1_3_0;
  double angle_1_3_4, angle_1_4_0;
  double angle_1_4_5, angle_1_5_0, angle_1_5_2;

  ok=1;
  nlon_input_cells  = input_grid->nx_fine;
  nlat_input_cells  = input_grid->nx_fine;
  nlon_output_cells = output_grid->nx;
  nlat_output_cells = output_grid->nx;
  n0 = j_out*nlon_output_cells + i_out;
  n1 = j_in*nlon_input_cells + i_in;
  n2 = j_in*nlon_input_cells + i_in+1;
  n3 = (j_in+1)*nlon_input_cells + i_in;
  v0[0] = output_grid->xt[n1];
  v0[1] = output_grid->yt[n1];
  v0[2] = output_grid->zt[n1];
  v1[0] = input_grid->xt[n1];
  v1[1] = input_grid->yt[n1];
  v1[2] = input_grid->zt[n1];
  v2[0] = input_grid->xt[n2];
  v2[1] = input_grid->yt[n2];
  v2[2] = input_grid->zt[n2];
  v3[0] = input_grid->xt[n3];
  v3[1] = input_grid->yt[n3];
  v3[2] = input_grid->zt[n3];
  angle_1_2_3 = spherical_angle_acc(v1, v2, v3);
  angle_1_2_0= spherical_angle_acc(v1, v2, v0);
  angle_1_3_0= spherical_angle_acc(v1, v3, v0);

  if (max(angle_1_2_0,angle_1_3_0)<angle_1_2_3) {
    index[0]=i_in;
    index[1]=j_in;
    index[2]=l_in;
  }
  else {
    n4 = j_in*nlon_input_cells + i_in-1;
    v4[0] = input_grid->xt[n4];
    v4[1] = input_grid->yt[n4];
    v4[2] = input_grid->zt[n4];
    angle_1_3_4 =spherical_angle_acc(v1, v3, v4);
    angle_1_3_0=spherical_angle_acc(v1, v3, v0);
    angle_1_4_0=spherical_angle_acc(v1, v4, v0);
    if (max(angle_1_3_0,angle_1_4_0)<angle_1_3_4) {
      index[0]=i_in-1;
      index[1]=j_in;
      index[2]=l_in;
    }
    else {
      n5 = (j_in-1)*nlon_input_cells + i_in;
      v5[0] = input_grid->xt[n5];
      v5[1] = input_grid->yt[n5];
      v5[2] = input_grid->zt[n5];
      angle_1_4_5 =spherical_angle_acc(v1, v4, v5);
      angle_1_5_0=spherical_angle_acc(v1, v5, v0);
      if (max(angle_1_4_0,angle_1_5_0)<angle_1_4_5 && i_in>1 && j_in>1) {
        index[0]=i_in-1;
        index[1]=j_in-1;
        index[2]=l_in;
      }
      else {
        angle_1_5_2=spherical_angle_acc(v1, v5, v2);
        if (max(angle_1_5_0,angle_1_2_0)<angle_1_5_2) {
          index[0]=i_in;
          index[1]=j_in-1;
          index[2]=l_in;
        }
        else {
          ok=0;
        }
      }
    }
  }
  return ok;

}; /* get_index_acc */


/*-------------------------------------------------------------
  determine lower left corner
  void get_closest_index_acc(int ig, int jg, int lg, int ok)
  --------------------------------------------------------------*/
int get_closest_index_acc(const Grid_config *input_grid, const Grid_config *output_grid, int *index,
                          int i_in, int j_in, int l_in, int i_out, int j_out)
{

  double angle_11, angle_11a, angle_11b;
  double angle_22, angle_22a, angle_22b;
  double angle_33, angle_33a, angle_33b;
  double angle_44, angle_44a, angle_44b;
  double angle_1,  angle_1a,  angle_1b;
  double angle_2,  angle_2a,  angle_2b;
  double angle_3,  angle_3a,  angle_3b;
  double angle_4,  angle_4a,  angle_4b;
  int n4, n5, n6, n7, n8;
  double v0[3], v1[3], v2[3], v3[3], v4[3], v5[3], v6[3], v7[3], v8[3];

  int found = 0;
  int nx_in  = input_grid->nx;
  int ny_in  = input_grid->ny;
  int nxd    = nx_in + 2;
  int nx_out = output_grid->nx_fine;
  int ny_out = output_grid->ny_fine;
  int n0     = j_out*nx_out+i_out;
  int n1     = j_in*nxd+i_in;
  int n2     = j_in*nxd+i_in+1;
  int n3     = (j_in+1)*nxd+i_in;

  get_vector(input_grid, n1, v1);
  get_vector(input_grid, n2, v2);
  get_vector(input_grid, n3, v3);
  get_vector(output_grid, n0, v0);

  angle_1 =spherical_angle_acc(v1, v2, v3);
  angle_1a=spherical_angle_acc(v1, v2, v0);
  angle_1b=spherical_angle_acc(v1, v3, v0);

  if (max(angle_1a,angle_1b) <= angle_1) {
    if (i_in==nx_in && j_in==ny_in) {
      angle_11 =spherical_angle_acc(v2, v3, v1);
      angle_11a=spherical_angle_acc(v2, v1, v0);
      angle_11b=spherical_angle_acc(v2, v3, v0);
    }
    else {
      n4 = (j_in+1)*nxd+i_in+1;
      get_vector(input_grid, n4, v4);
      angle_11 =spherical_angle_acc(v4, v3, v2);
      angle_11a=spherical_angle_acc(v4, v2, v0);
      angle_11b=spherical_angle_acc(v4, v3, v0);
    }
    if (max(angle_11a,angle_11b)<=angle_11) {
      found = 1;
      index[0]=i_in;
      index[1]=j_in;
      index[2]=l_in;
    }
  }

  else {
    n4 = j_in*nxd+i_in-1;
    get_vector(input_grid, n4, v4);
    angle_2 =spherical_angle_acc(v1,v3,v4);
    angle_2a=angle_1b;
    angle_2b=spherical_angle_acc(v1,v4,v0);
    if (max(angle_2a,angle_2b)<=angle_2) {
      if (i_in==1 && j_in==ny_in) {
        angle_22 =spherical_angle_acc(v3, v1, v4);
        angle_22a=spherical_angle_acc(v3, v4, v0);
        angle_22b=spherical_angle_acc(v3, v1, v0);
      }
      else {
        n5 = (j_in+1)*nxd+i_in-1;
        n6 = j_in    *nxd+i_in-1;
        get_vector(input_grid, n5, v5);
        get_vector(input_grid, n6, v6);
        angle_22 =spherical_angle_acc(v5, v3, v6);
        angle_22a=spherical_angle_acc(v5, v6, v0);
        angle_22b=spherical_angle_acc(v5, v3, v0);
      }
      if (max(angle_22a,angle_22b)<=angle_22) {
        found=1;
        index[0]=i_in-1;
        index[1]=j_in;
        index[2]=l_in;
      }
    }
    else {
      n5 = j_in*nxd+i_in-1;
      n6 = (j_in-1)*nxd+i_in;
      get_vector(input_grid, n5, v5);
      get_vector(input_grid, n6, v6);
      angle_3 =spherical_angle_acc(v1, v5, v6);
      angle_3a=angle_2b;
      angle_3b=spherical_angle_acc(v1, v6, v0);
      if (max(angle_3a,angle_3b)<=angle_3 && i_in>1 && j_in>1) {
        n7 = (j_in-1)*nxd+i_in-1;
        get_vector(input_grid, n7, v7);
        angle_33 =spherical_angle_acc(v7, v6, v5);
        angle_33a=spherical_angle_acc(v7, v5, v0);
        angle_33b=spherical_angle_acc(v7, v6, v0);
        if (max(angle_33a,angle_33b)<=angle_33) {
          found=1;
          index[0]=i_in-1;
          index[1]=j_in-1;
          index[2]=l_in;
        }
      }
      else {
        angle_4 =spherical_angle_acc(v1, v6, v2);
        angle_4a=angle_3b;
        angle_4b=spherical_angle_acc(v1, v2, v0);
        if (max(angle_4a,angle_4b)<=angle_4) {
          if (i_in==nx_in && j_in==1) {
            angle_44 =spherical_angle_acc(v2, v1, v6);
            angle_44a=spherical_angle_acc(v2, v6, v0);
            angle_44b=spherical_angle_acc(v2, v1, v0);
          }
          else {
            n8 = (j_in-1)*nxd+i_in+1;
            get_vector(input_grid, n8, v8);
            angle_44 =spherical_angle_acc(v8, v2, v6);
            angle_44a=spherical_angle_acc(v8, v6, v0);
            angle_44b=spherical_angle_acc(v8, v2, v0);
          }
          if (max(angle_44a,angle_44b)<=angle_44) {
            found=1;
            index[0]=i_in;
            index[1]=j_in-1;
            index[2]=l_in;
          }
        }
      }
    }
  }

  return found;


}; /* get_closest_index_acc */


/*--------------------------------------------------------------------------

calculate normalized great circle distance between v1 and v2
double normalize_great_circle_distance_acc(v1, v2)
---------------------------------------------------------------------------*/
double normalize_great_circle_distance_acc(const double *v1, const double *v2)
{
  double dist;

  dist=(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
    /sqrt((v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])
	  *(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]));
  dist = sign(min(1.,fabs(dist)),dist);
  dist = acos(dist);
  return dist;

}; /* normalize_great_circle_distance_acc */

/*---------------------------------------------------------------------
  double dist2side_acc(v1, v2, point)
  calculate shortest normalized distance on sphere
  from point to straight line defined by v1 and v2
  ------------------------------------------------------------------*/
double dist2side_acc(const double *v1, const double *v2, const double *point)
{
  double angle, side;

  angle = spherical_angle_acc(v1, v2, point);
  side  = normalize_great_circle_distance_acc(v1, point);

  return (asin(sin(side)*sin(angle)));

};/* dist2side_acc */


int max_weight_index_acc( double *var, int nvar)
{

  int ind, i;

  ind = 0;

  for(i=1; i<nvar; i++) {
    if(var[i]>var[ind]) ind = i;
  }

  return ind;
}

/*------------------------------------------------------------------------------
  void do_latlon_coarsening_acc(var_latlon, ylat, nlon, nlat, nz,
                            var_latlon_crs, nlon_crs, nlat_crs,
                            finer_steps, misval, varmisval)

  calculate variable on coarser latlon grid
  by doubling spatial resolution and preserving volume means
  ---------------------------------------------------------------------------*/
void do_latlon_coarsening_acc(const double *var_latlon, const double *ylat, int nlon, int nlat, int nz,
			  double *var_latlon_crs, int finer_steps, int has_missing, double missvalue)
{

  double  *var_latlon_old, *ylat_old, *var_latlon_new;
  double  dlat;
  int     nlon_old, nlat_old, nlon_new, nlat_new, steps, i, j;
  int     nlon_crs, nlat_crs;

  nlon_crs=nlon/pow(2,finer_steps);
  nlat_crs=(nlat-1)/pow(2,finer_steps)+1;
  switch (finer_steps) {
  case 0:
    if (nlon_crs !=nlon || nlat_crs != nlat) mpp_error("bilinear_interp(do_latlon_coarsening_acc): grid dimensions don't match");
    for(i=0; i<nlon*nlat; i++) var_latlon_crs[i] = var_latlon[i];
    break;

  case 1:
    redu2x_acc(var_latlon, ylat, nlon, nlat, var_latlon_crs, nlon_crs, nlat_crs, has_missing, missvalue);
    break;

  default:
    nlon_new=nlon;
    nlat_new=nlat;
    for(steps=1; steps<=finer_steps; steps++) {
      nlon_old=nlon_new;
      nlat_old=nlat_new;
      var_latlon_old = (double *)malloc(nlon_old*nlat_old*sizeof(double));
      ylat_old       = (double *)malloc(nlat_old*sizeof(double));
      if (steps==1) {
        for(i=0; i<nlat; i++) ylat_old[i] = ylat[i];
        for(i=0; i<nlon_old*nlat_old; i++) var_latlon_old[i] = var_latlon[i];
      }
      else {
        dlat=M_PI/(nlat_old-1);
        ylat_old[0]=-0.5*M_PI;
        ylat_old[nlat_old-1]= 0.5*M_PI;
        for(j=1; j<nlat_old-1; j++) ylat_old[j] = ylat_old[0] + j*dlat;
        for(i=0; i<nlon_old*nlat_old*nz; i++) var_latlon_old[i] = var_latlon_new[i];
        free(var_latlon_new);
      }

      nlon_new=nlon_new/2;
      nlat_new=(nlat_new-1)/2+1;
      var_latlon_new = (double *)malloc(nlon_new*nlat_new*nz*sizeof(double));
      redu2x_acc(var_latlon_old, ylat_old, nlon_old, nlat_old, var_latlon_new, nlon_new, nlat_new, has_missing, missvalue);
      free(var_latlon_old);
      free(ylat_old);
    }
    for(i=0; i<nlon_new*nlat_new*nz; i++) var_latlon_crs[i] = var_latlon_new[i];
    free(var_latlon_new);
  }
}; /* do_latlon_coarsening_acc */

/*------------------------------------------------------------------------------
  void redu2x_acc(varfin, yfin, nxfin, nyfin, varcrs, nxcrs, nycrs)
  this routine is for reducing fvccm data by a factor of 2
  volume averaging for all data except at the poles
  original developer: S.-J. Lin
  ----------------------------------------------------------------------------*/
void redu2x_acc(const double *varfin, const double *yfin, int nxfin, int nyfin, double *varcrs,
                int nxcrs, int nycrs, int has_missing, double missvalue)
{
  double  *cosp, *acosp, *vartmp;
  int     i, j, n1;

#pragma acc enter data copyin( varcrs[:nxcrs*nycrs] )

  /*------------------------------------------------------------------
    calculate cosine of latitude
    trick in cosp needed to maintain a constant field
    ----------------------------------------------------------------*/
  cosp  = (double *)malloc(nyfin*sizeof(double));
  acosp = (double *)malloc(nyfin*sizeof(double));
  vartmp = (double *)malloc(nxcrs*nyfin*sizeof(double));

  cosp[0]=0.;
  cosp[nyfin-1]=0.;
#pragma acc enter data copyin(acosp[:nyfin], vartmp[:nxcrs*nyfin])

#pragma acc parallel loop present(cosp[:nyfin], acosp[:nyfin])
  for(int j=1; j<nyfin-1; j++) {
    cosp[j]  = cos(yfin[j]);
    acosp[j] = 1./(cosp[j]+0.5*(cosp[j-1]+cosp[j+1]));
  }

  /*----------------------------------------------------------------
    x-sweep
    ----------------------------------------------------------------*/
  if(has_missing) {
#pragma acc parallel loop present(varfin[:nxfin*nyfin], vartmp[:nxcrs*nyfin])
    for(int j=1; j<nyfin-1; j++) {
      int n1 = j*nxfin, n2 = j*nxcrs;
      if (varfin[n1+nxfin-1] == missvalue || varfin[n1] == missvalue || varfin[n1+1] == missvalue)
        vartmp[n2] = missvalue;
      else vartmp[n2] = 0.25*(varfin[n1+nxfin-1]+2.*varfin[n1]+varfin[n1+1]);
      for(int i=2; i<nxfin-1; i+=2) {
        if (varfin[n1+i-1] == missvalue || varfin[n1+i] == missvalue || varfin[n1+i+1] == missvalue)
          vartmp[n2+i/2] = missvalue;
        else vartmp[n2+i/2] = 0.25*(varfin[n1+i-1]+2.*varfin[n1+i]+varfin[n1+i+1]);
      }
    }
  }
  else {
#pragma acc parallel loop present(varfin[:nxfin*nyfin], vartmp[:nxcrs*nyfin])
    for(int j=1; j<nyfin-1; j++) {
      int n1 = j*nxfin, n2 = j*nxcrs;
      vartmp[n2] = 0.25*(varfin[n1+nxfin-1]+2.*varfin[n1]+varfin[n1+1]);
      for(int i=2; i<nxfin-1; i+=2) {
        vartmp[n2+i/2] = 0.25*(varfin[n1+i-1]+2.*varfin[n1+i]+varfin[n1+i+1]);
      }
    }
  }
  /*---------------------------------------------------------------------
    poles:
    this code segment works for both the scalar and vector fields.
    Winds at poles are wave-1; the follwoing is quick & dirty yet the correct way
    The skipping method. A more rigorous way is to
    recompute the wave-1 components for the coarser grid.
    --------------------------------------------------------------------*/
#pragma acc parallel loop present(varcrs[:nxcrs*nycrs], varfin[:nxfin*nyfin])
  for(int i=0; i<nxcrs; i++) {
    int i2 = i*2;
    varcrs[i] = varfin[i2];
    varcrs[(nycrs-1)*nxcrs+i] = varfin[(nyfin-1)*nxfin+i2];
  }
  /*----------------------------------------------------------------
    y-sweep
    ----------------------------------------------------------------*/
  if (has_missing) {
#pragma acc parallel loop collapse(2) present(vartmp[:nxcrs*nyfin], cosp[:nyfin])
    for(j=1; j<nyfin-1; j++) {
      for(i=0; i<nxcrs; i++) {
        if (vartmp[n1] /= missvalue) vartmp[j*nxcrs+i] *= cosp[j];
      }
    }
#pragma acc parallel loop collapse(2) present(vartmp[:nxcrs*nyfin], varcrs[:nxcrs*nycrs], acosp[:nyfin])
    for(j=2; j<nyfin-2; j+=2) {
      for(i=0; i<nxcrs; i++) {
        if (vartmp[i+j*nxcrs] == missvalue || vartmp[i+(j-1)*nxcrs] == missvalue
            || vartmp[i+(j+1)*nxcrs] == missvalue )
          varcrs[i+j*nxcrs/2] = missvalue;
        else
          varcrs[i+j*nxcrs/2] = acosp[j]*(vartmp[i+j*nxcrs] + 0.5*(vartmp[i+(j-1)*nxcrs]+ vartmp[i+(j+1)*nxcrs]));
      }
    }
  }
  else {
#pragma acc parallel loop collapse(2) present(vartmp[:nxcrs*nyfin], cosp[:nyfin])
    for(j=1; j<nyfin-1; j++) {
      for(i=0; i<nxcrs; i++) {
        vartmp[j*nxcrs+i] *= cosp[j];
      }
    }
#pragma acc parallel loop collapse(2) present(vartmp[:nxcrs*nyfin], acosp[:nyfin], varcrs[:nxcrs*nycrs])
    for(j=2; j<nyfin-2; j+=2) {
      for(i=0; i<nxcrs; i++) {
        varcrs[i+j*nxcrs/2] = acosp[j]*(vartmp[i+j*nxcrs] + 0.5*(vartmp[i+(j-1)*nxcrs]+ vartmp[i+(j+1)*nxcrs]));
      }
    }
  }

#pragma acc exit data delete(vartmp[:nxcrs*nyfin], cosp[:nyfin], acosp[:nyfin])

  free(cosp);
  free(acosp);
  free(vartmp);

}; /*redu2x_acc*/

void read_remap_file(const Grid_config *output_grid, const Interp_config_acc *interp_acc)
{
  // check the size of the grid matching the size in remapping file
  printf("NOTE: reading index and weight for bilinear interpolation from file.\n");

  int nlon_output_cells = output_grid->nx_fine;
  int nlat_output_cells = output_grid->ny_fine;

  int fid = mpp_open(interp_acc->remap_file, MPP_READ);
  int nx2 = mpp_get_dimlen(fid, "nlon");
  int ny2 = mpp_get_dimlen(fid, "nlat");

  printf("grid size is nx=%d, ny=%d, remap file size is nx=%d, ny=%d.\n",
         nlon_output_cells, nlat_output_cells, nx2, ny2);
  if(nx2 != nlon_output_cells || ny2 != nlat_output_cells )
    mpp_error("bilinear_interp: size mismatch between grid size and remap file size");

  int vid = mpp_get_varid(fid, "index");
  mpp_get_var_value(fid, vid, interp_acc->index);

  vid = mpp_get_varid(fid, "weight");
  mpp_get_var_value(fid, vid, interp_acc->weight);

  mpp_close(fid);

#pragma acc enter data copyin(interp_acc[:1])
#pragma acc enter data copyin(interp_acc->index[:3*nlon_output_cells*nlat_output_cells], \
                              interp_acc->weight[:4*nlon_output_cells*nlat_output_cells])

}


void get_interp_index_test( const int input_ntiles, const int iter, const Grid_config *input_grid,
                            const Grid_config *output_grid, const double dlon_in, const double dlat_in,
                            const double lonbegin, const double latbegin, Interp_config_acc *interp_acc)
{

  int total_found = 0;
  int nlon_output_cells = output_grid->nx_fine;
  int nlat_output_cells = output_grid->ny_fine;
  int ncells_output = nlon_output_cells * nlat_output_cells;

  int nlon_input_cells = input_grid->nx;
  int nlat_input_cells = input_grid->ny;
  int nxd = nlon_input_cells + 2;
  int nyd = nlat_input_cells + 2;

  double dlon = dlon_in/nlon_output_cells;
  double dlat = dlat_in/(nlat_output_cells-1);

  int *i2_min=NULL; i2_min = acc_malloc(nxd*nyd*sizeof(int));
  int *i2_max=NULL; i2_max = acc_malloc(nxd*nyd*sizeof(int));
  int *j2_min=NULL; j2_min = acc_malloc(nxd*nyd*sizeof(int));
  int *j2_max=NULL; j2_max = acc_malloc(nxd*nyd*sizeof(int));
  int *i1_min=NULL; i1_min = acc_malloc(ncells_output*sizeof(int));
  int *i1_max=NULL; i1_max = acc_malloc(ncells_output*sizeof(int));
  int *j1_min=NULL; j1_min = acc_malloc(ncells_output*sizeof(int));
  int *j1_max=NULL; j1_max = acc_malloc(ncells_output*sizeof(int));
  int *found=NULL; found = acc_malloc(ncells_output*sizeof(int));
  double *shortest=NULL; shortest = acc_malloc(ncells_output*sizeof(double));

  float time[input_ntiles];
  for(int i=0 ; i<input_ntiles ; i++) time[i] = 0.0;

  for(int iiter=1 ; iiter<iter ; iiter++ ) {

#pragma acc parallel loop deviceptr(shortest)
    for(int i=0; i<ncells_output; i++) {
      shortest[i] = M_PI + M_PI;
    }

    for(int itile=0; itile<input_ntiles; itile++) {

      clock_t time_start, time_end;
      time_start = clock();

      int i2_start=ncells_output, i2_end=0, j2_start=ncells_output, j2_end=0;

#pragma acc parallel loop collapse(2) reduction(min: i2_start) reduction(min: j2_start) \
                                      reduction(max: i2_end) reduction(max: j2_end)     \
                                      present(input_grid[itile].xt[:nxd*nyd], \
                                              input_grid[itile].yt[:nxd*nyd], \
                                              input_grid[itile].zt[:nxd*nyd]) \
                                      copyin(input_grid[itile].latt[:nxd*nyd],\
                                             input_grid[itile].lont[:nxd*nyd])\
                                      copyout(i2_start, i2_end, j2_start, j2_end) \
                                      deviceptr(i2_min, i2_max, j2_min, j2_max)
      for(int jc=1; jc<=nlat_input_cells; jc++) {
        for(int ic=1; ic<=nlon_input_cells; ic++) {
          /*------------------------------------------------------
            guess bounding latlon output cells for given input cubed sphere cell
            ------------------------------------------------------*/
          int n1 = jc*nxd+ic;
          int n2 = (jc+1)*nxd+ic+1;

          int outgrid_j_min, outgrid_j_max;
          int outgrid_i_min, outgrid_i_max;
          double input_latt = (input_grid[itile].latt[n1] - latbegin);
          double input_lont = (input_grid[itile].lont[n1] - lonbegin);
          double dcub, v1[3], v2[3];

          get_vector(input_grid+itile, n1, v1);
          get_vector(input_grid+itile, n2, v2);
          dcub=iter*normalize_great_circle_distance_acc(v1, v2);

          // get approximate output grid indices in the vicinity
          outgrid_j_min=max(   1,  floor((input_latt-dcub)/dlat)-iiter+1 );
          outgrid_j_max=min(nlat_output_cells, ceil((input_latt+dcub)/dlat)+iiter-1 );

          outgrid_i_min=max(   1,  floor((input_lont-dcub)/dlon-iiter+1));
          outgrid_i_max=min(nlon_output_cells,ceil((input_lont+dcub)/dlon+iiter-1));

          if(outgrid_j_min==1 ) {
            outgrid_i_min = 1;
            outgrid_i_max = nlon_output_cells;
          }
          if(outgrid_j_max==nlat_output_cells ) {
            outgrid_i_min = 1;
            outgrid_i_max = nlon_output_cells;
          }

          i2_min[n1] = outgrid_i_min-1;
          i2_max[n1] = outgrid_i_max;
          j2_min[n1] = outgrid_j_min-1;
          j2_max[n1] = outgrid_j_max;

          i2_start = min(i2_start, outgrid_i_min-1);
          j2_start = min(j2_start, outgrid_j_min-1);
          i2_end = max(i2_end, outgrid_i_max);
          j2_end = max(j2_end, outgrid_j_max);

        } // input ic
      } // input jc


#pragma acc parallel loop collapse(2) deviceptr(i1_min, i1_max, j1_min, j1_max, \
                                                i2_min, i2_max, j2_min, j2_max, shortest, found)\
                                      copyin( i2_start, j2_start, i2_end, j2_end )
      for( int j2=j2_start ; j2<j2_end ; j2++) {
        for(int i2=i2_start ; i2<i2_end ; i2++) {
          int n2 = j2*nlon_output_cells + i2;
          int i1_min_tmp = nxd ; int i1_max_tmp = -10;
          int j1_min_tmp = nyd ; int j1_max_tmp = -10;
#pragma acc loop collapse(2) reduction(min:i1_min_tmp) reduction(min:j1_min_tmp) \
                             reduction(max:i1_max_tmp) reduction(max:j1_max_tmp)
          for( int j1=1 ; j1<=nlat_input_cells ; j1++) {
            for( int i1=1 ; i1<=nlon_input_cells ; i1++) {
              int n1 = j1*nxd + i1;
              if( j2 < j2_min[n1] ) continue;
              if( i2 < i2_min[n1] ) continue;
              if( j2 > j2_max[n1] ) continue;
              if( i2 > i2_max[n1] ) continue;
              i1_min_tmp = min(i1_min_tmp, i1);
              i1_max_tmp = max(i1_max_tmp, i1);
              j1_min_tmp = min(j1_min_tmp, j1);
              j1_max_tmp = max(j1_max_tmp, j1);
            }
          }
          i1_min[n2] = i1_min_tmp;
          i1_max[n2] = i1_max_tmp;
          j1_min[n2] = j1_min_tmp;
          j1_max[n2] = j1_max_tmp;
        }
      }


#pragma acc parallel loop collapse(2) present(input_grid[itile].xt[:nxd*nyd],   \
                                              input_grid[itile].yt[:nxd*nyd], \
                                              input_grid[itile].zt[:nxd*nyd], \
                                              output_grid->xt[:ncells_output], \
                                              output_grid->yt[:ncells_output], \
                                              output_grid->zt[:ncells_output], \
                                              interp_acc->index[:3*ncells_output]) \
                                       deviceptr(shortest, found, i1_min, j1_min, i1_max, j1_max)\
                                       copyin(i2_start, j2_start, i2_end, j2_end)
      for( int j2=j2_start ; j2<j2_end; j2++ ) {
        for( int i2=i2_start ; i2<i2_end ; i2++ ) {
          int n2 = j2*nlon_output_cells + i2;
          int i1_start = i1_min[n2];
          int j1_start = j1_min[n2];
          int i1_end = i1_max[n2];
          int j1_end = j1_max[n2];
          double distance, v1[3], v2[3];
          get_vector(output_grid, n2, v2);
#pragma acc loop seq
          for( int j1=j1_start ; j1<=j1_end ; j1++) {
#pragma acc loop seq
            for( int i1=i1_start ; i1<=i1_end ; i1++) {
              int n1 = j1*nxd + i1;
              get_vector(input_grid+itile, n1, v1);
              distance = normalize_great_circle_distance_acc(v1, v2);
              if( distance > shortest[n2]) continue;
              if( get_closest_index_acc(input_grid+itile, output_grid, interp_acc->index+3*n2,
                                        i1, j1, itile, i2, j2) ) {
                found[n2] = 1;
                shortest[n2] = distance;
              }
            }
          }
        }
      }

      time_end = clock();
      time[itile] = time[itile] + (float)(time_end-time_start)/CLOCKS_PER_SEC;

    } // itile

    total_found = 0;
#pragma acc parallel loop deviceptr(found) copy(total_found)
    for(int n2=0; n2<ncells_output ; n2++) {
      total_found = total_found + found[n2];
    }
    if(total_found == ncells_output) break;

  } //iiter

for(int itile=0 ; itile<input_ntiles ; itile++) printf("ITILE %d %f\n", itile, time[itile]);

  acc_free(i2_min);
  acc_free(i2_max);
  acc_free(j2_min);
  acc_free(j2_max);
  acc_free(i1_min);
  acc_free(i1_max);
  acc_free(j1_min);
  acc_free(j1_max);
  acc_free(found);
  acc_free(shortest);

  if(total_found != ncells_output) {
    printf("did not find all neighboring points for each output cell %d %d", total_found, ncells_output);
    exit(1);
  }

}


void get_interp_weights(const int input_ntiles, const int output_ntiles, const Grid_config *output_grid,
                        const Grid_config *input_grid, Interp_config_acc *interp_acc)
{

  int nlon_output_cells = output_grid->nx_fine;
  int nlat_output_cells = output_grid->ny_fine;
  int nlon_input_cells = input_grid->nx;
  int nlat_input_cells = input_grid->ny;
  int nxd = nlon_input_cells + 2;
  int nyd = nlat_input_cells + 2;

#pragma acc parallel loop collapse(2) present(interp_acc[:1], input_grid[:input_ntiles], output_grid[:output_ntiles])
  for(int j=0; j<nlat_output_cells; j++) {
    for(int i=0; i<nlon_output_cells; i++) {

      int n0 = j*nlon_output_cells + i;
      int m0 = 3*n0;
      int m1 = 4*n0;
      int ic = interp_acc->index[m0];
      int jc = interp_acc->index[m0+1];
      int itile = interp_acc->index[m0+2];
      int n1, n2, n3, n4;
      double v0[3], v1[3], v2[3], v3[3], v4[3];
      double dist1, dist2, dist3, dist4;
      double sum;

      get_vector(output_grid, n0, v0);

      //The section of the code that double checks that found == 0 for
      //all output cells been removed because it looked very buggy.

      /*------------------------------------------------------------
        calculate shortest distance to each side of rectangle
        formed by cubed sphere cell centers
        special corner treatment
        ------------------------------------------------------------*/

      if (ic==nlon_input_cells && jc==nlat_input_cells) {
        /*------------------------------------------------------------
          calculate weights for bilinear interpolation near corner
          ------------------------------------------------------------*/
        n1 = jc*nxd+ic;     get_vector(input_grid+itile, n1, v1);
        n2 = jc*nxd+ic+1;   get_vector(input_grid+itile, n2, v2);
        n3 = (jc+1)*nxd+ic; get_vector(input_grid+itile, n3, v3);

        interp_acc->weight[m1]  =dist2side_acc(v2, v3, v0);  //ic,   jc    weight
        interp_acc->weight[m1+1]=dist2side_acc(v2, v1, v0);  //ic,   jc+1  weight
        interp_acc->weight[m1+2]=0.;                         //ic+1, jc+1  weight
        interp_acc->weight[m1+3]=dist2side_acc(v1, v3, v0);  //ic+1, jc    weight
        sum=interp_acc->weight[m1] + interp_acc->weight[m1+1] + interp_acc->weight[m1+2] + interp_acc->weight[m1+3];
        interp_acc->weight[m1]  /=sum;
        interp_acc->weight[m1+1]/=sum;
        interp_acc->weight[m1+2]/=sum;
        interp_acc->weight[m1+3]/=sum;
        continue;
      }
      if (ic==0 && jc==nlat_input_cells) {
        // calculate weights for bilinear interpolation near corner
        n1 = jc*nxd+ic;       get_vector(input_grid+itile, n1, v1);
        n2 = jc*nxd+ic+1;     get_vector(input_grid+itile, n2, v2);
        n3 = (jc+1)*nxd+ic+1; get_vector(input_grid+itile, n3, v3);

        interp_acc->weight[m1]=dist2side_acc(v3, v2, v0);   // ic,   jc    weight
        interp_acc->weight[m1+1]=0.;                        // ic,   jc+1  weight
        interp_acc->weight[m1+2]=dist2side_acc(v2, v1, v0); // ic+1, jc+1  weight
        interp_acc->weight[m1+3]=dist2side_acc(v3, v1, v0); // ic+1, jc    weight
        sum=interp_acc->weight[m1] + interp_acc->weight[m1+1] + interp_acc->weight[m1+2] + interp_acc->weight[m1+3];
        interp_acc->weight[m1]  /=sum;
        interp_acc->weight[m1+1]/=sum;
        interp_acc->weight[m1+2]/=sum;
        interp_acc->weight[m1+3]/=sum;
        continue;
      }
      if (jc==0 && ic==nlon_input_cells) {
        // calculate weights for bilinear interpolation near corner
        n1 = jc*nxd+ic;       get_vector(input_grid+itile, n1, v1);
        n2 = (jc+1)*nxd+ic;   get_vector(input_grid+itile, n2, v2);
        n3 = (jc+1)*nxd+ic+1; get_vector(input_grid+itile, n3, v3);

        interp_acc->weight[m1]  =dist2side_acc(v2, v3, v0); // ic,   jc    weight
        interp_acc->weight[m1+1]=dist2side_acc(v1, v3, v0); // ic,   jc+1  weight
        interp_acc->weight[m1+2]=dist2side_acc(v1, v2, v0); // ic+1, jc+1  weight
        interp_acc->weight[m1+3]=0.;                        // ic+1, jc    weight
        sum=interp_acc->weight[m1] + interp_acc->weight[m1+1] + interp_acc->weight[m1+2] + interp_acc->weight[m1+3];
        interp_acc->weight[m1]  /=sum;
        interp_acc->weight[m1+1]/=sum;
        interp_acc->weight[m1+2]/=sum;
        interp_acc->weight[m1+3]/=sum;
        continue;
      }

      // calculate weights for bilinear interpolation if no corner
      n1 = jc*nxd+ic;       get_vector(input_grid+itile, n1, v1);
      n2 = jc*nxd+ic+1;     get_vector(input_grid+itile, n2, v2);
      n3 = (jc+1)*nxd+ic;   get_vector(input_grid+itile, n3, v3);
      n4 = (jc+1)*nxd+ic+1; get_vector(input_grid+itile, n4, v4);

      dist1=dist2side_acc(v1, v3, v0);
      dist2=dist2side_acc(v3, v4, v0);
      dist3=dist2side_acc(v4, v2, v0);
      dist4=dist2side_acc(v2, v1, v0);

      interp_acc->weight[m1]  =dist2*dist3; // ic,   jc    weight
      interp_acc->weight[m1+1]=dist3*dist4; // ic,   jc+1  weight
      interp_acc->weight[m1+2]=dist4*dist1; // ic+1, jc+1  weight
      interp_acc->weight[m1+3]=dist1*dist2; // ic+1, jc    weight

      sum=interp_acc->weight[m1] + interp_acc->weight[m1+1] + interp_acc->weight[m1+2] + interp_acc->weight[m1+3];
      interp_acc->weight[m1]  /=sum;
      interp_acc->weight[m1+1]/=sum;
      interp_acc->weight[m1+2]/=sum;
      interp_acc->weight[m1+3]/=sum;
    }
  }

}


void get_vector(const Grid_config *grid, const int icell, double *vector)
{

  vector[0] = grid->xt[icell];
  vector[1] = grid->yt[icell];
  vector[2] = grid->zt[icell];

}

void write_bilinear_remap_file(const int nlon_output_cells, const int nlat_output_cells, Interp_config_acc *interp_acc)
{

  int dims[3];
  int fid = mpp_open( interp_acc->remap_file, MPP_WRITE);
  int dim_nlon = mpp_def_dim(fid, "nlon", nlon_output_cells);
  int dim_nlat = mpp_def_dim(fid, "nlat", nlat_output_cells);
  int dim_three = mpp_def_dim(fid, "three", 3);
  int dim_four  = mpp_def_dim(fid, "four", 4);

  dims[0] = dim_three; dims[1] = dim_nlat; dims[2] = dim_nlon;
  int fld_index = mpp_def_var(fid, "index", NC_INT, 3, dims, 0);

  dims[0] = dim_four; dims[1] = dim_nlat; dims[2] = dim_nlon;
  int fld_weight = mpp_def_var(fid, "weight", NC_DOUBLE, 3, dims, 0);

  mpp_end_def(fid);

  mpp_put_var_value(fid, fld_index, interp_acc->index);
  mpp_put_var_value(fid, fld_weight, interp_acc->weight);
  mpp_close(fid);

}

void write_bilinear_interp_in_not_zhi_format(const int ncells, Interp_config_acc *interp_acc)
{

  int fid = mpp_open( "test_remap.nc", MPP_WRITE);
  int dim_ncells = mpp_def_dim(fid, "ncells", ncells);
  int fld_index1 = mpp_def_var(fid, "i_index", NC_INT, 1, &dim_ncells, 0);
  int fld_index2 = mpp_def_var(fid, "j_index", NC_INT, 1, &dim_ncells, 0);
  int fld_index3 = mpp_def_var(fid, "tile_index", NC_INT, 1, &dim_ncells, 0);

  mpp_end_def(fid);

  int *array1 ; array1 = (int *)malloc(ncells*sizeof(int));
  int *array2 ; array2 = (int *)malloc(ncells*sizeof(int));
  int *array3 ; array3 = (int *)malloc(ncells*sizeof(int));

  for(int i=0 ; i<ncells ; i++) {
    array1[i] = interp_acc->index[3*i];
    array2[i] = interp_acc->index[3*i+1];
    array3[i] = interp_acc->index[3*i+2];
  }

  mpp_put_var_value(fid, fld_index1, array1);
  mpp_put_var_value(fid, fld_index2, array2);
  mpp_put_var_value(fid, fld_index3, array3);
  mpp_close(fid);

  free(array1);
  free(array2);
  free(array3);

}
