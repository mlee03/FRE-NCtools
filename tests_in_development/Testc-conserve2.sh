#!/usr/bin/bash

source run_tools.sh
echo "Testc second order conservative remapping from C384 to lonlat 1152x720"

dir_in=/home/Mikyung.Lee/FRE-NCTools/TESTS_INPUT/Testc-input
dir_out=$PWD/tests_in_development/Testc-conserve2-output

if [ -d $dir_out ] ; then rm -rf $dir_out ; fi
mkdir -p $dir_out && cd $dir_out

data_file=00010101.land_month_cmip.tile1

scalar_field=\
"cSoil,nep,fLuc,cLand,fAnthDisturb,fFireNat,fProductDecomp,netAtmosLandCO2Flux,
treeFracBdlDcd,treeFracBdlEvg,treeFracNdlDcd,treeFracNdlEvg,vegFrac,vegHeight,
snc,snd,snw,snm,hfdsn,sbl,mrsos,mrso,mrfso,mrros,mrro,prveg,evspsblveg,evspsblsoi,
tran,mrlsl,tsl,treeFrac,grassFrac,cropFrac,pastureFrac,residualFrac,cVeg,cProduct,
lai,gpp,ra,npp,rh,fFire,fGrazing,fHarvest,nbp"

input_mosaic=$(run_generate_cubed_sphere 384)
output_mosaic=$(run_generate_lonlat 1152 720)

run_fregrid_lonlat_conserve2 $input_mosaic \
                             $output_mosaic \
                             $dir_in \
                             $dir_out \
                             $data_file \
                             $scalar_field

run_fregrid_lonlat_conserve2 $input_mosaic \
                             $output_mosaic \
                             $dir_in \
                             $dir_out \
                             $data_file \
                             $scalar_field
