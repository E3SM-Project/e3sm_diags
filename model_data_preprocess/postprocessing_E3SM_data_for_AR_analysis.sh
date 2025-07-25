#!/bin/bash

#===============================================================================
# E3SM Atmospheric River (AR) Detection and Analysis Script
#===============================================================================
#
# Purpose: Post-process E3SM 6-hourly (h2) instantaneous output to detect
#          atmospheric rivers while filtering out tropical cyclones
#
# Requirements:
#   - TempestRemap and TempestExtremes (available in e3sm-unified v1.5.0+)
#   - E3SM h2 (6-hourly) output files with TUQ/TVQ variables
#
# Usage: bash postprocessing_E3SM_data_for_AR_analysis.sh
#
# Author: Jill Zhang,
# Based on Paul Ullrich et al. 2021: TempestExtremes v2.1: a community framework for feature detection, tracking, and analysis in large datasets, and
# Reed et al. 2023: Evaluating the Simulation of CONUS Precipitation by Storm Type in E3SM.
# Date: July 2025
# Version: 1.0
#===============================================================================

#===============================================================================
# CONFIGURATION SECTION
#===============================================================================

# Load E3SM unified environment
# Uncomment the appropriate line for your system:
# source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh  # Perlmutter
source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh            # LCRC Chrysalis

# Simulation parameters
start="1850"                                    # Start year
end="1850"                                      # End year
caseid="v3.LR.historical_0051"                 # Case ID
atm_name="eam"                                  # Atmosphere component name (eam/cam)

# Grid configuration
res=30                                          # Grid resolution (30 for ne30, 120 for ne120)
pg2=true                                        # Use pg2 grids (true for v2/v3, false for v1)

# Directory paths
drc_in="/lcrc/group/e3sm/ac.forsyth2/E3SMv3/${caseid}/${caseid}/archive/atm/hist"
result_dir="/lcrc/group/e3sm/ac.zhang40/tests/ar_analysis/ar-${caseid}.${start}_${end}/"

# Derived variables
file_name="${caseid}_${start}_${end}"

echo "=============================================================================="
echo "E3SM Atmospheric River Detection Script"
echo "=============================================================================="
echo "Case ID: ${caseid}"
echo "Time period: ${start} - ${end}"
echo "Grid resolution: ne${res}"
echo "Input directory: ${drc_in}"
echo "Output directory: ${result_dir}"
echo "=============================================================================="

#===============================================================================
# SETUP AND INITIALIZATION
#===============================================================================

# Create output directory structure
echo "Setting up output directories..."
mkdir -p "${result_dir}"
mkdir -p "${result_dir}AR_mask_nofilt"
mkdir -p "${result_dir}AR_mask_filt"
mkdir -p "${result_dir}TVQ_PRECT_AR_mask"
mkdir -p "${result_dir}TC_masks"
mkdir -p "${result_dir}ETC_masks"

#===============================================================================
# GRID GENERATION
#===============================================================================

echo "Generating grid and connectivity files..."

# Generate appropriate mesh files based on grid type
if ${pg2}; then
    echo "  Using pg2 grids (E3SM v2/v3)..."
    GenerateCSMesh \
        --res ${res} \
        --alt \
        --file "${result_dir}outCSMeshne${res}.g"

    GenerateVolumetricMesh \
        --in "${result_dir}outCSMeshne${res}.g" \
        --out "${result_dir}outCSne${res}.g" \
        --np 2 \
        --uniform

    out_type="FV"
else
    echo "  Using np4 grids (E3SM v1)..."
    GenerateCSMesh \
        --res ${res} \
        --alt \
        --file "${result_dir}outCSne${res}.g"

    out_type="CGLL"
fi

# Generate connectivity file
echo "  Grid type: ${out_type}"
GenerateConnectivityFile \
    --in_mesh "${result_dir}outCSne${res}.g" \
    --out_type ${out_type} \
    --out_connect "${result_dir}connect_CSne${res}_v2.dat"

#===============================================================================
# FILE LIST GENERATION
#===============================================================================

echo "Generating input file lists..."

# Clean up existing file lists
rm -f "${result_dir}inputfile_${file_name}.txt"
rm -f "${result_dir}ar_nofilt_files_out.txt"
rm -f "${result_dir}ar_filt_files_out.txt"
rm -f "${result_dir}TVQ_PRECT_AR_mask_files_in.txt"
rm -f "${result_dir}TVQ_PRECT_AR_mask_files_out.txt"
rm -f "${result_dir}tc_mask_files_out.txt"
rm -f "${result_dir}etc_mask_files_out.txt"

# Process input files and create systematic output filenames
file_count=0
echo "  Searching pattern: ${drc_in}/${caseid}.${atm_name}.h2.*{${start}..${end}}*.nc"

for f in $(eval echo "${drc_in}/${caseid}.${atm_name}.h2.*{${start}..${end}}*.nc"); do
    if [ -f "$f" ]; then
        g=$(basename "$f")

        date_part="${g#${caseid}.${atm_name}.h2.}"
        date_part="${date_part%.nc}"

        ar_nofilt_file="${result_dir}AR_mask_nofilt/${caseid}.${atm_name}.h2.${date_part}.AR_mask_nofilt.nc"
        ar_filt_file="${result_dir}AR_mask_filt/${caseid}.${atm_name}.h2.${date_part}.AR_mask_filt.nc"
        tvq_prect_ar_file="${result_dir}TVQ_PRECT_AR_mask/${caseid}.${atm_name}.h2.${date_part}.TVQ_PRECT_AR_mask.nc"
        tc_mask_file="${result_dir}TC_masks/${caseid}.${atm_name}.h2.${date_part}.TC_mask.nc"
        etc_mask_file="${result_dir}ETC_masks/${caseid}.${atm_name}.h2.${date_part}.ETC_mask.nc"
        echo "$f" >> "${result_dir}inputfile_${file_name}.txt"
        echo "${ar_nofilt_file}" >> "${result_dir}ar_nofilt_files_out.txt"
        echo "${ar_filt_file}" >> "${result_dir}ar_filt_files_out.txt"
        echo "${ar_filt_file};$f" >> "${result_dir}TVQ_PRECT_AR_mask_files_in.txt"
        echo "${tvq_prect_ar_file}" >> "${result_dir}TVQ_PRECT_AR_mask_files_out.txt"
        echo "${tc_mask_file}" >> "${result_dir}tc_mask_files_out.txt"
        echo "${etc_mask_file}" >> "${result_dir}etc_mask_files_out.txt"

        ((file_count++))
    fi
done

echo "  Found ${file_count} input files"

# Change to results directory
cd "${result_dir}"

#===============================================================================
# STEP 1: TROPICAL CYCLONE DETECTION
#===============================================================================

echo "Step 1: Detecting tropical cyclones for filtering..."
echo "  Method: Sea level pressure minima with temperature criteria"
echo "  SLP decrease: 300 Pa within 4° radius"
echo "  Temperature decrease: 0.6 K at 200/500 hPa levels"

# Resolution-dependent parameters for closed contour command
if [ ${res} -eq 120 ]; then
    echo "  Using ne120 parameters..."
    DetectNodes \
        --in_connect "${result_dir}connect_CSne${res}_v2.dat" \
        --closedcontourcmd "PSL,300.0,4.0,0;_AVG(T200,T500),-0.6,4,0.30" \
        --mergedist 6.0 \
        --searchbymin PSL \
        --outputcmd "PSL,min,0;_VECMAG(UBOT,VBOT),max,2" \
        --timestride 1 \
        --in_data_list "${result_dir}inputfile_${file_name}.txt" \
        --out "${result_dir}out.dat"
elif [ ${res} -eq 30 ]; then
    echo "  Using ne30 parameters..."
    DetectNodes \
        --in_connect "${result_dir}connect_CSne${res}_v2.dat" \
        --closedcontourcmd "PSL,300.0,4.0,0;_AVG(T200,T500),-0.6,4,1.0" \
        --mergedist 6.0 \
        --searchbymin PSL \
        --outputcmd "PSL,min,0;_VECMAG(UBOT,VBOT),max,2" \
        --timestride 1 \
        --in_data_list "${result_dir}inputfile_${file_name}.txt" \
        --out "${result_dir}out.dat"
else
    echo "  ERROR: Grid resolution ${res} not supported!"
    echo "  Supported resolutions: 30 (ne30), 120 (ne120)"
    exit 1
fi

echo "  Concatenating TC detection files..."
cat "${result_dir}"out.dat0*.dat > "${result_dir}cyclones_${file_name}.txt"

echo "  Creating TC tracks..."
echo "  Max distance: 6.0°, Min duration: 6 time steps, Max gap: 1 time step"
echo "  Wind speed: >= 17.5 m/s, Latitude: -40° to 40°"

StitchNodes \
    --in_fmt "lon,lat,slp,wind" \
    --in_connect "${result_dir}connect_CSne${res}_v2.dat" \
    --range 6.0 \
    --mintime 6 \
    --maxgap 1 \
    --in "${result_dir}cyclones_${file_name}.txt" \
    --out "${result_dir}cyclones_stitch_${file_name}.dat" \
    --threshold "wind,>=,17.5,6;lat,<=,40.0,6;lat,>=,-40.0,6"

echo "  Calculating TC radial wind profiles and circulation size..."
echo "  Radial bins: 0.5° for ne30, 0.125° for ne120"

# Set radial bin size based on resolution
if [ ${res} -eq 120 ]; then
    radial_bin="0.125"
else
    radial_bin="0.5"
fi

NodeFileEditor \
    --in_connect "${result_dir}connect_CSne${res}_v2.dat" \
    --in_data_list "${result_dir}inputfile_${file_name}.txt" \
    --in_nodefile "${result_dir}cyclones_stitch_${file_name}.dat" \
    --out_nodefile "${result_dir}cyclones_radprof_${file_name}.dat" \
    --in_fmt "lon,lat,slp,wind" \
    --calculate "rprof=radial_wind_profile(U10,V10,159,${radial_bin});rsize=lastwhere(rprof,>,8)" \
    --out_fmt "lon,lat,rsize,rprof"

echo "  TC radial wind profile calculation completed"

# Create TC mask files for precipitation identification
echo "  Creating TC mask files..."
echo "  Binary mask: 1 within TC circulation radius, 0 elsewhere"


NodeFileFilter \
    --in_connect "${result_dir}connect_CSne${res}_v2.dat" \
    --in_nodefile "${result_dir}cyclones_radprof_${file_name}.dat" \
    --in_fmt "lon,lat,rsize,rprof" \
    --in_data_list "${result_dir}inputfile_${file_name}.txt" \
    --out_data_list "${result_dir}tc_mask_files_out.txt" \
    --bydist "rsize" \
    --maskvar "storm_mask"

echo "  TC mask creation completed"

#===============================================================================
# STEP 2: EXTRATROPICAL CYCLONE DETECTION AND TRACKING
#===============================================================================

echo "Step 2: Detecting extratropical cyclones..."
echo "  Method: Sea level pressure minima with closed contour criteria"
echo "  SLP decrease: 200 Pa within 6° radius"
echo "  Warm-core check: T200/T500 average to eliminate TCs"
echo "  Merge distance: 6°"

# Resolution-dependent parameters for ETC warm-core check
if [ ${res} -eq 120 ]; then
    echo "  Using ne120 parameters..."
    DetectNodes \
        --in_data_list "${result_dir}inputfile_${file_name}.txt" \
        --out "${result_dir}etc_out.dat" \
        --timefilter "6hr" \
        --in_connect "${result_dir}connect_CSne${res}_v2.dat" \
        --closedcontourcmd "PSL,200.0,6.0,0" \
        --noclosedcontourcmd "_AVG(T200,T500),-0.6,4,0.30" \
        --mergedist 6.0 \
        --searchbymin PSL \
        --outputcmd "PSL,min,0"
elif [ ${res} -eq 30 ]; then
    echo "  Using ne30 parameters..."
    DetectNodes \
        --in_data_list "${result_dir}inputfile_${file_name}.txt" \
        --out "${result_dir}etc_out.dat" \
        --timefilter "6hr" \
        --in_connect "${result_dir}connect_CSne${res}_v2.dat" \
        --closedcontourcmd "PSL,200.0,6.0,0" \
        --noclosedcontourcmd "_AVG(T200,T500),-0.6,4,1.0" \
        --mergedist 6.0 \
        --searchbymin PSL \
        --outputcmd "PSL,min,0"
else
    echo "  ERROR: Grid resolution ${res} not supported for ETC detection!"
    echo "  Supported resolutions: 30 (ne30), 120 (ne120)"
    exit 1
fi

# Concatenate ETC detection results
echo "  Concatenating ETC detection files..."
cat "${result_dir}"etc_out.dat0*.dat > "${result_dir}etc_candidates_${file_name}.txt"

# Create ETC tracks
echo "  Creating ETC tracks..."
echo "  Duration: >= 60 hours, Max gap: 18 hours"
echo "  Min trajectory distance: 12°"
echo "  Range: 6°"

StitchNodes \
    --in_connect "${result_dir}connect_CSne${res}_v2.dat" \
    --in_fmt "lon,lat,slp" \
    --range 6.0 \
    --mintime "60h" \
    --maxgap 18 \
    --in "${result_dir}etc_candidates_${file_name}.txt" \
    --out "${result_dir}etc_tracks_${file_name}.dat" \
    --min_endpoint_dist 12.0

echo "  ETC tracking completed"

# Create ETC mask files
echo "  Creating ETC mask files..."
echo "  Binary mask: 1 within 6° of ETC center, 0 elsewhere"


NodeFileFilter \
    --in_connect "${result_dir}connect_CSne${res}_v2.dat" \
    --in_nodefile "${result_dir}etc_tracks_${file_name}.dat" \
    --in_fmt "lon,lat,slp" \
    --in_data_list "${result_dir}inputfile_${file_name}.txt" \
    --out_data_list "${result_dir}etc_mask_files_out.txt" \
    --bydist 6.0 \
    --maskvar "storm_mask"

echo "  ETC mask creation completed"
echo "  ETC detection, tracking, and masking completed"

#===============================================================================
# STEP 3: ATMOSPHERIC RIVER DETECTION (UNFILTERED)
#===============================================================================

echo "Step 3: Detecting atmospheric rivers..."
echo "  Method: Laplacian of integrated water vapor transport (IVT)"
echo "  Variables: TUQ (zonal IVT), TVQ (meridional IVT)"
echo "  Threshold: <= -20,000"
echo "  Min latitude: 15°"
echo "  Min area: 4e5 km²"

DetectBlobs \
    --in_connect "${result_dir}connect_CSne${res}_v2.dat" \
    --in_data_list "${result_dir}inputfile_${file_name}.txt" \
    --out_list "${result_dir}ar_nofilt_files_out.txt" \
    --thresholdcmd "_LAPLACIAN{8,10}(_VECMAG(TUQ,TVQ)),<=,-20000,0" \
    --minabslat 15 \
    --geofiltercmd "area,>=,4e5km2" \
    --tagvar "ar_mask"

echo "  AR detection completed"

#===============================================================================
# STEP 4: FILTER ARs TO REMOVE TC INFLUENCE
#===============================================================================

echo "Step 4: Filtering ARs to remove tropical cyclone influence..."
echo "  Exclusion distance: 8.0° from any TC center"
echo "  Method: Remove AR detections near TC tracks"

NodeFileFilter \
    --in_connect "${result_dir}connect_CSne${res}_v2.dat" \
    --in_nodefile "${result_dir}cyclones_stitch_${file_name}.dat" \
    --in_fmt "lon,lat,slp,wind" \
    --in_data_list "${result_dir}ar_nofilt_files_out.txt" \
    --out_data_list "${result_dir}ar_filt_files_out.txt" \
    --var "ar_mask" \
    --bydist 8.0 \
    --invert

echo "  AR filtering completed"

#===============================================================================
# STEP 5: APPLY AR MASK TO VAPOR TRANSPORT AND PRECT FIELD
#===============================================================================

echo "Step 5: Applying AR mask to vapor transport and precipitation field..."
echo "  Creating AR-tagged and non-AR vapor transport products"

VariableProcessor \
    --in_data_list "${result_dir}TVQ_PRECT_AR_mask_files_in.txt" \
    --out_data_list "${result_dir}TVQ_PRECT_AR_mask_files_out.txt" \
    --var "_PROD(ar_mask,TVQ);_PROD(_DIFF(1,ar_mask),TVQ);_PROD(ar_mask,PRECT)" \
    --varout "TVQ_AR,TVQ_NONAR,PRECT_AR" \
    --in_connect "${result_dir}connect_CSne${res}_v2.dat"

echo "  Vapor transport and precipitation masking completed"

#===============================================================================
# COMPLETION
#===============================================================================

echo "=============================================================================="
echo "Processing completed successfully!"
echo "=============================================================================="
echo "Output files:"
echo "  Unfiltered AR masks: ${result_dir}AR_mask_nofilt/"
echo "  TC-filtered AR masks: ${result_dir}AR_mask_filt/"
echo "  AR-masked TVQ/PRECT: ${result_dir}TVQ_PRECT_AR_mask/"
echo "  TC tracks: ${result_dir}cyclones_stitch_${file_name}.dat"
echo "  TC radial profiles: ${result_dir}cyclones_radprof_${file_name}.dat"
echo "  TC masks: ${result_dir}TC_masks/"
echo "  ETC tracks: ${result_dir}etc_tracks_${file_name}.dat"
echo "  ETC masks: ${result_dir}ETC_masks/"
echo "=============================================================================="

# Optional cleanup (uncomment if needed)
# echo "Cleaning up intermediate files..."
# rm -f "${result_dir}"out.dat0*.dat
# rm -f "${result_dir}cyclones_${file_name}.txt"
# rm -f "${result_dir}"*.txt

echo "Script completed at $(date)"
