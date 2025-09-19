#!/usr/bin/env python3
"""
Generate Land Configuration File with Actual NetCDF Data
========================================================

1. Keep original variables from lat_lon_land_model_vs_model.cfg (preserve contour levels)
2. Update case_id for original variables with CSV group information  
3. Add new variables from CSV with proper unit conversions and 5th/95th percentile contours from actual data
"""

import csv
import re
import numpy as np
import xarray as xr
from pathlib import Path

def read_csv_data(csv_path):
    """Read CSV file and return variable lookup dictionary."""
    csv_vars = {}
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            csv_vars[row['Variable']] = row
    return csv_vars

def parse_original_config(config_path):
    """Parse original config file and extract variable sections."""
    sections = []
    current_section = []
    current_var = None
    
    with open(config_path, 'r') as f:
        for line in f:
            line = line.rstrip()
            
            if line.startswith('[#]'):
                # Start of new section
                if current_section:
                    sections.append((current_var, '\n'.join(current_section)))
                current_section = [line]
                current_var = None
            elif line.startswith('variables = '):
                # Extract variable name
                match = re.search(r'variables = \["([^"]+)"\]', line)
                if match:
                    current_var = match.group(1)
                current_section.append(line)
            else:
                current_section.append(line)
        
        # Add final section
        if current_section and current_var:
            sections.append((current_var, '\n'.join(current_section)))
    
    return sections

def update_case_id_in_section(section_text, new_case_id):
    """Update case_id in a configuration section."""
    lines = section_text.split('\n')
    updated_lines = []
    
    for line in lines:
        if line.startswith('case_id = '):
            updated_lines.append(f'case_id = "{new_case_id}"')
        else:
            updated_lines.append(line)
    
    return '\n'.join(updated_lines)

def get_proper_conversion_factor(original_units, csv_scale_factor):
    """
    Get proper conversion factor based on comprehensive_unit_conversion_reference.txt
    Returns: (conversion_factor, target_units)
    """
    csv_factor = float(csv_scale_factor.replace('E', 'e'))
    original_units = original_units.strip()
    
    # Carbon flux variables: CSV factor 3.15360E-02 -> use 86400 for lat_lon maps  
    if abs(csv_factor - 3.15360E-02) < 1e-15:
        return 86400.0, original_units.replace('/s', '/day')  # gC/m²/s to gC/m²/day
    
    # Nitrogen/Phosphorus flux: CSV factor 3.15360E+01 -> use 86400 for lat_lon maps
    elif abs(csv_factor - 3.15360E+01) < 1e-10:
        return 86400.0, original_units.replace('/s', '/day')  # gN/m²/s to gN/m²/day, gP/m²/s to gP/m²/day
    
    # Carbon state variables: CSV factor 1.00000E-09 -> use 1/1000 for lat_lon maps
    elif abs(csv_factor - 1.00000E-09) < 1e-15:
        return 1.0/1000.0, original_units.replace('gC', 'kgC')  # gC/m² to kgC/m²
    
    # Nitrogen/Phosphorus state: CSV factor 1.00000E-06 -> use 1.0 for lat_lon maps
    elif abs(csv_factor - 1.00000E-06) < 1e-15:
        return 1.0, original_units  # Keep gN/m², gP/m²
    
    # Water flux: CSV factor 86400.0 -> already correct for mm/s to mm/day
    elif abs(csv_factor - 86400.0) < 1e-10:
        return 86400.0, 'mm/day'  # mm/s to mm/day
    
    # No conversion needed: CSV factor 1.0
    elif abs(csv_factor - 1.0) < 1e-10:
        return 1.0, original_units
    
    # Unknown - use as-is with warning
    else:
        print(f"WARNING: Unknown conversion pattern for {original_units} with CSV factor {csv_factor}")
        return 1.0, original_units

def calculate_percentile_contours_from_data(data_array, var_name, conversion_factor):
    """
    Calculate 5th/95th percentile-based contour levels from actual netCDF data.
    """
    try:
        print(f"  Calculating percentiles for {var_name} (conversion factor: {conversion_factor})")
        
        # Take temporal mean first to reduce to 2D (lat, lon)
        if 'time' in data_array.dims:
            data_2d = data_array.mean(dim='time', skipna=True)
        else:
            data_2d = data_array
        
        # Apply unit conversion
        if abs(conversion_factor - 1.0) > 1e-10:
            data_2d_converted = data_2d * conversion_factor
            print(f"    Applied conversion factor {conversion_factor}")
        else:
            data_2d_converted = data_2d
            print(f"    No conversion needed (factor = 1)")
        
        # Remove NaN values and flatten
        valid_data = data_2d_converted.values.flatten()
        valid_data = valid_data[~np.isnan(valid_data)]
        
        if len(valid_data) == 0:
            print(f"    WARNING: No valid data for {var_name}")
            return [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
        
        # Calculate percentiles
        p5 = np.percentile(valid_data, 5)
        p95 = np.percentile(valid_data, 95)
        p25 = np.percentile(valid_data, 25) 
        p75 = np.percentile(valid_data, 75)
        median = np.percentile(valid_data, 50)
        
        # Generate 16 contour levels using 5th-95th percentile range
        if abs(p95 - p5) < 1e-15:  # Nearly constant data
            delta = max(abs(median) * 0.01, 1e-10)
            levels = [median + (i-7.5) * delta/7.5 for i in range(16)]
        else:
            # Use 5th to 95th percentile range with some extension
            min_level = p5
            max_level = p95
            
            # Generate 16 levels spanning 5th to 95th percentiles
            levels = []
            for i in range(16):
                frac = i / 15.0
                level = min_level + frac * (max_level - min_level)
                levels.append(level)
        
        # Round to appropriate precision
        max_abs = max(abs(min(levels)), abs(max(levels)))
        if max_abs < 1e-10:
            precision = 15
        elif max_abs < 1e-6:
            precision = 10
        elif max_abs < 1e-3:
            precision = 8
        elif max_abs < 1:
            precision = 6
        elif max_abs < 1000:
            precision = 4
        else:
            precision = 2
        
        levels = [round(x, precision) for x in levels]
        
        print(f"    Data range: p5={p5:.4g}, median={median:.4g}, p95={p95:.4g}")
        print(f"    Contour levels: [{levels[0]:.4g}, ..., {levels[-1]:.4g}]")
        
        return levels
        
    except Exception as e:
        print(f"    ERROR calculating percentiles for {var_name}: {e}")
        # Return default levels
        return [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]

def generate_new_variable_section(var_name, csv_row, data_vars):
    """Generate configuration section for a new variable from CSV using actual data."""
    
    group = csv_row['Group']
    original_units = csv_row['Original Units'].strip()
    csv_scale_factor = csv_row['Scale Factor']
    long_name = csv_row['Long name'].strip()
    
    # Get proper conversion factor
    conversion_factor, target_units = get_proper_conversion_factor(original_units, csv_scale_factor)
    
    # Calculate contour levels from actual data
    if var_name in data_vars:
        contour_levels = calculate_percentile_contours_from_data(data_vars[var_name], var_name, conversion_factor)
    else:
        print(f"  WARNING: {var_name} not found in netCDF data, using default levels")
        # Use default levels based on conversion factor
        if conversion_factor == 86400.0:  # Flux variables
            contour_levels = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        elif conversion_factor == 1.0/1000.0:  # Carbon state
            contour_levels = [0, 1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25]
        else:  # Default
            contour_levels = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    
    contour_str = ",".join([f"{x:.6g}" for x in contour_levels])
    
    section_lines = [
        "[#]",
        'sets = ["lat_lon_land"]',
        f'case_id = "{group}"',
        f'variables = ["{var_name}"]',
        'seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]',
        'regions = ["global"]',
        'test_colormap = "WhiteBlueGreenYellowRed.rgb"',
        'reference_colormap = "WhiteBlueGreenYellowRed.rgb"',
        'diff_colormap = "BrBG"',
        f'contour_levels = [{contour_str}]'
    ]
    
    # Add metadata comments
    section_lines.extend([
        f'# group: {group}',
        f'# original_units: {original_units}',
        f'# target_units: {target_units}',
        f'# conversion_factor: {conversion_factor}',
        f'# csv_scale_factor: {csv_scale_factor}',
        f'# long_name: {long_name}',
        f'# contours: 5th/95th percentiles from PI control data'
    ])
    
    return '\n'.join(section_lines)

def main():
    """Main function to generate the land configuration file."""
    
    print("Generate Land Configuration File with Actual NetCDF Data")
    print("=" * 60)
    
    # File paths
    csv_path = "/tmp/zppy_land_fields.csv"
    original_config_path = "/gpfs/fs1/home/ac.zhang40/e3sm_diags/e3sm_diags/driver/default_diags/lat_lon_land_model_vs_model.cfg"
    netcdf_path = "/lcrc/group/e3sm2/ac.zhang40/E3SMv3/v3.LR.piControl_land_ilamb/post/lnd/native/clim/50yr/v3.LR.piControl_ANN_000101_005012_climo.nc"
    output_path = "/gpfs/fs1/home/ac.zhang40/e3sm_diags/auxiliary_tools/debug/999-add-all-land-var/lat_lon_land_model_vs_model.cfg"
    
    # Read CSV data
    csv_vars = read_csv_data(csv_path)
    print(f"Loaded {len(csv_vars)} variables from CSV")
    
    # Parse original configuration
    original_sections = parse_original_config(original_config_path)
    original_vars = {var: section for var, section in original_sections}
    print(f"Found {len(original_vars)} variables in original config")
    
    # Read netCDF data
    print(f"Reading netCDF file: {netcdf_path}")
    try:
        ds = xr.open_dataset(netcdf_path, decode_times=False)
        print(f"Successfully opened dataset with dimensions: {dict(ds.dims)}")
        
        # Get all data variables (exclude coordinate variables)
        coord_vars = {'lon', 'lat', 'time', 'area', 'topo', 'landfrac', 'landmask', 'pftmask'}
        data_vars = {}
        for var_name, var in ds.data_vars.items():
            if var_name not in coord_vars and 'time' in var.dims and len(var.dims) >= 2:
                data_vars[var_name] = var
        
        print(f"Found {len(data_vars)} data variables in netCDF")
        
    except Exception as e:
        print(f"ERROR reading netCDF file: {e}")
        print("Continuing without netCDF data - will use default contour levels")
        ds = None
        data_vars = {}
    
    # Categorize variables
    csv_var_names = set(csv_vars.keys())
    original_var_names = set(original_vars.keys())
    netcdf_var_names = set(data_vars.keys())
    
    new_vars = csv_var_names & netcdf_var_names - original_var_names  # New vars that exist in both CSV and netCDF
    
    print(f"Variables to keep from original: {len(original_var_names)}")
    print(f"New variables to add: {len(new_vars)}")
    
    # Generate output
    output_lines = []
    
    # Header
    output_lines.extend([
        "# ============================================================================",
        "# COMPLETE ELM LAND CONFIGURATION WITH ACTUAL DATA",
        "# ============================================================================",
        "# Generated from:",
        f"# - Original: {original_config_path}",
        f"# - CSV data: {csv_path}",
        f"# - NetCDF data: {netcdf_path}",
        "# - Unit conversions: comprehensive_unit_conversion_reference.txt",
        "#",
        "# Structure:",
        "# 1. Original variables (preserved settings + updated case_id with CSV group)",
        "# 2. New variables from CSV (proper unit conversions + 5th/95th percentile contours)",
        "# ============================================================================",
        ""
    ])
    
    # SECTION 1: Original variables with updated case_id
    output_lines.extend([
        "# ============================================================================",
        "# ORIGINAL VARIABLES (settings preserved, case_id updated with CSV group)",
        "# ============================================================================",
        ""
    ])
    
    for var_name in sorted(original_var_names):
        if var_name in csv_vars:
            # Update case_id with CSV group
            group = csv_vars[var_name]['Group']
            updated_section = update_case_id_in_section(original_vars[var_name], group)
            output_lines.append(updated_section)
            output_lines.append(f'# Updated case_id to: {group}')
        else:
            group = "Additional Variables"
            updated_section = update_case_id_in_section(original_vars[var_name], group)
            output_lines.append(updated_section)
            output_lines.append('# Original variable (no CSV group info)')
        
        output_lines.append("")
    
    # SECTION 2: New variables from CSV with actual data
    output_lines.extend([
        "# ============================================================================", 
        "# NEW VARIABLES FROM CSV (proper unit conversions + data-driven contours)",
        "# ============================================================================",
        ""
    ])
    
    # Group new variables by category
    vars_by_group = {}
    for var_name in sorted(new_vars):
        csv_row = csv_vars[var_name]
        group = csv_row['Group']
        if group not in vars_by_group:
            vars_by_group[group] = []
        vars_by_group[group].append((var_name, csv_row))
    
    for group_name in sorted(vars_by_group.keys()):
        group_vars = vars_by_group[group_name]
        output_lines.append(f"# {group_name.upper()} ({len(group_vars)} new variables)")
        output_lines.append("# " + "="*60)
        output_lines.append("")
        
        for var_name, csv_row in sorted(group_vars):
            print(f"Processing new variable: {var_name}")
            section = generate_new_variable_section(var_name, csv_row, data_vars)
            output_lines.append(section)
            output_lines.append("")
    
    # Write output file
    with open(output_path, 'w') as f:
        f.write('\n'.join(output_lines))
    
    # Close dataset if opened
    if ds is not None:
        ds.close()
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY:")
    print(f"- Original variables: {len(original_var_names)}")
    print(f"- New variables: {len(new_vars)}")
    print(f"- Total variables: {len(original_var_names) + len(new_vars)}")
    print(f"- New variable groups: {len(vars_by_group)}")
    print(f"\nOutput file: {output_path}")
    print("="*60)

if __name__ == "__main__":
    main()
