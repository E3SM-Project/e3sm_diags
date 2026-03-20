#!/usr/bin/env python3
"""
Simple validation script comparing original vs COSP mp_partition algorithms.
Uses first time step of 1850 data for quick validation.
"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from e3sm_diags.driver.mp_partition_driver import compute_lcf, compute_lcf_cosp

# Data paths
DATA_DIR = "/lcrc/group/e3sm2/ac.zhang40/E3SMv3/v3.LR.historical_0101/post/atm/180x360_aave/ts/monthly/5yr"

def load_validation_data():
    """Load first time step of 1850 data for both COSP and original variables."""
    print("=" * 60)
    print("LOADING 1850 VALIDATION DATA")
    print("=" * 60)

    data = {}

    # Load COSP variables
    print("Loading COSP variables...")
    try:
        # COSP cloud ice
        cice_cosp_file = f"{DATA_DIR}/CLD_CAL_TMPICE_185001_185412.nc"
        cice_cosp_ds = xr.open_dataset(cice_cosp_file)
        data['cice_cosp'] = cice_cosp_ds['CLD_CAL_TMPICE'].isel(time=0)  # First time step
        data['cosp_temp'] = cice_cosp_ds['cosp_temp']

        # COSP cloud liquid
        cliq_cosp_file = f"{DATA_DIR}/CLD_CAL_TMPLIQ_185001_185412.nc"
        cliq_cosp_ds = xr.open_dataset(cliq_cosp_file)
        data['cliq_cosp'] = cliq_cosp_ds['CLD_CAL_TMPLIQ'].isel(time=0)  # First time step

        print(f"✅ COSP data loaded: {data['cice_cosp'].shape}")
        print(f"   Temperature bins: {len(data['cosp_temp'])}")
        data['cosp_loaded'] = True

    except Exception as e:
        print(f"❌ Failed to load COSP variables: {e}")
        data['cosp_loaded'] = False

    # Load original variables
    print("Loading original variables...")
    try:
        # Original cloud ice
        cice_orig_file = f"{DATA_DIR}/CLDICE_185001_185412.nc"
        cice_orig_ds = xr.open_dataset(cice_orig_file)
        data['cice_orig'] = cice_orig_ds['CLDICE'].isel(time=0)  # First time step

        # Original cloud liquid
        cliq_orig_file = f"{DATA_DIR}/CLDLIQ_185001_185412.nc"
        cliq_orig_ds = xr.open_dataset(cliq_orig_file)
        data['cliq_orig'] = cliq_orig_ds['CLDLIQ'].isel(time=0)  # First time step

        # Temperature
        temp_orig_file = f"{DATA_DIR}/T_185001_185412.nc"
        temp_orig_ds = xr.open_dataset(temp_orig_file)
        data['temp_orig'] = temp_orig_ds['T'].isel(time=0)  # First time step

        print(f"✅ Original data loaded: {data['cice_orig'].shape}")
        data['orig_loaded'] = True

    except Exception as e:
        print(f"❌ Failed to load original variables: {e}")
        data['orig_loaded'] = False

    # Load land fraction
    print("Loading land fraction...")
    try:
        landfrac_file = f"{DATA_DIR}/LANDFRAC_185001_185412.nc"
        landfrac_ds = xr.open_dataset(landfrac_file)
        data['landfrac'] = landfrac_ds['LANDFRAC'].isel(time=0)  # First time step
        print(f"✅ Land fraction loaded: {data['landfrac'].shape}")
    except Exception as e:
        print(f"❌ Failed to load LANDFRAC: {e}")
        data['landfrac'] = None

    return data

def compare_data_ranges(data):
    """Compare data ranges and characteristics."""
    print("\n" + "=" * 60)
    print("DATA CHARACTERISTICS")
    print("=" * 60)

    if data['cosp_loaded']:
        print("COSP VARIABLES:")
        print(f"  CLD_CAL_TMPICE range: {data['cice_cosp'].min().values:.6f} - {data['cice_cosp'].max().values:.6f}")
        print(f"  CLD_CAL_TMPLIQ range: {data['cliq_cosp'].min().values:.6f} - {data['cliq_cosp'].max().values:.6f}")
        print(f"  cosp_temp range (C):  {data['cosp_temp'].min().values:.1f} - {data['cosp_temp'].max().values:.1f}")
        print(f"  cosp_temp range (K):  {(data['cosp_temp'] + 273.15).min().values:.1f} - {(data['cosp_temp'] + 273.15).max().values:.1f}")

    if data['orig_loaded']:
        print("ORIGINAL VARIABLES:")
        print(f"  CLDICE range:         {data['cice_orig'].min().values:.6f} - {data['cice_orig'].max().values:.6f}")
        print(f"  CLDLIQ range:         {data['cliq_orig'].min().values:.6f} - {data['cliq_orig'].max().values:.6f}")
        print(f"  T range (K):          {data['temp_orig'].min().values:.1f} - {data['temp_orig'].max().values:.1f}")

    if data['cosp_loaded'] and data['orig_loaded']:
        print("COMPARISON:")
        print(f"  Data shapes:")
        print(f"    COSP:     {data['cice_cosp'].shape} (cosp_temp, lat, lon)")
        print(f"    Original: {data['cice_orig'].shape} (lev, lat, lon)")

        # Valid point comparison
        cosp_total = (data['cice_cosp'] + data['cliq_cosp'])
        orig_total = (data['cice_orig'] + data['cliq_orig'])

        print(f"  Valid points (>1e-9):")
        print(f"    COSP:     {(cosp_total > 1e-9).sum().values:,}")
        print(f"    Original: {(orig_total > 1e-9).sum().values:,}")

        print(f"  Valid points (>0.001):")
        print(f"    COSP:     {(cosp_total > 0.001).sum().values:,}")
        print(f"    Original: {(orig_total > 0.001).sum().values:,}")

def run_validation(data):
    """Run both algorithms and compare."""
    print("\n" + "=" * 60)
    print("RUNNING ALGORITHMS")
    print("=" * 60)

    results = {}

    # Add time dimension back for algorithms (they expect time dimension)
    if data['cosp_loaded']:
        cice_cosp = data['cice_cosp'].expand_dims('time')
        cliq_cosp = data['cliq_cosp'].expand_dims('time')

        print("Running COSP algorithm...")
        try:
            temp_centers_cosp, lcf_cosp = compute_lcf_cosp(
                cice_cosp, cliq_cosp, data['cosp_temp'], data['landfrac']
            )
            results['cosp'] = {
                'temp_centers': temp_centers_cosp,
                'lcf': lcf_cosp,
                'success': True
            }
            print(f"✅ COSP completed - LCF range: {lcf_cosp.min():.4f} - {lcf_cosp.max():.4f}")
            print(f"   Non-zero bins: {np.sum(lcf_cosp > 0)}/20")

        except Exception as e:
            print(f"❌ COSP failed: {e}")
            results['cosp'] = {'success': False, 'error': str(e)}

    if data['orig_loaded']:
        cice_orig = data['cice_orig'].expand_dims('time')
        cliq_orig = data['cliq_orig'].expand_dims('time')
        temp_orig = data['temp_orig'].expand_dims('time')

        print("Running original algorithm...")
        try:
            temp_centers_orig, lcf_orig = compute_lcf(
                cice_orig, cliq_orig, temp_orig, data['landfrac']
            )
            results['orig'] = {
                'temp_centers': temp_centers_orig,
                'lcf': lcf_orig,
                'success': True
            }
            print(f"✅ Original completed - LCF range: {lcf_orig.min():.4f} - {lcf_orig.max():.4f}")
            print(f"   Non-zero bins: {np.sum(lcf_orig > 0)}/20")

        except Exception as e:
            print(f"❌ Original failed: {e}")
            results['orig'] = {'success': False, 'error': str(e)}

    return results

def plot_results(results):
    """Plot comparison of results."""
    print("\n" + "=" * 60)
    print("GENERATING PLOTS")
    print("=" * 60)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Plot LCF curves
    if results.get('orig', {}).get('success'):
        temp_c = results['orig']['temp_centers'] - 273.15
        ax1.plot(temp_c, results['orig']['lcf'], 'o-',
                label='Original (CLDICE/CLDLIQ)', linewidth=2, markersize=6, color='blue')

    if results.get('cosp', {}).get('success'):
        temp_c = results['cosp']['temp_centers'] - 273.15
        ax1.plot(temp_c, results['cosp']['lcf'], 's-',
                label='COSP (CLD_CAL)', linewidth=2, markersize=6, color='red')

    ax1.set_xlabel('Temperature (°C)')
    ax1.set_ylabel('Liquid Condensate Fraction')
    ax1.set_title('Algorithm Comparison - 1850-01 (First Time Step)')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_xlim(-60, 10)
    ax1.set_ylim(0, 1)

    # Plot difference if both succeeded
    if (results.get('orig', {}).get('success') and
        results.get('cosp', {}).get('success')):

        # Interpolate original to COSP temperature grid
        lcf_orig_interp = np.interp(
            results['cosp']['temp_centers'],
            results['orig']['temp_centers'],
            results['orig']['lcf']
        )
        diff = results['cosp']['lcf'] - lcf_orig_interp
        temp_c = results['cosp']['temp_centers'] - 273.15

        ax2.plot(temp_c, diff, 'go-', linewidth=2, markersize=6)
        ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5)
        ax2.set_xlabel('Temperature (°C)')
        ax2.set_ylabel('LCF Difference (COSP - Original)')
        ax2.set_title('Difference Plot')
        ax2.grid(True, alpha=0.3)

        # Statistics
        rmse = np.sqrt(np.mean(diff**2))
        bias = np.mean(diff)
        ax2.text(0.05, 0.95, f'RMSE: {rmse:.4f}\nBias: {bias:.4f}',
                transform=ax2.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    else:
        ax2.text(0.5, 0.5, 'Cannot compute difference:\nOne or both algorithms failed',
                ha='center', va='center', transform=ax2.transAxes)
        ax2.set_title('Difference Not Available')

    plt.tight_layout()

    # Save plot to specified directory
    import os
    output_dir = "/lcrc/group/e3sm/public_html/diagnostic_output/ac.zhang40/tests/"
    os.makedirs(output_dir, exist_ok=True)

    output_file = os.path.join(output_dir, "mp_partition_validation_1850.png")
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {output_file}")
    plt.show()

    # Print numerical comparison
    if (results.get('orig', {}).get('success') and
        results.get('cosp', {}).get('success')):

        print("\n📊 NUMERICAL COMPARISON:")
        print("Temperature (°C) | Original LCF | COSP LCF | Difference")
        print("-" * 55)

        # Interpolate for comparison
        lcf_orig_interp = np.interp(
            results['cosp']['temp_centers'],
            results['orig']['temp_centers'],
            results['orig']['lcf']
        )

        for i, temp_k in enumerate(results['cosp']['temp_centers']):
            temp_c = temp_k - 273.15
            orig_lcf = lcf_orig_interp[i]
            cosp_lcf = results['cosp']['lcf'][i]
            diff = cosp_lcf - orig_lcf
            print(f"{temp_c:8.1f}     | {orig_lcf:8.4f}     | {cosp_lcf:8.4f} | {diff:8.4f}")

def main():
    """Main validation function."""
    print("🔍 MP PARTITION VALIDATION - 1850 Data")
    print("Using first time step for quick comparison")

    # Load data
    data = load_validation_data()

    # Compare data characteristics
    compare_data_ranges(data)

    # Run validation
    results = run_validation(data)

    # Plot results
    plot_results(results)

    print("\n🎉 Validation complete!")

if __name__ == "__main__":
    main()