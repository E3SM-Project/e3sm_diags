"""
This aerosol budget set is requested by the E3SM Aerosol Working Group. The
script is integrated in e3sm_diags by Jill Zhang, with input from Kai Zhang,
Taufiq Hassan, Xue Zheng, Ziming Ke, Susannah Burrows, and Naser Mahfouz.
"""
from __future__ import annotations

import csv
import json
import os
from typing import TYPE_CHECKING

import numpy as np
import xarray as xr
import xcdat as xc

import e3sm_diags
from e3sm_diags.driver.utils.climo_xr import ClimoFreq
from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.io import _get_output_dir
from e3sm_diags.driver.utils.regrid import _hybrid_to_pressure
from e3sm_diags.driver.utils.type_annotations import MetricsDict
from e3sm_diags.logger import custom_logger

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter


logger = custom_logger(__name__)

# km units
REARTH = 6.37122e6
# kg/s to Tg/yr units.
UNITS_CONV = 86400.0 * 365.0 * 1e-9

SPECIES_NAMES = {
    "bc": "Black Carbon",
    "dst": "Dust",
    "mom": "Marine Organic Matter",
    "ncl": "Sea Salt",
    "pom": "Primary Organic Matter",
    "so4": "Sulfate",
    "soa": "Secondary Organic Aerosol",
}
MISSING_VALUE = 999.999


def run_diag(parameter: CoreParameter) -> CoreParameter:
    """Run the aerosol budget diagnostics.

    Parameters
    ----------
    parameter : CoreParameter
        The parameter for the diagnostic.

    Returns
    -------
    CoreParameter
        The parameter for the diagnostic with the result (completed or failed).

    Raises
    ------
    ValueError
        If the run type is not valid.
    """
    variables = parameter.variables[0].split(", ")
    run_type = parameter.run_type
    seasons = parameter.seasons

    for season in seasons:
        metrics_dict_test = {}
        metrics_dict_ref = {}
        test_ds = Dataset(parameter, data_type="test")
        parameter.test_name_yrs = test_ds.get_name_yrs_attr(season)
        parameter.ref_name_yrs = "Aerosol Global Benchmarks (Present Day)"

        for aerosol in variables:
            logger.info("Variable: {}".format(aerosol))
            metrics_dict_test[aerosol] = _create_metrics_dict(test_ds, aerosol, season)

            if run_type == "model_vs_model":
                ref_ds = Dataset(parameter, data_type="ref")
                parameter.ref_name_yrs = ref_ds.get_name_yrs_attr(season)

                for aerosol in variables:
                    metrics_dict_ref[aerosol] = _create_metrics_dict(
                        ref_ds, aerosol, season
                    )

            elif run_type == "model_vs_obs":
                parameter.ref_name = parameter.ref_name_yrs
                ref_data_path = os.path.join(
                    e3sm_diags.INSTALL_PATH,
                    "control_runs",
                    "aerosol_global_metrics_benchmarks.json",
                )

                with open(ref_data_path, "r") as myfile:
                    ref_file = myfile.read()

                metrics_ref = json.loads(ref_file)

                for aerosol in variables:
                    metrics_dict_ref[aerosol] = metrics_ref[aerosol]
            else:
                raise ValueError("Invalid run_type={}".format(run_type))

            parameter.output_file = f"{parameter.test_name}-{season}-budget-table"
            fnm = os.path.join(
                _get_output_dir(parameter),
                parameter.output_file + ".csv",
            )

            with open(fnm, "w") as table_csv:
                writer = csv.writer(
                    table_csv,
                    delimiter=",",
                    quotechar="'",
                    quoting=csv.QUOTE_MINIMAL,
                    lineterminator="\n",
                )
                writer.writerow([f"Test: {parameter.test_name_yrs}"])
                writer.writerow([f"Ref: {parameter.ref_name_yrs}"])
                writer.writerow([" ", "Test", "Ref"])

                for key, values in metrics_dict_test.items():
                    writer.writerow([SPECIES_NAMES[key]])

                    for value in values:
                        line = []
                        line.append(value)
                        line.append(values[value])  # type: ignore
                        line.append(metrics_dict_ref[key][value])  # type: ignore
                        writer.writerows([line])

                    writer.writerows([""])

    logger.info(f"Metrics saved in {fnm}")

    return parameter


def _create_metrics_dict(
    ds_test: Dataset, aerosol: str, season: ClimoFreq
) -> MetricsDict:
    wetdep = ds_test.get_climo_dataset(f"{aerosol}_SFWET", season)[f"{aerosol}_SFWET"]
    drydep = ds_test.get_climo_dataset(f"{aerosol}_DDF", season)[f"{aerosol}_DDF"]
    srfemis = ds_test.get_climo_dataset(f"SF{aerosol}", season)[f"SF{aerosol}"]
    area = ds_test.get_climo_dataset(f"{aerosol}_DDF", season)["area"]

    area_m2 = area * REARTH**2

    burden = _calc_column_integral(ds_test, aerosol, season)
    # Convert kg to Tg.
    burden_total = global_integral(burden, area_m2) * 1e-9

    srfemis_gi = global_integral(srfemis, area_m2) * UNITS_CONV
    sink_gi = global_integral((drydep - wetdep), area_m2) * UNITS_CONV
    drydep_gi = global_integral(drydep, area_m2) * UNITS_CONV
    wetdep_gi = global_integral(wetdep, area_m2) * UNITS_CONV

    if aerosol in ["bc", "pom", "so4"]:
        var_key = f"{aerosol}_CLXF"
        elvemis = ds_test.get_climo_dataset(var_key, season)[var_key]

        elvemis_gi = global_integral(elvemis, area_m2) * UNITS_CONV
    else:
        elvemis_gi = 0.0

    metrics_dict: MetricsDict = {
        "Surface Emission (Tg/yr)": f"{srfemis_gi:.2f}",
        "Elevated Emission (Tg/yr)": f"{elvemis_gi:.2f}",
        "Sink (Tg/yr)": f"{sink_gi:.2f}",
        "Dry Deposition (Tg/yr)": f"{drydep_gi:.2f}",
        "Wet Deposition (Tg/yr)": f"{wetdep_gi:.2f}",
        "Burden (Tg)": f"{burden_total:.2f}",
        "Lifetime (Days)": f"{burden_total/sink_gi*365:.2f}",
    }

    return metrics_dict


def _calc_column_integral(
    test_ds: Dataset, aerosol: str, season: ClimoFreq
) -> xr.DataArray:
    """Calculate column integrated mass

    Take aerosol and change it to the appropriate string.
      - ncl -> SEASALT, dst -> DUST, rest1 -> REST1
      - only needed if ABURDEN terms are available

    Parameters
    ----------
    test_ds : Dataset
        The e3sm_diags test dataset object.
    aerosol : str
        The aerosol key.
    season : ClimoFreq.
        The climo season.

    Returns
    -------
    xr.DataArray
        The column integrated mass.
    """
    if aerosol == "ncl":
        aerosol_name = "SEASALT"
    elif aerosol == "dst":
        aerosol_name = "DUST"
    else:
        aerosol_name = aerosol.upper()

    try:
        # if ABURDEN terms are available, use them
        var_key = f"ABURDEN{aerosol_name}"
        burden = test_ds.get_climo_dataset(var_key, season)[var_key]
    except IOError:
        # if not, use the Mass_ terms and integrate over the column
        var_key = f"Mass_{aerosol}"
        ds_mass = test_ds.get_climo_dataset(var_key, season)
        mass = ds_mass[var_key]

        if mass.attrs["units"] != "kg/kg":
            raise RuntimeError(
                f"ERROR in aerosol_budget_driver/calc_column_integral!"
                f"Mass_{aerosol} at season {season} has units {mass.units}."
                f"But kg/kg units were expected."
            )

        pressure_levs = _hybrid_to_pressure(
            ds_mass, var_key, p0=100000.0, a_key="hyai", b_key="hybi"
        )
        # Preserve the original units.
        pressure_levs.attrs["units"] = ds_mass[var_key].attrs["units"]

        # (72,lat,lon)
        plev_key = xc.get_dim_keys(pressure_levs, axis="Z")
        delta_p = pressure_levs.diff(plev_key)

        # Perform numpy index-based arithmetic instead of xarray label-based
        # arithmetic because the Z dims of `mass` and `delta_p` can have
        # different names ("lev" vs. "ilev") with slightly different values. The
        # CDAT version of this code uses numpy index-based arithmetic too.
        with xr.set_options(keep_attrs=True):
            # mass density * mass air   kg/m2
            mass_3d = mass.copy()
            mass_3d.values = mass.values * delta_p.values / 9.8

        # kg/m2
        lev_key = xc.get_dim_keys(mass_3d, axis="Z")
        burden = mass_3d.sum(dim=lev_key, keep_attrs=True)

    return burden


def global_integral(var: xr.DataArray, area_m2: xr.DataArray) -> float:
    """Compute global integral of 2 dimensional properties

    Parameters
    ----------
    var : xr.DataArray
        The variable.
    area_m2 : xr.DataArray
        The area in m2 units.

    Returns
    -------
    float
        The global integral of 2 dimensional properities
    """
    result = np.sum(np.sum(abs(var) * area_m2, axis=0), axis=0)

    return result
