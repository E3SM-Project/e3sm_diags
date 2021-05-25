from __future__ import print_function

from typing import Any, Dict

import cdms2
import cdutil
import MV2

from acme_diags.driver import utils
from acme_diags.plot import plot


def create_annual_cycle(dataset, variable):
    month_list = [f"{x:02}" for x in list(range(1, 13))]

    for imon, month in enumerate(month_list):
        # fin = cdms2.open(file_path)
        # var = fin(variable)
        # fin.close()
        var = dataset.get_climo_variable(variable, month)
        if month == "01":
            var_ac = MV2.zeros([12] + list(var.shape)[:])
            var_ac.id = var.id
            var_ac.long_name = var.long_name
            var_ac.units = var.units
            time = cdms2.createAxis(range(1, 13))
            time.id = "time"
            var_ac.setAxis(0, time)
            time.designateTime()
            var_ac.setAxis(1, var.getAxis(0))
            var_ac.setAxis(2, var.getAxis(1))

            var_ac[
                0,
            ] = var
        else:
            var_ac[
                imon,
            ] = var

    return var_ac


def run_diag(parameter):
    variables = parameter.variables
    ref_name = getattr(parameter, "ref_name", "")

    test_data = utils.dataset.Dataset(parameter, test=True)
    ref_data = utils.dataset.Dataset(parameter, ref=True)

    parameter.test_name_yrs = utils.general.get_name_and_yrs(parameter, test_data, "01")
    parameter.ref_name_yrs = utils.general.get_name_and_yrs(parameter, ref_data, "01")

    for var in variables:
        test_ac = create_annual_cycle(test_data, var)
        ref_ac = create_annual_cycle(ref_data, var)

        test_ac_reg, ref_ac_reg = utils.general.regrid_to_lower_res(
            test_ac,
            ref_ac,
            parameter.regrid_tool,
            parameter.regrid_method,
        )

        test_ac_zonal_mean = cdutil.averager(test_ac, axis="x", weights="generate")
        ref_ac_zonal_mean = cdutil.averager(ref_ac, axis="x", weights="generate")
        test_ac_reg_zonal_mean = cdutil.averager(
            test_ac_reg, axis="x", weights="generate"
        )
        ref_ac_reg_zonal_mean = cdutil.averager(
            ref_ac_reg, axis="x", weights="generate"
        )

        diff_ac = test_ac_reg_zonal_mean - ref_ac_reg_zonal_mean

        parameter.var_id = var
        parameter.output_file = "-".join([ref_name, var, "Annual-Cycle"])
        parameter.main_title = str(" ".join([var, "Zonel Mean Annual Cycle"]))

        parameter.viewer_descr[var] = (
            test_ac.long_name
            if hasattr(test_ac, "long_name")
            else "No long_name attr in test data."
        )

        metrics_dict: Dict[str, Any] = {}

        plot(
            parameter.current_set,
            ref_ac_zonal_mean,
            test_ac_zonal_mean,
            diff_ac,
            metrics_dict,
            parameter,
        )
        utils.general.save_ncfiles(
            parameter.current_set,
            ref_ac_zonal_mean,
            test_ac_zonal_mean,
            diff_ac,
            parameter,
        )

    return parameter
