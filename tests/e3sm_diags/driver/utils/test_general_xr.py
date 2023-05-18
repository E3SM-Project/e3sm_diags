import numpy as np
import xarray as xr

from e3sm_diags.driver.utils.general_xr import has_z_axis


class TestHasZAxis:
    def test_returns_true_if_both_test_and_ref_arrays_have_z_axis(self):

        # Has Z axis
        z_axis1 = xr.DataArray(data=np.array([0]), attrs={"axis": "Z"})
        dv1 = xr.DataArray(data=[0], coords=[z_axis1])

        assert has_z_axis(dv1)

        # Has "positive" attribute equal to "up"
        z_axis2 = xr.DataArray(data=np.array([0]), attrs={"positive": "up"})
        dv2 = xr.DataArray(data=[0], coords=[z_axis2])

        assert has_z_axis(dv2)

        # Has "positive" attribute equal to "down"
        z_axis3 = xr.DataArray(data=np.array([0]), attrs={"positive": "down"})
        dv3 = xr.DataArray(data=[0], coords=[z_axis3])

        assert has_z_axis(dv3)

        # Has name equal to "lev"
        z_axis4 = xr.DataArray(name="lev", dims=["lev"], data=np.array([0]))
        dv4 = xr.DataArray(data=[0], coords={"lev": z_axis4})

        assert has_z_axis(dv4)

        # Has name equal to "plev"
        z_axis5 = xr.DataArray(name="plev", dims=["plev"], data=np.array([0]))
        dv5 = xr.DataArray(data=[0], coords=[z_axis5])

        assert has_z_axis(dv5)

        # Has name equal to "depth"
        z_axis6 = xr.DataArray(name="depth", dims=["depth"], data=np.array([0]))
        dv6 = xr.DataArray(data=[0], coords=[z_axis6])

        assert has_z_axis(dv6)

    def test_returns_false_if_either_test_or_ref_arrays_does_not_have_z_axis(self):
        dv1 = xr.DataArray(data=[0])

        assert not has_z_axis(dv1)
