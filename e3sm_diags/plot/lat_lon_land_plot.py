from e3sm_diags.plot.lat_lon_plot import plot as base_plot


def plot(reference, test, diff, metrics_dict, parameter):
    return base_plot(reference, test, diff, metrics_dict, parameter)
