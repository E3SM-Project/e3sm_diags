import numpy as np
import numpy.ma as ma
import xarray as xr
import xcdat as xc

from e3sm_diags.driver.utils.climo_xr import CLIMO_CYCLE_MAP, ClimoFreq
from e3sm_diags.logger import _setup_child_logger

logger = _setup_child_logger(__name__)

SEASON_IDX = {
    "01": [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "02": [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "03": [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "04": [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    "05": [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    "06": [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    "07": [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    "08": [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    "09": [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    "10": [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    "11": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    "12": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    "DJF": [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    "MAM": [0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
    "JJA": [0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0],
    "SON": [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0],
    "ANN": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
}


def composite_diurnal_cycle(
    ds: xr.Dataset, var_key: str, season: ClimoFreq, fft: bool = True
) -> (
    tuple[ma.MaskedArray, np.ndarray] | tuple[xr.DataArray, xr.DataArray, xr.DataArray]
):
    """Compute the composite diurnal cycle for a variable for a given season.

    TODO: Add unit tests for this function.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset containing the variable.
    var_key : str
        The key of the variable.
    season : ClimoFreq
        The season for the climatology.
    fft : bool, optional
        Calculate using Fast Fourier transform, by default True.

    Returns
    -------
    tuple[ma.MaskedArray, np.ndarray] | tuple[xr.DataArray, xr.DataArray, xr.DataArray]
        Either a tuple containing the masked array for the diurnal cycle of the variable
        and the time coordinates as LST, or a tuple of three DataArrays for mean,
        amplitudes, and times-of-maximum of the first Fourier harmonic component
        (if ``fft=True``).
    """
    var = ds[var_key].copy()

    lat, lon = _get_lat_and_lon(var)
    time = _get_time(ds, var_key)
    time_freq, start_time = _get_time_freq_and_start_time(time)

    site = lat is None and lon is None
    if site:
        nlat = 1
        nlon = 1

        lat = [ds.lat.values]  # type: ignore
        lon = [ds.lon.values]  # type: ignore
    else:
        nlat = len(lat)  # type: ignore
        nlon = len(lon)  # type: ignore

    var_diurnal = _calc_var_diurnal(var, season, time, time_freq, site)

    # Convert GMT to LST
    nt = time_freq
    lst = np.zeros((nt, nlat, nlon))
    for it, itime in enumerate(np.arange(0, 24, 24 / nt)):
        for ilon in range(nlon):
            lst[it, :, ilon] = (itime + start_time + lon[ilon] / 360 * 24) % 24  # type: ignore

    # Compute mean, amplitude and max time of the first three Fourier components.
    if not fft:
        return var_diurnal, lst
    else:
        cycmean, maxvalue, tmax = _fft_all_grid(var_diurnal, lst)

        amplitude_data = np.zeros((nlat, nlon))
        amplitude_data[:, :] = maxvalue[0]
        amplitude = xr.DataArray(
            name="PRECT_diurnal_amplitude",
            data=amplitude_data,
            coords={lat.name: lat, lon.name: lon},  # type: ignore
            attrs={
                "longname": "Amplitude of diurnal cycle of PRECT",
                "units": var.units,
            },
        )

        maxtime_data = np.zeros((nlat, nlon))
        maxtime_data[:, :] = tmax[0]
        maxtime = xr.DataArray(
            name="PRECT_diurnal_phase",
            data=maxtime_data,
            coords={lat.name: lat, lon.name: lon},  # type: ignore
            attrs={
                "longname": "Phase of diurnal cycle of PRECT",
                "units": "hour",
            },
        )

        cmean_data = np.zeros((nlat, nlon))
        cmean_data[:, :] = cycmean
        cmean = xr.DataArray(
            name="PRECT_diurnal_cycmean",
            data=cmean_data,
            coords={lat.name: lat, lon.name: lon},  # type: ignore
            attrs={
                "longname": "Mean of diurnal cycle of PRECT",
                "units": "hour",
            },
        )

        return cmean, amplitude, maxtime


def _calc_var_diurnal(
    var: xr.DataArray, season: str, time: xr.DataArray, time_freq: int, site: bool
) -> ma.MaskedArray:
    cycle = CLIMO_CYCLE_MAP.get(season, [season])
    ncycle = len(cycle)

    time_coords_months = time.dt.month.values

    # var_diurnal has shape i.e. (ncycle, ntimesteps, [lat,lon]) for lat lon data
    var_diurnal = ma.zeros([ncycle] + [time_freq] + list(np.shape(var))[1:])

    for n in range(ncycle):
        # Get time index for each month/season.
        # Using a list comprehension to make looping faster, also
        # to have time_coords_months an array gets more speedup.
        time_idxs = np.array(
            [
                SEASON_IDX[cycle[n]][time_coords_months[i] - 1]
                for i in range(len(time_coords_months))
            ],
            dtype=np.int64,
        ).nonzero()

        var_reshape = np.reshape(
            var[time_idxs].values,
            (int(var[time_idxs].shape[0] / time_freq), time_freq)
            + var[time_idxs].shape[1:],
        )
        var_diurnal[n,] = ma.average(var_reshape, axis=0)

    if not site:
        var_diurnal = np.squeeze(var_diurnal)

    return var_diurnal


def _get_time(ds: xr.Dataset, var_key: str) -> xr.DataArray:
    ds_decoded = xr.decode_cf(ds, decode_times=True, use_cftime=True)

    try:
        time = xc.get_dim_coords(ds_decoded, axis="T")
    except (ValueError, KeyError) as err:
        raise KeyError(
            f"This variable ({var_key}) does not have a time axis. "
            "Climatology cannot be run on this variable."
        ) from err

    return time


def _get_time_freq_and_start_time(time: xr.DataArray) -> tuple[int, np.ndarray]:
    time_0 = time[0].dt.hour + time[0].dt.minute / 60 + time[0].dt.second / 3600
    time_1 = time[1].dt.hour + time[1].dt.minute / 60 + time[1].dt.second / 3600

    time_freq = int(24 / (time_1 - time_0))
    start_time = time_0

    logger.info(f"start_time {time.values[0]} {start_time.item()}")
    logger.info(f"var_time_freq={time_freq}")

    return time_freq, start_time.values


def _get_lat_and_lon(
    var: xr.DataArray,
) -> tuple[xr.DataArray | None, xr.DataArray | None]:
    lat = None
    lon = None

    try:
        lat = xc.get_dim_coords(var, axis="Y")
    except (ValueError, KeyError):
        pass

    try:
        lon = xc.get_dim_coords(var, axis="X")
    except (ValueError, KeyError):
        pass

    return lat, lon


def _fft_all_grid(
    var_diurnal: ma.MaskedArray, lst_time: np.ndarray
) -> tuple[ma.MaskedArray, np.ndarray, np.ndarray]:
    """Calculate Fast Fourier transform.

    This version of fastFT does all gridpoints at once. It uses a numerical
    Python function to compute a FAST Fourier transform, which should give the
    same result as a simple SLOW Fourier integration via the trapezoidal rule.

    Do NOT detrend the time series first, in order to retain the "sawtooth"
    frequency implied by the inumpy.t length of the time series (e.g. the
    24-hour period from a composite-diurnal cycle).

    On inumpy.t: x[k,i,j] = values at each gridpoint (i,j) for N times (k),
      - e.g. N = 8 for a 3-hr composite-diurnal cycle t[k,i,j] = timepoints
        at each gridpoint (i,j) for N times (k), e.g. Local Standard Times

    On output: c[i,j] = mean value at each gridpoint (i,j) in the time series
    ("zeroth" term in Fourier series)
      - maxvalue[n,i,j] = amplitude at each gridpoint (i,j) for each
                          Fourier harmonic (n)
      - tmax    [n,i,j] = time of maximum at each gridpoint (i,j) for each
                          Fourier harmonic (n)

    Source: Curt Covey, PCMDI/LLNL (December 2016)

    Parameters
    ----------
    var_diurnal : ma.MaskedArray
        The diurnal cycle for the variable.
    lst_time : np.ndarray
        A numpy array of LST time values.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        A tuple of numpy arrays for mean, amplitudes, and times-of-maximum of
        the first three Fourier harmonic components of the diurnal cycle of a
        variable.
    """
    # Creating output arrays
    if len(var_diurnal.shape) == 1:
        nx = 1
        ny = 1
    else:
        nx = var_diurnal.shape[1]
        ny = var_diurnal.shape[2]

    # time  of maximum for nth component (n=0 => diurnal, n=1 => semi...)
    tmax = np.zeros((3, nx, ny))
    # value of maximum for nth component (= 1/2 peak-to-peak amplitude)
    maxvalue = np.zeros((3, nx, ny))

    logger.info(
        "Calling numpy FFT function and converting from complex-valued FFT to real-valued amplitude and phase"
    )
    X = np.fft.ifft(var_diurnal, axis=0)
    logger.info("FFT output shape={}".format(X.shape))

    # Converting from complex-valued FFT to real-valued amplitude and phase
    a = X.real
    b = X.imag
    S = np.sqrt(a**2 + b**2)
    c = S[0]  # Zeroth harmonic = mean-value "constant term" in Fourier series.
    for n in range(3):
        # Adding first + last terms, second + second-to-last, ...
        maxvalue[n] = S[n + 1] + S[-n - 1]
        tmax[n] = np.arctan2(b[n + 1], a[n + 1])
        tmax[n] = tmax[n] * 12.0 / (np.pi * (n + 1))  # Radians to hours
        tmax[n] = tmax[n] + lst_time[0]  # GMT to LST
        tmax[n] = tmax[n] % (24 / (n + 1))

    return c, maxvalue, tmax
