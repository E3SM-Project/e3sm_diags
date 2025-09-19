import logging

import numpy as np
import xarray as xr
from scipy.signal import detrend

logging.basicConfig(level=logging.INFO)


def helper():
    """Prints all the functions that are included in this module."""
    f = [
        "decompose2SymAsym(arr)",
        "rmvAnnualCycle(data, spd, fCrit)",
        "smoothFrq121(data,nsmth_iter=1)",
        "smoothBackground_wavefreq(data)",
        "resolveWavesHayashi( varfft: xr.DataArray, nDayWin: int, spd: int ) -> xr.DataArray",
        "split_hann_taper(series_length, fraction)",
        "spacetime_power(data, segsize=96, noverlap=60, spd=1, latitude_bounds=None, dosymmetries=False, rmvLowFrq=False)",
        "genDispersionCurves(nWaveType=6, nPlanetaryWave=50, rlat=0, Ahe=[50, 25, 12])",
    ]
    [print(fe) for fe in f]


def getNearestInd1D(arr, val):
    """Given a specified value, find the index of an array that most closely matches that value.
    arr: numpy array or xarray DataArray
    return: Integer
    Example:  If arr=[0.2, 0.5, -0.9, 1.3] and val=0.4, function would return 1.
    Note:  Uses python's 'argmin', which returns first occurrence if multiple indices "tie"
           for closest to the specified value.
    """
    difference_array = np.absolute(arr - val)
    i_closest = int(difference_array.argmin())
    return i_closest


def decompose2SymAsym(arr):
    """Mimic NCL function to decompose into symmetric and asymmetric parts.

    arr: xarray DataArray

    return: DataArray with symmetric in SH, asymmetric in NH

    Note:
        This function produces indistinguishable results from NCL version.
    """
    lat_dim = arr.dims.index("lat")
    # flag to follow NCL convention and put symmetric component in SH
    # & asymmetric in NH
    # method: use flip to reverse latitude, put in DataArray for coords, use loc/isel
    # to assign to negative/positive latitudes (exact equator is left alone)
    data_sym = 0.5 * (arr.values + np.flip(arr.values, axis=lat_dim))
    data_asy = 0.5 * (arr.values - np.flip(arr.values, axis=lat_dim))
    data_sym = xr.DataArray(data_sym, dims=arr.dims, coords=arr.coords)
    data_asy = xr.DataArray(data_asy, dims=arr.dims, coords=arr.coords)
    out = arr.copy()  # might not be best to copy, but is safe
    out.loc[{"lat": arr["lat"][arr["lat"] < 0]}] = data_sym.isel(lat=data_sym.lat < 0)
    out.loc[{"lat": arr["lat"][arr["lat"] > 0]}] = data_asy.isel(lat=data_asy.lat > 0)
    return out


def rmvAnnualCycle(data, spd, fCrit):
    """remove frequencies less than fCrit from data.

    data: xarray DataArray
    spd: sampling frequency in samples-per-day
    fCrit: frequency threshold; remove frequencies < fCrit

    return: xarray DataArray, shape of data

    Note: fft/ifft preserves the mean because z = fft(x), z[0] is the mean.
          To keep the mean here, we need to keep the 0 frequency.

    Note: This function reproduces the results from the NCL version.

    Note: Two methods are available, one using fft/ifft and the other rfft/irfft.
          They both produce output that is indistinguishable from NCL's result.
    """
    dimz = data.sizes
    ntim = dimz["time"]
    time_ax = list(data.dims).index("time")
    # Method 1: Uses the complex FFT, returns the negative frequencies too, but they
    # should be redundant b/c they are conjugate of positive ones.
    cf = np.fft.fft(data.values, axis=time_ax)
    freq = np.fft.fftfreq(ntim, spd)
    cf[(freq != 0) & (np.abs(freq) < fCrit), ...] = 0.0  # keeps the mean
    z = np.fft.ifft(cf, n=ntim, axis=0)
    # Method 2: Uses the real FFT. In this case,
    # cf = np.fft.rfft(data.values, axis=time_ax)
    # freq = np.linspace(1, (ntim*spd)//2, (ntim*spd)//2) / ntim
    # fcrit_ndx = np.argwhere(freq < fCrit).max()
    # if fcrit_ndx > 1:
    #     cf[1:fcrit_ndx+1, ...] = 0.0
    # z = np.fft.irfft(cf, n=ntim, axis=0)
    da_z = xr.DataArray(z.real, dims=data.dims, coords=data.coords)
    return da_z


def smoothFrq121(data, nsmth_iter=1):
    """Following what was used in the NCL routine 'wkSpaceTime', smooth input array
    [nsmth_iter] times along the frequency dimension, for a sub-section of wavenumbers
    (as the smoothing is cosmetic for plotting, only smooth for a subset of wavenumbers
    to avoid unnecessary smoothing of non-plotted values).

    Do not use 0 frequency when smoothing.  Uses weights that sum to 1 to ensure
    smoothing is conservative.
    """
    assert isinstance(data, xr.DataArray)
    print(
        "\nFrom smoothFrq121:  Frequency smoothing "
        + str(nsmth_iter)
        + " times for subset of wavenumbers (pos frqs only)"
    )
    nfrq = len(data["frequency"])
    i_frq0 = int(
        np.where(data["frequency"] == 0)[0]
    )  # index of 'frequency' coord where it equals 0
    for smth_iter in range(nsmth_iter):
        data[..., i_frq0 + 1] = (
            0.75 * data[..., i_frq0 + 1] + 0.25 * data[..., i_frq0 + 2]
        )
        data[..., -1] = 0.25 * data[..., -2] + 0.75 * data[..., -1]
        for i in range(i_frq0 + 2, nfrq - 1):
            data[..., i] = (
                0.25 * data[..., i - 1] + 0.5 * data[..., i] + 0.25 * data[..., i + 1]
            )
    return data


def smoothBackground_wavefreq(data):
    """From Wheeler and Kiladis (1999, doi: ):
    "[Smooth] many times with a 1-2-1 filter in frequency and wavenumber. The number of
    passes of the 1-2-1 filter we have used is 10 in frequency throughout, and from 10
    to 40 in wavenumber, being 10 at low frequencies and 40 at higher frequencies
    increasing in two different steps."

    Here we do the following smoothing along the wavenumber dimension (note: we only
    smooth in the positive frequency domain as these will be the values plotted!):
    0 < frq < 0.1   :   5 1-2-1 smoothing passes
    0.1 <= frq < 0.2:  10 1-2-1 smoothing passes
    0.2 <= frq < 0.3:  20 1-2-1 smoothing passes
    0.3 <= frq      :  40 1-2-1 smoothing passes
    This follows what was used in the NCL routine 'wkSpaceTime'.
    """
    assert isinstance(data, xr.DataArray)
    nfrq = len(data["frequency"])
    i_frq0 = int(
        np.where(data["frequency"] == 0)[0]
    )  # index of 'frequency' coord where it equals 0
    #    i_frq03 = getNearestInd1D(data['frequency'], 0.3)   # index of 'frequency' coord closest to 0.3
    i_minwav4smth = int(
        np.where(data["wavenumber"] == -27)[0]
    )  # index of 'wavenumber' coord where it equals -27
    i_maxwav4smth = int(
        np.where(data["wavenumber"] == 27)[0]
    )  # index of 'wavenumber' coord where it equals 27

    # Looping over positive frequencies, smooth in wavenumber.  The number of smoothing
    #   passes in wavenumber is dependent on the frequency, with more smoothing at high
    #   frequencies and less smoothing at low frequencies
    print("\nSmoothing background spectrum in wavenumber (pos frq only)...")
    for i in range(i_frq0 + 1, nfrq):  # Loop over all positive frequencies
        if data["frequency"][i] < 0.1:
            nsmth_iter = 5
            print(
                "  Wavenumber smoothing "
                + str(nsmth_iter)
                + " times for freq: "
                + str(float(data["frequency"][i]))
            )
            for smth_iter in range(nsmth_iter):
                for j in range(i_minwav4smth, i_maxwav4smth + 1):
                    data[j, i] = (
                        0.25 * data[j - 1, i] + 0.5 * data[j, i] + 0.25 * data[j + 1, i]
                    )

        if data["frequency"][i] >= 0.1 and data["frequency"][i] < 0.2:
            nsmth_iter = 10
            print(
                "  Wavenumber smoothing "
                + str(nsmth_iter)
                + " times for freq: "
                + str(float(data["frequency"][i]))
            )
            for smth_iter in range(10):
                for j in range(i_minwav4smth, i_maxwav4smth + 1):
                    data[j, i] = (
                        0.25 * data[j - 1, i] + 0.5 * data[j, i] + 0.25 * data[j + 1, i]
                    )

        if data["frequency"][i] >= 0.2 and data["frequency"][i] < 0.3:
            nsmth_iter = 20
            print(
                "  Wavenumber smoothing "
                + str(nsmth_iter)
                + " times for freq: "
                + str(float(data["frequency"][i]))
            )
            for smth_iter in range(20):
                for j in range(i_minwav4smth, i_maxwav4smth + 1):
                    data[j, i] = (
                        0.25 * data[j - 1, i] + 0.5 * data[j, i] + 0.25 * data[j + 1, i]
                    )

        if data["frequency"][i] >= 0.3:
            nsmth_iter = 40
            print(
                "  Wavenumber smoothing "
                + str(nsmth_iter)
                + " times for freq: "
                + str(float(data["frequency"][i]))
            )
            for smth_iter in range(40):
                for j in range(i_minwav4smth, i_maxwav4smth + 1):
                    data[j, i] = (
                        0.25 * data[j - 1, i] + 0.5 * data[j, i] + 0.25 * data[j + 1, i]
                    )

    # For all wavenumbers, smooth in frequency:  10 passes
    #   Do not use 0 frequency when smoothing
    #   Use weights that sum to 1 to ensure smoothing is conservative
    nsmth_iter = 10
    print(
        "Frequency smoothing "
        + str(nsmth_iter)
        + " times for all wavenumbers (pos frqs only)"
    )
    for smth_iter in range(10):
        data[:, i_frq0 + 1] = 0.75 * data[:, i_frq0 + 1] + 0.25 * data[:, i_frq0 + 2]
        data[:, -1] = 0.25 * data[:, -2] + 0.75 * data[:, -1]
        for i in range(i_frq0 + 2, nfrq - 1):
            data[:, i] = (
                0.25 * data[:, i - 1] + 0.5 * data[:, i] + 0.25 * data[:, i + 1]
            )
    return data


def resolveWavesHayashi(varfft: xr.DataArray, nDayWin: int, spd: int) -> xr.DataArray:
    """This is a direct translation from the NCL routine to python/xarray.
    input:
        varfft : expected to have rightmost dimensions of wavenumber and frequency.
        varfft : expected to be an xarray DataArray with coordinate variables.
        nDayWin : integer that is the length of the segments in days.
        spd : the sampling frequency in `timesteps` per day.

    returns:
        a DataArray that is reordered to have correct westward & eastward propagation.

    """
    # -------------------------------------------------------------
    # Special reordering to resolve the Progressive and Retrogressive waves
    # Reference: Hayashi, Y.
    #    A Generalized Method of Resolving Disturbances into
    #    Progressive and Retrogressive Waves by Space and
    #    Fourier and TimeCross Spectral Analysis
    #    J. Meteor. Soc. Japan, 1971, 49: 125-128.
    # -------------------------------------------------------------

    # in NCL varfft is dimensioned (2,mlon,nSampWin), but the first dim doesn't matter b/c python supports complex numbers.
    #
    # Create array PEE(NL+1,NT+1) which contains the (real) power spectrum.
    # all the following assume indexing starting with 0
    # In this array (PEE), the negative wavenumbers will be from pn=0 to NL/2-1 (left).
    # The positive wavenumbers will be for pn=NL/2+1 to NL (right).
    # Negative frequencies will be from pt=0 to NT/2-1 (left).
    # Positive frequencies will be from pt=NT/2+1 to NT  (right).
    # Information about zonal mean will be for pn=NL/2 (middle).
    # Information about time mean will be for pt=NT/2 (middle).
    # Information about the Nyquist Frequency is at pt=0 and pt=NT
    #

    # In PEE, define the
    # WESTWARD waves to be either
    #          positive frequency and negative wavenumber
    #          OR
    #          negative freq and positive wavenumber.
    # EASTWARD waves are either positive freq and positive wavenumber
    #          OR negative freq and negative wavenumber.

    # Note that frequencies are returned from fftpack are ordered like so
    #    input_time_pos [ 0    1   2    3     4      5    6   7  ]
    #    ouput_fft_coef [mean 1/7 2/7  3/7 nyquist -3/7 -2/7 -1/7]
    #                    mean,pos freq to nyq,neg freq hi to lo
    #
    # Rearrange the coef array to give you power array of freq and wave number east/west
    # Note east/west wave number *NOT* eq to fft wavenumber see Hayashi '71
    # Hence, NCL's 'cfftf_frq_reorder' can *not* be used.
    # BPM: This goes for np.fft.fftshift
    #
    # For ffts that return the coefficients as described above, here is the algorithm
    # coeff array varfft(2,n,t)   dimensioned (2,0:numlon-1,0:numtim-1)
    # new space/time pee(2,pn,pt) dimensioned (2,0:numlon  ,0:numtim  )
    #
    # NOTE: one larger in both freq/space dims
    # the initial index of 2 is for the real (indx 0) and imag (indx 1) parts of the array
    #
    #
    #    if  |  0 <= pn <= numlon/2-1    then    | numlon/2 <= n <= 1
    #        |  0 <= pt < numtim/2-1             | numtim/2 <= t <= numtim-1
    #
    #    if  |  0         <= pn <= numlon/2-1    then    | numlon/2 <= n <= 1
    #        |  numtime/2 <= pt <= numtim                | 0        <= t <= numtim/2
    #
    #    if  |  numlon/2  <= pn <= numlon    then    | 0  <= n <= numlon/2
    #        |  0         <= pt <= numtim/2          | numtim/2 <= t <= 0
    #
    #    if  |  numlon/2   <= pn <= numlon    then    | 0        <= n <= numlon/2
    #        |  numtim/2+1 <= pt <= numtim            | numtim-1 <= t <= numtim/2

    # local variables : dimvf, numlon, N, varspacetime, pee, wave, freq

    # bpm: if varfft is a numpy array, then we need to know which dim is longitude
    #      if it is an xr.DataArray, then we can just use that directly. This is
    #      reason enough to insist on a DataArray.
    #      varfft should have a last dimension of "segments" of size N; should make a convention for the name of that dimension and insist on it here.
    logging.debug(f"[Hayashi] nDayWin: {nDayWin}, spd: {spd}")
    dimnames = varfft.dims
    dimvf = varfft.shape
    mlon = len(varfft["wavenumber"])
    N = dimvf[-1]
    logging.info(f"[Hayashi] input dims is {dimnames}, {dimvf}")
    logging.info(f"[Hayashi] input coords is {varfft.coords}")
    if len(dimnames) != len(varfft.coords):
        logging.error("The size of varfft.coords is incorrect.")
        raise ValueError("STOP")

    nshape = list(dimvf)
    nshape[-2] += 1
    nshape[-1] += 1
    logging.debug(f"[Hayashi] The nshape ends up being {nshape}")
    # this is a reordering, use Ellipsis to allow arbitrary number of dimensions,
    # but we insist that the wavenumber and frequency dims are rightmost.
    # we will fill the new array in increasing order (arbitrary choice)
    logging.debug("allocate the re-ordered array")
    varspacetime = np.full(nshape, np.nan, dtype=type(varfft))
    # first two are the negative wavenumbers (westward), second two are the positive wavenumbers (eastward)
    logging.debug(
        f"[Hayashi] Assign values into array. Notable numbers: mlon//2={mlon // 2}, N//2={N // 2}"
    )
    varspacetime[..., 0 : mlon // 2, 0 : N // 2] = varfft[
        ..., mlon // 2 : 0 : -1, N // 2 :
    ]  # neg.k, pos.w
    varspacetime[..., 0 : mlon // 2, N // 2 :] = varfft[
        ..., mlon // 2 : 0 : -1, 0 : N // 2 + 1
    ]  # pos.k,
    varspacetime[..., mlon // 2 :, 0 : N // 2 + 1] = varfft[
        ..., 0 : mlon // 2 + 1, N // 2 :: -1
    ]  # assign eastward & neg.freq.
    varspacetime[..., mlon // 2 :, N // 2 + 1 :] = varfft[
        ..., 0 : mlon // 2 + 1, -1 : N // 2 - 1 : -1
    ]  # assign eastward & pos.freq.
    print(varspacetime.shape)
    #  Create the real power spectrum pee = sqrt(real^2+imag^2)^2
    logging.debug("calculate power")
    pee = (
        (np.abs(varspacetime)) ** 2
    )  # JJB: abs(a+bi) = sqrt(a**2 + b**2), which is what is done in NCL's resolveWavesHayashi
    logging.debug("put into DataArray")
    # add meta data for use upon return
    wave = np.arange(-mlon // 2, (mlon // 2) + 1, 1, dtype=int)
    freq = (
        np.linspace(-1 * nDayWin * spd / 2, nDayWin * spd / 2, (nDayWin * spd) + 1)
        / nDayWin
    )

    print(f"freq size is {freq.shape}.")
    odims = list(dimnames)
    odims[-2] = "wavenumber"
    odims[-1] = "frequency"
    ocoords = {}
    for c in varfft.coords:
        logging.debug(f"[hayashi] working on coordinate {c}")
        if (c != "wavenumber") and (c != "frequency"):
            ocoords[c] = varfft[c]
        elif c == "wavenumber":
            # FIXME: mypy error: Incompatible types in assignment (expression has type "ndarray[tuple[int], dtype[Any]]", target has type "DataArray")  [assignment]
            ocoords["wavenumber"] = wave  # type: ignore
        elif c == "frequency":
            # FIXME: mypy error: Incompatible types in assignment (expression has type "ndarray[tuple[Any, ...], dtype[float64]]", target has type "DataArray")  [assignment]
            ocoords["frequency"] = freq  # type: ignore
    pee = xr.DataArray(pee, dims=odims, coords=ocoords)
    return pee


def split_hann_taper(series_length, fraction):
    """Implements `split cosine bell` taper of length series_length where only fraction of points are tapered (combined on both ends).

    This returns a function that tapers to zero on the ends. To taper to the mean of a series X:
    XTAPER = (X - X.mean())*series_taper + X.mean()
    """
    npts = int(np.rint(fraction * series_length))  # total size of taper
    taper = np.hanning(npts)
    series_taper = np.ones(series_length)
    series_taper[0 : npts // 2 + 1] = taper[0 : npts // 2 + 1]
    series_taper[-npts // 2 + 1 :] = taper[npts // 2 + 1 :]
    return series_taper


def spacetime_power(
    data,
    segsize=96,
    noverlap=60,
    spd=1,
    latitude_bounds=None,
    dosymmetries=False,
    rmvLowFrq=False,
):
    """Perform space-time spectral decomposition and return power spectrum following Wheeler-Kiladis approach.

    data: an xarray DataArray to be analyzed; needs to have (time, lat, lon) dimensions.
    segsize: integer denoting the size of time samples that will be decomposed (typically about 96)
    noverlap: integer denoting the number of days of overlap from one segment to the next (typically about segsize-60 => 2-month overlap)
    spd: sampling rate, in "samples per day" (e.g. daily=1, 6-houry=4)

    latitude_bounds: a tuple of (southern_extent, northern_extent) to reduce data size.

    dosymmetries: if True, follow NCL convention of putting symmetric component in SH, antisymmetric in NH
                  If True, the function returns a DataArray with a `component` dimension.

    rmvLowFrq: if True, remove frequencies < 1/segsize from data.

    Method
    ------
        1. Subsample in latitude if latitude_bounds is specified.
        2. Detrend the data (but keeps the mean value, as in NCL)
        3. High-pass filter if rmvLowFrq is True
        4. Construct symmetric/antisymmetric array if dosymmetries is True.
        5. Construct overlapping window view of data.
        6. Detrend the segments (strange enough, removing mean).
        7. Apply taper in time dimension of windows (aka segments).
        8. Fourier transform
        9. Apply Hayashi reordering to get propagation direction & convert to power.
       10. return DataArray with power.

    Notes
    -----
        Upon returning power, this should be comparable to "raw" spectra.
        Next step would be be to smooth with `smooth_wavefreq`,
        and divide raw spectra by smooth background to obtain "significant" spectral power.

    """

    segsize_steps = spd * segsize  # Segment length, in time step units
    noverlap_steps = spd * noverlap  # Segment overlap, in time step units
    stride_steps = (
        segsize_steps - noverlap_steps
    )  # Stride for overlapping windows, in time step units

    if latitude_bounds is not None:
        assert isinstance(latitude_bounds, tuple)
        data = data.sel(
            lat=slice(*latitude_bounds)
        )  # CAUTION: is this a mutable argument?
        logging.info(f"Data reduced by latitude bounds. Size is {data.sizes}")
        slat = latitude_bounds[0]
        nlat = latitude_bounds[1]
    else:
        slat = data["lat"].min().item()
        nlat = data["lat"].max().item()

    # Remove the *long-term* temporal linear trend
    #   Using scipy.signal.detrend will remove the mean as well, but we will add the time
    #   mean back into the detrended data to be consistent with the approach used in the
    #   NCL version (https://www.ncl.ucar.edu/Document/Functions/Diagnostics/wkSpaceTime.shtml):
    xmean = data.mean(dim="time").load()
    xdetr = detrend(data.values, axis=0, type="linear")
    xdetr = xr.DataArray(xdetr, dims=data.dims, coords=data.coords)
    xdetr += xmean  # put the mean back in
    # --> Tested and confirmed that this approach gives same answer as NCL

    # filter low-frequencies
    if rmvLowFrq:
        data = rmvAnnualCycle(xdetr, spd, 1 / segsize_steps)
    # --> Tested and confirmed that this function gives same answer as NCL

    # NOTE: at this point "data" has been modified to have its long-term linear trend and
    #       low frequencies (lower than 1./segsize_steps) removed

    dimsizes = data.sizes  # dict
    lon_size = dimsizes["lon"]
    lat_size = dimsizes["lat"]
    lat_dim = data.dims.index("lat")
    if dosymmetries:
        data = decompose2SymAsym(data)
    # testing: pass -- Gets the same result as NCL.

    # Windowing with the xarray "rolling" operation, and then limit overlap with `construct` to produce a new dataArray.
    # JJB: Windowing segment bounds have been validated with test prints.
    x_roll = data.rolling(
        time=segsize_steps, min_periods=segsize_steps
    )  # WK99 use 96-day window
    x_win = x_roll.construct("segments", stride=stride_steps).dropna(
        "time"
    )  # WK99 say "2-month" overlap; dims: (nSegments,nlat,nlon,segment_length_steps)
    logging.debug(f"[spacetime_power] x_win shape is {x_win.shape}")

    # For each segment, remove the temporal mean and linear trend:
    if np.logical_not(np.any(np.isnan(x_win))):
        logging.info("No missing, so use simplest segment detrend.")
        x_win_detr = detrend(
            x_win.values, axis=-1, type="linear"
        )  # <-- missing data makes this not work; axis=-1 for segment window dimension
        x_win = xr.DataArray(x_win_detr, dims=x_win.dims, coords=x_win.coords)
    else:
        logging.warning(
            "EXTREME WARNING -- This method to detrend with missing values present does not quite work, probably need to do interpolation instead."
        )
        logging.warning(
            "There are missing data in x_win, so have to try to detrend around them."
        )
        x_win_cp = x_win.values.copy()
        logging.info(
            f"[spacetime_power] x_win_cp windowed data has shape {x_win_cp.shape} \n \t It is a numpy array, copied from x_win which has dims: {x_win.sizes} \n \t ** about to detrend this in the rightmost dimension."
        )
        x_win_cp[np.logical_not(np.isnan(x_win_cp))] = detrend(
            x_win_cp[np.logical_not(np.isnan(x_win_cp))]
        )
        x_win = xr.DataArray(x_win_cp, dims=x_win.dims, coords=x_win.coords)

    # 3. Taper in time to reduce spectral "leakage" and make the segment periodic for FFT analysis
    # taper = np.hanning(segsize)  # WK seem to use some kind of stretched out hanning window; unclear if it matters
    taper = split_hann_taper(
        segsize_steps, 0.1
    )  # try to replicate NCL's approach of tapering 10% of segment window
    x_wintap = x_win * taper  # would do XTAPER = (X - X.mean())*series_taper + X.mean()
    # But since we have removed the mean, taper going to 0 is equivalent to taper going to the mean.

    # Do the transform using 2D FFT
    #   JJB:  As written below, this does not give the same answer as the 2-step approach,
    #         not sure why but didn't look into it more.  Will use 2-step approach to be
    #         safe (also mimics NCL version better)
    # z = np.fft.fft2(x_wintap, axes=(2,3)) / (lon_size * segsize)

    # Or do the transform with 2 steps
    #   Note:  x_wintap is shape (nSegments, nlat, nlon, segment_length_steps)
    z = (
        np.fft.fft(x_wintap, axis=2) / lon_size
    )  # note that np.fft.fft() produces same answers as NCL cfftf;  axis=2 to do fft along dim=lon
    z = (
        np.fft.fft(z, axis=3) / segsize_steps
    )  # axis=3 to do fft along dim=segment_length_steps

    z = xr.DataArray(
        z,
        dims=("time", "lat", "wavenumber", "frequency"),
        coords={
            "time": x_wintap["time"],
            "lat": x_wintap["lat"],
            "wavenumber": np.fft.fftfreq(lon_size, 1 / lon_size),
            "frequency": np.fft.fftfreq(segsize_steps, 1 / spd),
        },
    )
    #
    # The FFT is returned following ``standard order`` which has negative frequencies in second half of array.
    #
    # IMPORTANT:
    # If this were typical 2D FFT, we would do the following to get the frequencies and reorder:
    #         z_k = np.fft.fftfreq(x_wintap.shape[-2], 1/lon_size)
    #         z_v = np.fft.fftfreq(x_wintap.shape[-1], 1)  # Assumes 1/(1-day) timestep
    # reshape to get the frequencies centered
    #         z_centered = np.fft.fftshift(z, axes=(2,3))
    #         z_k_c = np.fft.fftshift(z_k)
    #         z_v_c = np.fft.fftshift(z_v)
    # and convert to DataArray as this:
    #         d1 = list(x_win.dims)
    #         d1[-2] = "wavenumber"
    #         d1[-1] = "frequency"
    #         c1 = {}
    #         for d in d1:
    #             if d in x_win.coords:
    #                 c1[d] = x_win[d]
    #             elif d == "wavenumber":
    #                 c1[d] = z_k_c
    #             elif d == "frequency":
    #                 c1[d] = z_v_c
    #         z_centered = xr.DataArray(z_centered, dims=d1, coords=c1)
    # BUT THAT IS INCORRECT TO GET THE PROPAGATION DIRECTION OF ZONAL WAVES
    # (in testing, it seems to end up opposite in wavenumber)
    # Apply reordering per Hayashi to get correct wave propagation convention
    #     this function is customized to expect z to be a DataArray
    z_pee = resolveWavesHayashi(z, segsize_steps // spd, spd)
    # z_pee is spectral power already.
    # z_pee is a DataArray w/ coordinate vars for wavenumber & frequency

    # average over all available segments and sum over latitude
    # OUTPUT DEPENDS ON SYMMETRIES
    if dosymmetries:
        # multipy by 2 b/c we only used one hemisphere
        # JJB: This summing of powers over latitude appears to be done a little bit earlier
        #      here when compared to NCL's version, but this should not matter.
        # JJB: I think the factor of 2 comes in here because the symmetric and antisymmetric
        #      data were packaged into either the Norther or Southern Hemis, so the 'sum'
        #      across latitude would need to be multiplied by 2.  For example, if the sym/asym data
        #      were -not- packaged into single hemispheres, no multiplicative factor would
        #      be needed.
        # JJB:  Also, 'squeeze' removed degenerate dimensions 'time' and 'lat'
        z_symmetric = (
            2.0
            * z_pee.isel(lat=z_pee.lat <= 0).mean(dim="time").sum(dim="lat").squeeze()
        )
        z_symmetric.name = "power"
        z_antisymmetric = (
            2.0
            * z_pee.isel(lat=z_pee.lat > 0).mean(dim="time").sum(dim="lat").squeeze()
        )
        z_antisymmetric.name = "power"
        z_final = xr.concat([z_symmetric, z_antisymmetric], "component")
        z_final = z_final.assign_coords({"component": ["symmetric", "antisymmetric"]})
        print("\nMetadata for z_final is:")
        print(z_final.dims)
        print(z_final.coords)
        print(z_final.attrs)
    else:
        lat = z_pee["lat"]
        lat_inds = np.argwhere(((lat <= nlat) & (lat >= slat)).values).squeeze()
        z_final = z_pee.isel(lat=lat_inds).mean(dim="time").sum(dim="lat").squeeze()
    return z_final


def genDispersionCurves(nWaveType=6, nPlanetaryWave=50, rlat=0, Ahe=[50, 25, 12]):  # noqa: C901
    """
    Function to derive the shallow water dispersion curves. Closely follows NCL version.

    input:
        nWaveType : integer, number of wave types to do
        nPlanetaryWave: integer
        rlat: latitude in radians (just one latitude, usually 0.0)
        Ahe: [50.,25.,12.] equivalent depths
              ==> defines parameter: nEquivDepth ; integer, number of equivalent depths to do == len(Ahe)

    returns: tuple of size 2
        Afreq: Frequency, shape is (nWaveType, nEquivDepth, nPlanetaryWave)
        Apzwn: Zonal savenumber, shape is (nWaveType, nEquivDepth, nPlanetaryWave)

    notes:
        The outputs contain both symmetric and antisymmetric waves. In the case of
        nWaveType == 6:
        0,1,2 are (ASYMMETRIC) "MRG", "IG", "EIG" (mixed rossby gravity, inertial gravity, equatorial inertial gravity)
        3,4,5 are (SYMMETRIC) "Kelvin", "ER", "IG" (Kelvin, equatorial rossby, inertial gravity)
    """
    nEquivDepth = len(Ahe)  # this was an input originally, but I don't know why.
    pi = np.pi
    radius = 6.37122e06  # [m]   average radius of earth
    g = 9.80665  # [m/s] gravity at 45 deg lat used by the WMO
    omega = 7.292e-05  # [1/s] earth's angular vel
    # U     = 0.0   # NOT USED, so Commented
    # Un    = 0.0   # since Un = U*T/L  # NOT USED, so Commented
    ll = 2.0 * pi * radius * np.cos(np.abs(rlat))
    Beta = 2.0 * omega * np.cos(np.abs(rlat)) / radius
    fillval = 1e20

    # NOTE: original code used a variable called del,
    #       I just replace that with `dell` because `del` is a python keyword.

    # Initialize the output arrays
    Afreq = np.empty((nWaveType, nEquivDepth, nPlanetaryWave))
    Apzwn = np.empty((nWaveType, nEquivDepth, nPlanetaryWave))

    for ww in range(1, nWaveType + 1):
        for ed, he in enumerate(Ahe):
            # this loops through the specified equivalent depths
            # ed provides index to fill in output array, while
            # he is the current equivalent depth
            # T = 1./np.sqrt(Beta)*(g*he)**(0.25) This is close to pre-factor of the dispersion relation, but is not used.
            c = np.sqrt(g * he)  # phase speed
            L = np.sqrt(
                c / Beta
            )  # was: (g*he)**(0.25)/np.sqrt(Beta), this is Rossby radius of deformation

            for wn in range(1, nPlanetaryWave + 1):
                s = -20.0 * (wn - 1) * 2.0 / (nPlanetaryWave - 1) + 20.0
                k = 2.0 * pi * s / ll
                kn = k * L

                # Anti-symmetric curves
                if ww == 1:  # MRG wave
                    if k < 0:
                        dell = np.sqrt(1.0 + (4.0 * Beta) / (k**2 * c))
                        deif = k * c * (0.5 - 0.5 * dell)

                    if k == 0:
                        deif = np.sqrt(c * Beta)

                    if k > 0:
                        deif = fillval

                if ww == 2:  # n=0 IG wave
                    if k < 0:
                        deif = fillval

                    if k == 0:
                        deif = np.sqrt(c * Beta)

                    if k > 0:
                        dell = np.sqrt(1.0 + (4.0 * Beta) / (k**2 * c))
                        deif = k * c * (0.5 + 0.5 * dell)

                if ww == 3:  # n=2 IG wave
                    n = 2.0
                    dell = Beta * c
                    deif = np.sqrt((2.0 * n + 1.0) * dell + (g * he) * k**2)
                    # do some corrections to the above calculated frequency.......
                    for i in range(1, 5 + 1):
                        deif = np.sqrt(
                            (2.0 * n + 1.0) * dell
                            + (g * he) * k**2
                            + g * he * Beta * k / deif
                        )

                # symmetric curves
                if ww == 4:  # n=1 ER wave
                    n = 1.0
                    if k < 0.0:
                        dell = (Beta / c) * (2.0 * n + 1.0)
                        deif = -Beta * k / (k**2 + dell)
                    else:
                        deif = fillval

                if ww == 5:  # Kelvin wave
                    deif = k * c

                if ww == 6:  # n=1 IG wave
                    n = 1.0
                    dell = Beta * c
                    deif = np.sqrt((2.0 * n + 1.0) * dell + (g * he) * k**2)
                    # do some corrections to the above calculated frequency
                    for i in range(1, 5 + 1):
                        deif = np.sqrt(
                            (2.0 * n + 1.0) * dell
                            + (g * he) * k**2
                            + g * he * Beta * k / deif
                        )

                eif = deif  # + k*U since  U=0.0
                P = 2.0 * pi / (eif * 24.0 * 60.0 * 60.0)  #  => PERIOD
                # dps  = deif/k  # Does not seem to be used.
                # R    = L #<-- this seemed unnecessary, I just changed R to L in Rdeg
                # Rdeg = (180.*L)/(pi*6.37e6) # And it doesn't get used.

                Apzwn[ww - 1, ed, wn - 1] = s
                if deif != fillval:
                    # P = 2.*pi/(eif*24.*60.*60.) # not sure why we would re-calculate now
                    Afreq[ww - 1, ed, wn - 1] = 1.0 / P
                else:
                    Afreq[ww - 1, ed, wn - 1] = fillval
    return Afreq, Apzwn
