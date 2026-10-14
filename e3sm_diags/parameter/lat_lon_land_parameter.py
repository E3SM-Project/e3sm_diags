from .core_parameter import CoreParameter


class LatLonLandParameter(CoreParameter):
    def __init__(self):
        super(LatLonLandParameter, self).__init__()

        # Work around the Cartopy >=0.25 masked ``contourf`` regression (#1026),
        # which floods NaN (ocean) regions or blanks valid land regions.
        self.transform_first: bool = True
