# The type annotation for the metrics dictionary. The key is the
# type of metrics and the value is a sub-dictionary of metrics (key is metrics
# type and value is float). There is also a "unit" key representing the
# units for the variable.
from e3sm_diags.driver.utils.climo_xr import ClimoFreq

UnitAttr = str
MetricsSubDict = dict[str, float | None | list[float]]
MetricsDict = dict[str, UnitAttr | MetricsSubDict]

# Type for time slice specification: individual time index for snapshot analysis
# Examples: "0", "5", "42"
TimeSlice = str

# Union type for time selection - can be either climatology season or time slice
TimeSelection = ClimoFreq | TimeSlice
