# The type annotation for the metrics dictionary. The key is the
# type of metrics and the value is a sub-dictionary of metrics (key is metrics
# type and value is float). There is also a "unit" key representing the
# units for the variable.
UnitAttr = str
MetricsSubDict = dict[str, float | None | list[float]]
MetricsDict = dict[str, UnitAttr | MetricsSubDict]
