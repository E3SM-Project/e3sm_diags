# Targeted Plot Baselines

These PNG baselines are expected outputs for the targeted image-regression tests.
They use deterministic synthetic inputs from
[plot_image_regression_case.py](/Users/vo13/Repositories/E3SM-Project/e3sm_diags/tests/integration/plot_image_regression_case.py:1)
and validate plot rendering without downloaded integration data.

Committed baselines and `baseline_metadata.json` should be refreshed with the
same `conda-env/ci.yml` and Python 3.13 environment used by the main GitHub
Actions Layer 2 visual-regression job. The repository also provides a manual
`Update Image Baselines` workflow that regenerates these files directly on
`main` using that same authority.

## What Each Plot Represents

- `lat_lon_plot/`: Regional latitude-longitude map over the western Pacific. Checks map layout, coastlines, contour rendering, and regional axis labeling for a stable non-dateline case.
- `polar_plot/`: Northern polar stereographic map. Checks polar projection rendering, circular crop, coastlines, and structured polar anomalies.
- `zonal_mean_2d_plot/`: Latitude-pressure zonal mean cross section. Checks vertical section contours, pressure-axis layout, and smooth meridional/vertical structure.
- `cosp_histogram_plot/`: Synthetic COSP cloud histogram. Checks histogram panel layout, 3x3 cloud-class summaries, and mixed-sign test-minus-reference differences.

## Baseline Files

- `*.png`: Full three-panel figure.
- `*.2.png`: Cropped difference panel used by the targeted regression suite.
