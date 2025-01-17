import os
from typing import Dict, List

from cdp.cdp_viewer import OutputViewer

from e3sm_diags.logger import custom_logger

logger = custom_logger(__name__)


def create_viewer(root_dir, parameters):
    """
    Given a set of parameters for the enso_diags set,
    create a single webpage.

    Return the title and url for this page.
    """
    viewer = OutputViewer(path=root_dir)

    # The name that's displayed on the viewer.
    display_name = "Tropical Subseasonal Variability Diagnostics"
    set_name = "tropical_subseasonal_diags"
    # The title of the colums on the webpage.
    # Appears in the second and third columns of the bolded rows.
    cols = ["Description", "Plot"]
    viewer.add_page(display_name, short_name=set_name, columns=cols)
    for plot_type in ["Wavenumber Frequency", "Lag correlation"]:
        # Appears in the first column of the bolded rows.
        viewer.add_group(plot_type.capitalize())

        for var in parameters.variables:
            # Appears in the first column of the non-bolded rows.
            # This should use param.case_id to match the output_dir determined by
            # get_output_dir in e3sm_diags/plot/cartopy/enso_diags_plot.py.
            # Otherwise, the plot image and the plot HTML file will have URLs
            # differing in the final directory name.
            for spec_type in ["norm_sym", "norm_sym_zoom", "norm_asy", "norm_asy_zoom", "raw_sym", "raw_asy", "background"]:
                viewer.add_row(f'{var} {spec_type} ref_name')
                # Adding the description for this var to the current row.
                # This was obtained and stored in the driver for this plotset.
                # Appears in the second column of the non-bolded rows.
                #viewer.add_col(param.viewer_descr[var])
                viewer.add_col(f'Long description for var')
                # Link to an html version of the plot png file.
                # Appears in the third column of the non-bolded rows.
                image_relative_path = f'{var}_{spec_type}_15N-15N.png'
                viewer.add_col(
                    image_relative_path,
                    is_file=True,
                    title="Plot",
                    #other_files=formatted_files,
                    #meta=create_metadata(param),
                )

    url = viewer.generate_page()

    return display_name, url
