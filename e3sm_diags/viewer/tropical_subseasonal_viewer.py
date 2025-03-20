import os

from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.viewer.core_viewer import OutputViewer
from e3sm_diags.viewer.default_viewer import create_metadata
from e3sm_diags.viewer.utils import add_header, h1_to_h3

logger = _setup_child_logger(__name__)


def create_viewer(root_dir, parameters):
    """
    Given a set of parameters for the enso_diags set,
    create a single webpage.

    Return the title and url for this page.
    """
    viewer = OutputViewer(path=root_dir)

    # The name that's displayed on the viewer.
    display_name = "Tropical Subseasonal Variability Diagnostics"
    set_name = "tropical_subseasonal"
    # The title of the colums on the webpage.
    # Appears in the second and third columns of the bolded rows.
    cols = ["Description", "Plot"]
    viewer.add_page(display_name, short_name=set_name, columns=cols)
    for param in parameters:
        # for plot_type in ["Wavenumber Frequency", "Lag correlation"]:
        # for plot_type in ["Wavenumber Frequency"]:
        # Appears in the first column of the bolded rows.
        plot_type = param.case_id

        for var in param.variables:
            viewer.add_group(f"{plot_type.capitalize()} - {var}")
            print("var,var", var)
            # Appears in the first column of the non-bolded rows.
            # This should use param.case_id to match the output_dir determined by
            # get_output_dir in e3sm_diags/plot/cartopy/enso_diags_plot.py.
            # Otherwise, the plot image and the plot HTML file will have URLs
            # differing in the final directory name.
            formatted_files = []
            for spec_type in [
                "norm_sym",
                "norm_sym_zoom",
                "norm_asy",
                "norm_asy_zoom",
                "raw_sym",
                "raw_asy",
                "background",
            ]:
                viewer.add_row(f"{var} {spec_type} {param.ref_name}")
                # Adding the description for this var to the current row.
                # This was obtained and stored in the driver for this plotset.
                # Appears in the second column of the non-bolded rows.
                # viewer.add_col(param.viewer_descr[var])
                descrip = f"{spec_type} power spectral for {var} 15N-15S"
                if "zoom" in spec_type:
                    descrip = f"{descrip} zoom in for MJO"
                viewer.add_col(descrip)
                # Link to an html version of the plot png file.
                # Appears in the third column of the non-bolded rows.
                ext = param.output_format[0]
                relative_path = os.path.join(
                    "..", set_name, param.case_id, f"{var}_{spec_type}_15N-15S"
                )
                if param.save_netcdf:
                    nc_files = []
                    for ifile in range(3):
                        nc_file_path = "{}_{}.nc".format(relative_path, ifile)
                        nc_files.append(nc_file_path)
                    formatted_files = [{"url": f, "title": f} for f in nc_files]
                image_relative_path = "{}.{}".format(relative_path, ext)

                # image_relative_path = f'{var}_{spec_type}_15N-15N.png'
                viewer.add_col(
                    image_relative_path,
                    is_file=True,
                    title="Plot",
                    other_files=formatted_files,
                    meta=create_metadata(param),
                )

    url = viewer.generate_page()
    add_header(root_dir, os.path.join(root_dir, url), parameters)
    h1_to_h3(os.path.join(root_dir, url))

    return display_name, url
