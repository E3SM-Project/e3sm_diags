import os

from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.viewer.core_viewer import OutputViewer
from e3sm_diags.viewer.default_viewer import create_metadata
from e3sm_diags.viewer.utils import add_header, h1_to_h3

logger = _setup_child_logger(__name__)


def create_viewer(root_dir, parameters):
    """
    Given a set of parameters for the precip_pdf set,
    create a single webpage.

    Return the title and url for this page.
    """
    viewer = OutputViewer(path=root_dir)

    # The name that's displayed on the viewer
    display_name = "Precipitation PDF Diagnostics"
    set_name = "precip_pdf"
    # The title of the columns on the webpage
    # Appears in the second and third columns of the bolded rows
    cols = ["Description", "Plot"]
    viewer.add_page(display_name, short_name=set_name, columns=cols)

    for param in parameters:
        # Appears in the first column of the bolded rows
        plot_type = param.case_id

        for var in param.variables:
            viewer.add_group(f"{plot_type.capitalize()} - {var}")

            # Loop through regions
            for region in param.regions:
                # Include ref_name in row title to distinguish different reference datasets
                viewer.add_row(f"{var} PDF {region} vs {param.ref_name}")

                # Adding the description for this var to the current row
                # Appears in the second column of the non-bolded rows
                descrip = (
                    f"Precipitation PDF for {var} ({region}) vs {param.reference_name}"
                )

                viewer.add_col(descrip)

                # Link to an html version of the plot png file
                # Appears in the third column of the non-bolded rows
                ext = param.output_format[0]
                # Include ref_name in path to match updated filename format
                relative_path = os.path.join(
                    "..",
                    set_name,
                    param.case_id,
                    f"{var}_PDF_{region}_{param.ref_name}",
                )

                formatted_files = []
                if param.save_netcdf:
                    nc_file_path = f"{relative_path}.nc"
                    formatted_files = [{"url": nc_file_path, "title": nc_file_path}]

                image_relative_path = f"{relative_path}.{ext}"

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
