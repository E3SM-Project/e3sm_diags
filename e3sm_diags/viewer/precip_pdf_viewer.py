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

    # Determine which seasons to display based on parameter settings
    # Check if any parameter has season_subset enabled
    has_season_subset = any(
        getattr(param, "season_subset", True) for param in parameters
    )

    if has_season_subset:
        # Show all seasons in separate columns
        seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
        cols = ["Description"] + seasons
    else:
        # Only show annual (all months)
        seasons = ["ANN"]
        cols = ["Description", "Plot"]

    viewer.add_page(display_name, short_name=set_name, columns=cols)

    for param in parameters:
        # Appears in the first column of the bolded rows
        plot_type = param.case_id

        # Normalize ref_name to list for consistent handling
        if isinstance(param.ref_name, str):
            ref_names = [param.ref_name] if param.ref_name else []
        else:
            ref_names = param.ref_name

        # Create string for display and filename
        ref_names_str = "_".join(ref_names) if ref_names else ""
        ref_names_display = " & ".join(ref_names) if ref_names else ""

        for var in param.variables:
            viewer.add_group(f"{plot_type.capitalize()} - {var}")

            # Loop through regions
            for region in param.regions:
                # Include ref_name in row title to distinguish different reference datasets
                viewer.add_row(f"{var} PDF {region} vs {ref_names_display}")

                # Adding the description for this var to the current row
                # Appears in the second column of the non-bolded rows
                if hasattr(param, "reference_name") and param.reference_name:
                    descrip = f"Precipitation PDF for {var} ({region}) vs {param.reference_name}"
                else:
                    descrip = (
                        f"Precipitation PDF for {var} ({region}) vs {ref_names_display}"
                    )

                viewer.add_col(descrip)

                # Add a column for each season
                ext = param.output_format[0]
                for season in seasons:
                    # Include ref_name and season in path to match updated filename format
                    relative_path = os.path.join(
                        "..",
                        set_name,
                        param.case_id,
                        f"{var}_PDF_{region}_{ref_names_str}_{season}",
                    )

                    formatted_files = []
                    if param.save_netcdf:
                        # Link to the season-specific netCDF cache file
                        nc_file_path = f"{var}_PDF_global_test_{param.test_name}_{param.test_start_yr}-{param.test_end_yr}_{season}.nc"
                        nc_relative_path = os.path.join(
                            "..", set_name, param.case_id, nc_file_path
                        )
                        formatted_files.append(
                            {
                                "url": nc_relative_path,
                                "title": f"Test NetCDF ({season})",
                            }
                        )

                        # Link all reference netCDF files (one for each ref dataset)
                        for ref_name in ref_names:
                            ref_nc_path = f"{var}_PDF_global_ref_{ref_name}_{param.ref_start_yr}-{param.ref_end_yr}_{season}.nc"
                            ref_nc_relative_path = os.path.join(
                                "..", set_name, param.case_id, ref_nc_path
                            )
                            formatted_files.append(
                                {
                                    "url": ref_nc_relative_path,
                                    "title": f"{ref_name} NetCDF ({season})",
                                }
                            )

                    image_relative_path = f"{relative_path}.{ext}"

                    viewer.add_col(
                        image_relative_path,
                        is_file=True,
                        title=season,
                        other_files=formatted_files,
                        meta=create_metadata(param),
                    )

    url = viewer.generate_page()
    add_header(root_dir, os.path.join(root_dir, url), parameters)
    h1_to_h3(os.path.join(root_dir, url))

    return display_name, url
