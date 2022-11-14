import os
import shutil

from cdp.cdp_viewer import OutputViewer

from .lat_lon_viewer import _cvs_to_html
from .utils import add_header, h1_to_h3


def create_viewer(root_dir, parameters):
    """
    Given a set of parameters for a the diff_diags set,
    create a single webpage.

    Return the title and url for this page.
    """
    # The name that's displayed on the viewer.
    display_name = "Aerosol Budget Tables"
    set_name = "aerosol_budget"
    csv_table_dir = os.path.join(
        root_dir, "..", f"{set_name}"
    )  # output_dir/viewer/table-data

    table_dir = os.path.join(
        root_dir, f"{set_name}-table"
    )  # output_dir/viewer/table-data
    if not os.path.exists(table_dir):
        os.mkdir(table_dir)
    seasons = parameters[0].seasons
    test_name = (
        parameters[0].short_test_name
        if parameters[0].short_test_name
        else parameters[0].test_name
    )
    ref_name = (
        parameters[0].short_ref_name
        if parameters[0].short_ref_name
        else parameters[0].ref_name
    )

    viewer = OutputViewer(path=root_dir)

    # The title of the colums on the webpage.
    viewer.add_page(display_name, columns=seasons)
    viewer.add_group("All-Species")
    viewer.add_row("All-Species")

    for season in seasons:
        csv_path = os.path.join(
            csv_table_dir, f"{parameters[0].test_name}-{season}-budget-table.csv"
        )
        shutil.copy(csv_path, table_dir)

        csv_path = os.path.join(
            table_dir, f"{parameters[0].test_name}-{season}-budget-table.csv"
        )
        print("csv_path", csv_path)
        html_path = _cvs_to_html(csv_path, season, test_name, ref_name)

        print("html_path", html_path)
        print("current working dir", os.getcwd())
        # change to relative path i.g. this: viewer/aerosol_budget/ANN_metrics_table.html
        html_path = "/".join(html_path.split("/")[-3:])
        print("html_path", html_path)
        viewer.add_col(html_path, is_file=True, title=season)

    url = viewer.generate_page()
    add_header(root_dir, os.path.join(root_dir, url), parameters)
    h1_to_h3(os.path.join(root_dir, url))

    return display_name, url
