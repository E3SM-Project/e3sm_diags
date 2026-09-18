"""Checks on the default diagnostic cfgs that ship with e3sm_diags."""

from collections import Counter
from pathlib import Path

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parser.core_parser import CoreParser

DEFAULT_DIAGS_DIR = (
    Path(__file__).resolve().parents[2] / "e3sm_diags" / "driver" / "default_diags"
)


def test_lat_lon_model_vs_obs_blocks_do_not_write_the_same_output():
    """No two blocks may produce the same output file.

    Blocks run in parallel, so two that share an output race: whichever finishes
    last overwrites the other's plot and metrics, and the result changes from
    run to run. A repeated key within a block causes this silently, since the
    last value wins -- e.g. a land block that also says ``regions = ["global"]``
    becomes a second global block.
    """
    parameters = CoreParser().get_cfg_parameters(
        files_to_open=[str(DEFAULT_DIAGS_DIR / "lat_lon_model_vs_obs.cfg")],
        check_values=False,
    )
    # A block that leaves a key unset gets the CoreParameter default at run time.
    defaults = CoreParameter()

    def value(param, name):
        return getattr(param, name, None) or getattr(defaults, name)

    outputs: Counter = Counter()
    for param in parameters:
        for variable in value(param, "variables"):
            for region in value(param, "regions"):
                for season in value(param, "seasons"):
                    for plev in value(param, "plevs") or [None]:
                        key = (
                            value(param, "case_id"),
                            value(param, "ref_name"),
                            variable,
                            region,
                            season,
                            plev,
                        )
                        outputs[key] += 1

    collisions = sorted(str(key) for key, count in outputs.items() if count > 1)
    assert collisions == []
