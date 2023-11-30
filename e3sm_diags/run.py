import copy
from itertools import chain
from typing import List

import e3sm_diags  # noqa: F401
from e3sm_diags.e3sm_diags_driver import get_default_diags_path, main
from e3sm_diags.logger import custom_logger, move_log_to_prov_dir
from e3sm_diags.parameter import SET_TO_PARAMETERS
from e3sm_diags.parameter.core_parameter import DEFAULT_SETS, CoreParameter
from e3sm_diags.parser.core_parser import CoreParser

logger = custom_logger(__name__)


class Run:
    """
    Used to run diagnostics.
    A class is needed because we often need to store some
    state regarding what we need to run, like the sets selected.
    """

    def __init__(self):
        self.parser = CoreParser()

        # The list of sets to run using parameter objects.
        self.sets_to_run = []

    @property
    def has_cfg_file_arg(self):
        """A property to check if `-d/--diags` was set to a `.cfg` filepath.

        Returns
        -------
        bool
            True if list contains more than one path, else False.
        """
        args = self.parser.view_args()

        return len(args.other_parameters) > 0

    def run_diags(
        self, parameters: List[CoreParameter], debug: bool = False
    ) -> List[CoreParameter]:
        """Run a set of diagnostics with a list of parameters.

        Parameters
        ----------
        parameters : List[CoreParameter]
            A list of parameters.
        debug : bool, optional
            Run in debug mode or not, by default False.

              * If True, only the parameters passed via ``parameters`` will be
                run. The sets to run are based on the sets defined by the
                parameters. This makes it easy to debug a few sets.
              * If False, run all sets using the list of parameters passed in
                this function and parameters defined in a .cfg file (if
                defined), or use the .cfg file(s) for default diagnostics. This
                is the default option.

        Returns
        -------
        List[CoreParameter]
            A list of parameter objects with their results.

        Raises
        ------
        RuntimeError
            If a diagnostic run using a parameter fails for any reason.
        """
        params = self.get_run_parameters(parameters, debug)

        if params is None or len(params) == 0:
            raise RuntimeError(
                "No parameters we able to be extracted. Please "
                "check the parameters you defined."
            )

        try:
            params_results = main(params)
        except Exception:
            logger.exception("Error traceback:", exc_info=True)

        move_log_to_prov_dir(params_results[0].results_dir)

        return params_results

    def get_run_parameters(self, parameters: List[CoreParameter], debug: bool = False):
        """
        Based on sets_to_run and the list of parameters, get the final list of
        paremeters to run the diags on.
        """
        self._validate_parameters(parameters)

        # FIXME: This line produces some unintended side-effects. For example,
        # let's say we have two objects: 1. CoreParameter, 2. ZonalMean2DParameter.
        # If object 1 has `plevs=[200]`, this value will get copied to object 2.
        # Object 2 has a check to make sure plevs has more than 1 value. This
        # breaks the diagnostic run as a result. The workaround is to loop
        # over `run_diags()` function and run one parameter at a time.
        self._add_parent_attrs_to_children(parameters)

        if not debug:
            run_params = self._get_cfg_parameters(parameters)
        else:
            run_params = self._get_debug_parameters(parameters)

        self.parser.check_values_of_params(run_params)

        return run_params

    def _validate_parameters(self, parameters: List[CoreParameter]):
        if parameters is None or not isinstance(parameters, list):
            raise RuntimeError("You must pass in a list of parameter objects.")

        param_types_list = [
            p.__class__ for p in parameters if p.__class__ != CoreParameter
        ]
        param_types_set = set(param_types_list)

        if len(param_types_set) != len(param_types_list):
            raise RuntimeError(
                "You passed in two or more non-CoreParameter objects of the same type."
            )

    def _get_cfg_parameters(
        self, parameters: List[CoreParameter]
    ) -> List[CoreParameter]:
        """Get the run parameters using all sets and a .cfg file if passed.

        Parameters
        ----------
        parameters : List[CoreParameter]
            A list of parameter objects.

        Returns
        -------
        List[CoreParameter]
            A list of parameter objects including ones defined in a .cfg file.
            Any non-CoreParameter objects will be replaced by a sub-class
            based on the set (``SETS_TO_PARAMETERS``).
        """
        run_params = []

        if self.has_cfg_file_arg:
            cfg_params = self._get_custom_params_from_cfg_file()
        else:
            run_type = parameters[0].run_type
            cfg_params = self._get_default_params_from_cfg_file(run_type)

        if len(self.sets_to_run) == 0:
            self.sets_to_run = DEFAULT_SETS

        for set_name in self.sets_to_run:
            param = self._get_instance_of_param_class(
                SET_TO_PARAMETERS[set_name], parameters
            )

            # Since each parameter will have lots of default values, we want to
            # remove them. Otherwise when calling get_parameters(), these
            # default values will take precedence over values defined in
            # other_params.
            self._remove_attrs_with_default_values(param)
            param.sets = [set_name]

            # # FIXME: Make a deep copy of cfg_params because there is some
            # buggy code in this method that changes parameter attributes in
            # place, which affects downstream operations. The original
            # cfg_params needs to be perserved for each iteration of this
            # for loop.
            params = self.parser.get_parameters(
                orig_parameters=param,
                other_parameters=copy.deepcopy(cfg_params),
                cmd_default_vars=False,
                argparse_vals_only=False,
            )

            # Makes sure that any parameters that are selectors will be in param.
            self._add_attrs_with_default_values(param)

            # The select() call in get_parameters() was made for the original
            # command-line way of using CDP. We just call it manually with the
            # parameter object param.
            params = self.parser.select(param, params)

            run_params.extend(params)

        return run_params

    def _get_custom_params_from_cfg_file(self) -> List[CoreParameter]:
        """Get parameters using the cfg file set by `-d`/`--diags`.

        Returns
        -------
        List[CoreParameter]
            A list of parameter objects.
        """
        params = self.parser.get_cfg_parameters(argparse_vals_only=False)

        params_final = self._convert_params_to_subclass(params)

        return params_final

    def _get_default_params_from_cfg_file(self, run_type: str) -> List[CoreParameter]:
        """Get parameters using the default diagnostic .cfg file(s).

        Parameters
        ----------
        run_type : str
            The run type used to check for which .cfg file(s) to reference.

        Returns
        -------
        List[CoreParameter]
            A list of parameter objects.
        """
        # Get the paths to the default diags .cfg file(s).
        paths = []
        for set_name in self.sets_to_run:
            path = get_default_diags_path(set_name, run_type, False)
            paths.append(path)

        # Convert the .cfg file(s) to parameter objects.
        params = self.parser.get_cfg_parameters(
            files_to_open=paths, argparse_vals_only=False
        )

        # Update parameter objects using subclass with default values.
        params_final = self._convert_params_to_subclass(params)

        return params_final

    def _convert_params_to_subclass(
        self,
        parameters: List[CoreParameter],
    ) -> List[CoreParameter]:
        new_params: List[CoreParameter] = []

        # For each of the params, add in the default values using the parameter
        # classes in SET_TO_PARAMETERS.
        for param in parameters:
            set_key = param.sets[0]
            new_param = SET_TO_PARAMETERS[set_key]() + param

            new_params.append(new_param)

        return new_params

    def _get_debug_parameters(
        self, parameters: List[CoreParameter]
    ) -> List[CoreParameter]:
        """Get the parameters explicitly defined in the Python script only.

        This method replaces CoreParameter objects with the related sub-class
        for the specified set.

        Parameters
        ----------
        parameters : List[CoreParameter]
            A list of parameter objects.

        Returns
        -------
        List[CoreParameter]
            A list of parameter objects. Any non-CoreParameter objects will be
            replaced by a sub-class based on the set (``SETS_TO_PARAMETERS``).
        """
        debug_params = []

        if len(self.sets_to_run) == 0:
            sets_to_run = [param.sets for param in parameters]
            self.sets_to_run = list(chain.from_iterable(sets_to_run))

        for set_name in self.sets_to_run:
            # For each of the set_names, get the corresponding parameter.
            api_param = self._get_instance_of_param_class(
                SET_TO_PARAMETERS[set_name], parameters
            )

            # Since each parameter will have lots of default values, we want to remove them.
            # Otherwise when calling get_parameters(), these default values
            # will take precedence over values defined in other_params.
            self._remove_attrs_with_default_values(api_param)
            api_param.sets = [set_name]

            # Makes sure that any parameters that are selectors will be in param.
            self._add_attrs_with_default_values(api_param)

            debug_params.append(api_param)

        self.parser.check_values_of_params(debug_params)

        return debug_params

    def _add_parent_attrs_to_children(self, parameters):
        """
        For any parameter class that's inherited from another, copy
        the attributes of the parent to the child.

        Ex: If the user wants to run set-specific parameters for
        'zonal_mean_2d', they'd pass in a ZonalMean2dParameter
        and a CoreParameter.
        But most of the important parameters are in CoreParameter,
        so copy them over to the ZonalMean2dParameter.
        """

        def get_parent(param):
            """
            From parameters, get any object that's a parent
            type to param.

            Ex: CoreParameter object is a parent of AreaMeanTimeSeriesParameter object
            """
            try:
                parent_class = param.__class__.__mro__[1]
                parent = self._get_instance_of_param_class(parent_class, parameters)
            except RuntimeError:
                parent = None

            return parent

        for i in range(len(parameters)):
            parent = get_parent(parameters[i])
            # Make sure that the new object is actually a parent.
            if not parent or type(parent) == type(parameters[i]):  # noqa:E721
                continue

            # Otherwise, add the the parent's attributes.
            # Since we're modifying this parent object (by
            # removing the default values before addition)
            # make a deepcopy first.
            parent = copy.deepcopy(parent)
            # find attributes that are not defaults

            nondefault_param_parent = self._find_attrs_with_nondefault_values(parent)
            nondefault_param_child = self._find_attrs_with_nondefault_values(
                parameters[i]
            )

            self._remove_attrs_with_default_values(parent)

            # Simply copy over all attribute from parent to children
            # parameters[i] += parent

            for attr in dir(parent):
                if not attr.startswith("_") and not hasattr(parameters[i], attr):
                    # This attr of parent is a user-defined one and does not
                    # already exist in the parameters[i] parameter object.
                    attr_value = getattr(parent, attr)
                    setattr(parameters[i], attr, attr_value)

            logger.info(
                list(set(nondefault_param_parent) - set(nondefault_param_child))
            )
            for attr in list(
                set(nondefault_param_parent) - set(nondefault_param_child)
            ):
                # 'seasons' is a corner case that don't need to get in to none-core sets, Ex. area mean time series
                if attr != "seasons":
                    attr_value = getattr(parent, attr)
                    setattr(parameters[i], attr, attr_value)

    def _add_attrs_with_default_values(self, param):
        """
        In the param, add any missing parameters
        with their default value.
        """
        new_instance = param.__class__()
        for attr in dir(new_instance):
            # Ignore any of the hidden attributes.
            if attr.startswith("_"):
                continue

            if not hasattr(param, attr):
                val = getattr(new_instance, attr)
                setattr(param, attr, val)

    def _remove_attrs_with_default_values(self, param):
        """
        In the param, remove any parameters that
        have their default value.
        """
        new_instance = param.__class__()
        for attr in dir(param):
            # Ignore any of the hidden attributes.
            if attr.startswith("_"):
                continue

            if hasattr(new_instance, attr) and getattr(new_instance, attr) == getattr(
                param, attr
            ):
                delattr(param, attr)

    def _find_attrs_with_nondefault_values(self, param):
        """
        In the param, find any parameters that
        have nondefault value.
        """
        nondefault_attr = []
        new_instance = param.__class__()
        for attr in dir(param):
            # Ignore any of the hidden attributes.
            if attr.startswith("_"):
                continue

            if hasattr(new_instance, attr) and getattr(new_instance, attr) != getattr(
                param, attr
            ):  # This is only valid when the attr values are lists not numpy array
                nondefault_attr.append(attr)
        return nondefault_attr

    def _get_instance_of_param_class(self, cls, parameters):
        """
        In the list of parameters, get the class for
        the parameter object corresponding to cls.
        """
        # Get the Method Resolution Order (MRO) for this class.
        # So get the list of the classes in the inheritance ordering.

        # Ex: For the 'zonal_mean_2d' set, ZonalMean2dParameter is
        # the parameter for it.
        # But if a user doesn't want to modify the set-specific
        # parameters for 'zonal_mean_2d', they can just pass in
        # a single CoreParameter object to run_diags().
        # Using this, we handle this use-case.
        class_types = cls.__mro__

        for cls_type in class_types:
            for p in parameters:
                if isinstance(p, cls_type):
                    return p

        msg = "There's weren't any class of types {} in your parameters."
        raise RuntimeError(msg.format(class_types))


runner = Run()
