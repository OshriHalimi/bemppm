# pylint: disable-msg=too-many-arguments

"""Definition of the modified Helmholtz boundary operators."""


def single_layer(domain, range_, dual_to_range,
                 wave_number,
                 label="SLP", symmetry='no_symmetry',
                 parameters=None):
    """Return the single-layer boundary operator."""

    import bempp
    from bempp.core.operators.boundary.modified_helmholtz import single_layer_ext
    from bempp.api.assembly import ElementaryBoundaryOperator

    if parameters is None:
        parameters = bempp.api.global_parameters

    return ElementaryBoundaryOperator( \
        single_layer_ext(parameters, domain, range_,
                         dual_to_range, wave_number, "", symmetry),
        parameters=parameters, label=label)


def double_layer(domain, range_, dual_to_range,
                 wave_number,
                 label="DLP", symmetry='no_symmetry',
                 parameters=None):
    """Return the double-layer boundary operator."""

    import bempp
    from bempp.core.operators.boundary.modified_helmholtz import double_layer_ext
    from bempp.api.assembly import ElementaryBoundaryOperator

    if parameters is None:
        parameters = bempp.api.global_parameters

    return ElementaryBoundaryOperator( \
        double_layer_ext(parameters, domain, range_,
                         dual_to_range, wave_number, "", symmetry),
        parameters=parameters, label=label)


def adjoint_double_layer(domain, range_, dual_to_range,
                         wave_number,
                         label="ADJ_DLP", symmetry='no_symmetry',
                         parameters=None):
    """Return the adjoint double-layer boundary operator."""

    import bempp
    from bempp.core.operators.boundary.modified_helmholtz import adjoint_double_layer_ext
    from bempp.api.assembly import ElementaryBoundaryOperator

    if parameters is None:
        parameters = bempp.api.global_parameters

    return ElementaryBoundaryOperator( \
        adjoint_double_layer_ext(parameters, domain, range_,
                                 dual_to_range, wave_number, "", symmetry),
        parameters=parameters, label=label)


def hypersingular(domain, range_, dual_to_range, wave_number,
                  label="HYP", symmetry='no_symmetry',
                  parameters=None, use_slp=False):
    """Return the hypersingular boundary operator."""

    import bempp
    from bempp.core.operators.boundary.modified_helmholtz import hypersingular_ext
    from bempp.api.assembly.boundary_operator import BoundaryOperator
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly import LocalBoundaryOperator

    if parameters is None:
        parameters = bempp.api.global_parameters

    if domain != dual_to_range and use_slp:
        print("Compound assembly based on slp operator requires 'domain' and 'dual_to_range' space to be identical." +
              " Switching to standard assembly.")
        use_slp = False

    if not use_slp:
        return ElementaryBoundaryOperator( \
            hypersingular_ext(parameters, domain, range_,
                              dual_to_range, wave_number, "", symmetry),
            parameters=parameters, label=label)
    else:

        if not isinstance(use_slp, BoundaryOperator):
            new_domain = domain.discontinuous_space
            new_dual_to_range = dual_to_range.discontinuous_space
            slp = single_layer(new_domain, range_, new_dual_to_range, wave_number, parameters=parameters)
        else:
            slp = use_slp

        # Now generate the compound operator

        test_local_ops = []
        trial_local_ops = []

        from bempp.api.assembly.boundary_operator import CompoundBoundaryOperator
        from bempp.core.operators.boundary.sparse import curl_value_ext
        from bempp.core.operators.boundary.sparse import value_times_normal_ext

        for index in range(3):
            # Definition of range_ does not matter in next operator
            test_local_op = LocalBoundaryOperator(curl_value_ext(slp.dual_to_range, range_, dual_to_range, index),
                label='CURL')
            test_local_ops.append(test_local_op)
            trial_local_ops.append(test_local_op.transpose(range_))  # Range parameter arbitrary

        term1 = CompoundBoundaryOperator(test_local_ops, slp, trial_local_ops, label=label + "_term1")

        test_local_ops = []
        trial_local_ops = []

        for index in range(3):
            # Definition of range_ does not matter in next operator
            test_local_op = LocalBoundaryOperator(
                value_times_normal_ext(slp.dual_to_range, range_, dual_to_range, index),
                    label='VALUE_TIMES_NORMAL')
            test_local_ops.append(test_local_op)
            trial_local_ops.append(test_local_op.transpose(range_))  # Range parameter arbitrary

        term2 = (wave_number * wave_number) * CompoundBoundaryOperator(test_local_ops, slp, trial_local_ops,
                                                                       label=label + "_term2")

        return term1 + term2

def multitrace_operator(grid, wave_number, parameters=None):

    def op(operator):
        if operator==hypersingular:
            def op_impl(domain, range_, dual_to_range, label="HYP", symmetry="no_symmetry",
                        parameters=None, use_slp=False):
                return hypersingular(domain, range_, dual_to_range, wave_number, label, symmetry, parameters,
                                     use_slp)
            return op_impl
        else:
            import inspect
            defaults = inspect.getargspec(operator).defaults
            def op_impl(domain, range_, dual_to_range, label=defaults[0], symmetry=defaults[1],
                        parameters=None):
                return operator(domain, range_, dual_to_range, wave_number, label, symmetry, parameters)
            return op_impl

    from bempp.api.operators.boundary import _common
    return _common.multitrace_operator_impl(grid, op(single_layer), op(double_layer),
                                            op(hypersingular), parameters)


def interior_calderon_projector(grid, wave_number, parameters=None):

    from .sparse import multitrace_identity

    return .5 * multitrace_identity(grid, parameters) + multitrace_operator(grid, wave_number, parameters)

def exterior_calderon_projector(grid, wave_number, parameters=None):

    from .sparse import multitrace_identity

    return .5 * multitrace_identity(grid, parameters) - multitrace_operator(grid, wave_number, parameters)
