# pylint: disable-msg=too-many-arguments

"""Definition of the Maxwell boundary operators."""


def electric_field(space,
                   wave_number,
                   label="EFIE", symmetry='no_symmetry',
                   parameters=None, use_slp=False):
    """Return the electric field boundary operator."""

    import bempp
    from bempp.core.operators.boundary.maxwell import electric_field_ext
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly.boundary_operator import BoundaryOperator
    from bempp.api.assembly import LocalBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractIntegralOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractLocalOperator

    if parameters is None:
        parameters = bempp.api.global_parameters

    if not use_slp:
        return ElementaryBoundaryOperator( \
                ElementaryAbstractIntegralOperator(
            electric_field_ext(parameters, space._impl, space._impl, space._impl,
                               wave_number, "", symmetry)),
            parameters=parameters, label=label)
    else:

        if not isinstance(use_slp, BoundaryOperator):

            new_space = space.discontinuous_space
            slp = bempp.api.operators.boundary.helmholtz.single_layer(new_space, new_space, new_space, wave_number,
                                                                  parameters=parameters)
        else:
            slp = use_slp

        test_local_ops = []
        trial_local_ops = []

        from bempp.api.assembly.boundary_operator import CompoundBoundaryOperator
        from bempp.core.operators.boundary.sparse import vector_value_times_scalar_ext
        from bempp.core.operators.boundary.sparse import div_times_scalar_ext

        kappa = -1.j * wave_number

        for index in range(3):
            # Definition of range_ does not matter in next operator
            test_local_op = LocalBoundaryOperator(ElementaryAbstractLocalOperator(
                vector_value_times_scalar_ext(slp.dual_to_range._impl, space._impl, space._impl, index)),
                    label='VECTOR_VALUE')
            test_local_ops.append(test_local_op)
            trial_local_ops.append(test_local_op.transpose(space))  # Range parameter arbitrary

        term1 = CompoundBoundaryOperator(test_local_ops, kappa * slp, trial_local_ops, label=label+"_term1")

        test_local_ops = []
        trial_local_ops = []

        div_op = LocalBoundaryOperator(ElementaryAbstractLocalOperator(div_times_scalar_ext(slp.dual_to_range._impl, space._impl, space._impl)),
            label='DIV')
        div_op_transpose = div_op.transpose(space) # Range space does not matter

        term2 = CompoundBoundaryOperator([div_op], (1. / kappa) * slp,
                                         [div_op_transpose], label=label+"_term2")

        return term1 + term2


def magnetic_field(space,
                   wave_number,
                   label="MFIE", symmetry='no_symmetry',
                   parameters=None):
    """Return the magnetic field boundary operator."""

    import bempp
    from bempp.core.operators.boundary.maxwell import magnetic_field_ext
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractIntegralOperator

    if parameters is None:
        parameters = bempp.api.global_parameters

    return ElementaryBoundaryOperator( \
            ElementaryAbstractIntegralOperator(
        magnetic_field_ext(parameters, space._impl, space._impl, space._impl,
                           wave_number, "", symmetry)),
        parameters=parameters, label=label)
