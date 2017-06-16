"""Validation tests for the FMM"""

from unittest import TestCase
import bempp.api
import numpy as np

#pylint: disable=invalid-name

TOL = 1E-2

class TestFmm(TestCase):
    """Test the FMM for the various potential operators."""

    def test_laplace(self):
        """Test the FMM for Laplace operators."""

        grid = bempp.api.shapes.regular_sphere(6)
        space = bempp.api.function_space(grid, "DP", 0)

        fmm_params = bempp.api.common.global_parameters()
        fmm_params.fmm.expansion_order = 4
        fmm_params.assembly.boundary_operator_assembly_type = 'fmm'

        hmat_params = bempp.api.common.global_parameters()
        hmat_params.assembly.boundary_operator_assembly_type = 'hmat'
        hmat_params.hmat.eps = 1E-4

        from bempp.api.operators.boundary import laplace as mode

        op_types = [mode.single_layer, mode.double_layer, 
                    mode.adjoint_double_layer]
        test_vec = np.random.randn(space.global_dof_count)
        diff = []
        
        for op_type in op_types:
            op_fmm = op_type(space, space, space, parameters=fmm_params)
            op_hmat = op_type(space, space, space, parameters=hmat_params)
            result_fmm = op_fmm.weak_form() * test_vec
            result_hmat = op_hmat.weak_form() * test_vec

            diff.append(np.linalg.norm(result_fmm - result_hmat) / 
                np.linalg.norm(result_hmat))
            print(diff[-1])

        self.assertTrue(np.all(np.array(diff) < TOL) )

    def test_modified_helmholtz(self):
        """Test the FMM for Laplace operators."""

        grid = bempp.api.shapes.regular_sphere(6)
        space = bempp.api.function_space(grid, "DP", 0)

        k = 1j

        fmm_params = bempp.api.common.global_parameters()
        fmm_params.assembly.boundary_operator_assembly_type = 'fmm'
        fmm_params.fmm.expansion_order = 4
 

        hmat_params = bempp.api.common.global_parameters()
        hmat_params.assembly.boundary_operator_assembly_type = 'hmat'
        hmat_params.hmat.eps = 1E-4

        from bempp.api.operators.boundary import modified_helmholtz as mode
        op_types = [mode.single_layer, mode.double_layer, 

                    mode.adjoint_double_layer]
        test_vec = np.random.randn(space.global_dof_count)
        diff = []
        
        for op_type in op_types:
            op_fmm = op_type(space, space, space, 1j, parameters=fmm_params)
            op_hmat = op_type(space, space, space, 1j, parameters=hmat_params)
            result_fmm = op_fmm.weak_form() * test_vec
            result_hmat = op_hmat.weak_form() * test_vec

            diff.append(np.linalg.norm(result_fmm - result_hmat) / 
                np.linalg.norm(result_hmat))

        self.assertTrue(np.all(np.array(diff) < TOL) )


if __name__ == "__main__":

    from unittest import main
    main()
                    

        

        
        

        

