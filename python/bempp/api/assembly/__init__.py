"""This module contains basic classes for the assembly of integral operators."""

from .discrete_boundary_operator import GeneralNonlocalDiscreteBoundaryOperator
from .discrete_boundary_operator import DenseDiscreteBoundaryOperator
from .discrete_boundary_operator import SparseDiscreteBoundaryOperator
from .discrete_boundary_operator import InverseSparseDiscreteBoundaryOperator
from .discrete_boundary_operator import ZeroDiscreteBoundaryOperator
from .discrete_boundary_operator import DiscreteRankOneOperator
from .discrete_boundary_operator import as_matrix
from .boundary_operator import BoundaryOperator
from .boundary_operator import LocalBoundaryOperator
from .boundary_operator import InverseLocalBoundaryOperator
from .boundary_operator import ElementaryBoundaryOperator
from .boundary_operator import BoundaryOperator
from .boundary_operator import ZeroBoundaryOperator
from .boundary_operator import RankOneBoundaryOperator
from .blocked_operator import BlockedOperator
from .blocked_operator import BlockedDiscreteOperator
from .grid_function import GridFunction
from .assembler import assemble_dense_block
from .potential_operator import PotentialOperator

from . import functors

