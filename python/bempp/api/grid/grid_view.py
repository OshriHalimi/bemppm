"""Implementation of a view onto a grid."""


class GridView(object):
    """Provides interfaces to query the entities of a grid."""

    def __init__(self, impl):
        """Constructor. Not to be called directly."""
        self._impl = impl
        self._index_set = None
        self._reverse_element_map = None
        self._elements = None
        self._vertices = None
        self._edges = None
        self._connectivity = None
        self._maximum_diameter = None
        self._minimum_diameter = None

    def _create_connectivity_matrices(self):
        """
        Create connectivity matrices between entities.

        This function creates sparse matrices that map vertices to
        elements and edges to elements.

        """
        import numpy as np
        from scipy.sparse import csc_matrix

        entity_count = {0: self.entity_count(0),
                        1: self.entity_count(1),
                        2: self.entity_count(2)}

        edge_indices = []
        edge_element_indices = []
        vertex_indices = []
        vertex_element_indices = []

        ind_set = self.index_set()

        for element in self.entity_iterator(0):

            index = ind_set.entity_index(element)

            for i in range(3):
                # Loop over vertices
                vertex_index = ind_set.sub_entity_index(element, i, 2)
                vertex_indices.append(vertex_index)
                vertex_element_indices.append(index)

            for j in range(3):
                # Loop over edges
                edge_index = ind_set.sub_entity_index(element, j, 1)
                edge_indices.append(edge_index)
                edge_element_indices.append(index)

        # Now create the matrices

        self._connectivity = {}

        self._connectivity['vertices'] = csc_matrix(
            (np.ones(len(vertex_indices), dtype='int64'),
             (vertex_indices, vertex_element_indices)),
            shape=(entity_count[2], entity_count[0]),
            dtype=np.int64)

        self._connectivity['edges'] = csc_matrix(
            (np.ones(len(edge_indices), dtype='int64'),
             (edge_indices, edge_element_indices)),
            shape=(entity_count[1], entity_count[0]),
            dtype=np.int64)

    def entity_count(self, codimension):
        """Return the number of entities of a given codimension."""
        return self._impl.entity_count(codimension)

    def index_set(self):
        """Return an IndexSet object for the GridView."""
        if self._index_set is None:
            from bempp.api.grid.index_set import IndexSet
            self._index_set = IndexSet(self._impl.index_set())
        return self._index_set

    def minimum_element_diameter(self):
        """Return minimum element diameter."""
        
        if self._minimum_diameter is None:
            self._minimum_diameter = self._impl.minimum_element_diameter()
        
        return self._minimum_diameter

    def maximum_element_diameter(self):
        """Return maximum element diameter."""
        
        if self._maximum_diameter is None:
            self._maximum_diameter = self._impl.maximum_element_diameter()
        
        return self._maximum_diameter



    def entity_iterator(self, codimension):
        """Return an entity iterator for a given codimension."""
        from bempp.api.grid.entity_iterator import EntityIterator
        if codimension == 0:
            if self._elements is None:
                self._elements = list(
                    EntityIterator(codimension, self._impl.entity_iterator(
                        codimension)))
            return iter(self._elements)
        if codimension == 1:
            if self._edges is None:
                self._edges = list(
                    EntityIterator(codimension, self._impl.entity_iterator(
                        codimension)))
            return iter(self._edges)
        if codimension == 2:
            if self._vertices is None:
                self._vertices = list(
                    EntityIterator(codimension, self._impl.entity_iterator(
                        codimension)))
            return iter(self._vertices)
        raise ValueError("Unknown codimension.")

    def element_from_index(self, index):
        """Map a given index to the associated element."""
        if self._reverse_element_map is None:
            self._reverse_element_map = [None] * self.entity_count(0)
            ind_set = self.index_set()
            for element in self.entity_iterator(0):
                elem_index = ind_set.entity_index(element)
                self._reverse_element_map[elem_index] = element

        return self._reverse_element_map[index]

    @property
    def vertex_to_element_matrix(self):
        """
        Return the vertex to element matrix.

        The vertex to element matrix is a sparse matrix A, where A[i, j]
        is 1 if vertex i is associated with element j, otherwise A[i, j] = 0.

        """
        if self._connectivity is None:
            self._create_connectivity_matrices()
        return self._connectivity['vertices']

    @property
    def edge_to_element_matrix(self):
        """
        Return the edge to element matrix.

        The dge to element matrix is a sparse matrix A, where A[i, j]
        is 1 if edge i is associated with element j, otherwise A[i, j] = 0.

        """
        if self._connectivity is None:
            self._create_connectivity_matrices()
        return self._connectivity['edges']

    @property
    def dim(self):
        """Return the dimension of the grid."""
        return self._impl.dim

    @property
    def dim_world(self):
        """Return the dimension of the space containing the grid."""
        return self._impl.dim_world

    @property
    def vertices(self):
        """Return a (3 x n_vertices) array with the vertices of the grid."""
        return self._impl.vertices

    @property
    def elements(self):
        """Return a (3 x n_elements) array with the elements of the grid."""
        return self._impl.elements

    @property
    def domain_indices(self):
        """Return a list of domain indices."""
        return self._impl.domain_indices
