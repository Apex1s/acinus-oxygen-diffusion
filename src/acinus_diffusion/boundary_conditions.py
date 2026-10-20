"""
Boundary condition descriptions for the diffusion equation.

Each condition writes the first-order finite-difference equation of a
boundary point into the rows of a (dense or LIL sparse) matrix. ``indices``
are the boundary points and ``neighbor_indices`` the adjacent interior points
along the inward normal, at distance ``h``.

These classes are a standalone toolkit for assembling custom systems. The
solvers in :mod:`acinus_diffusion.stationary` build their matrix directly (in
COO form, which these row-wise writers do not target) and implement the same
equations; ``tests/test_boundary_conditions.py`` checks the two agree.
"""


def _neighbors(indices, neighbor_indices, h):
    if neighbor_indices is None or h is None:
        raise ValueError("this boundary condition needs neighbor_indices and h")
    if len(neighbor_indices) != len(indices):
        raise ValueError("indices and neighbor_indices must have the same length")
    if h <= 0:
        raise ValueError("h must be positive")
    return zip(indices, neighbor_indices)


class BoundaryCondition:
    """Base class for boundary conditions."""

    def __init__(self, value=0.0):
        self.value = value

    def apply(self, matrix, rhs, indices, neighbor_indices=None, h=None):
        raise NotImplementedError("Subclasses must implement apply method")


class DirichletBC(BoundaryCondition):
    """Dirichlet boundary condition: u = value."""

    def apply(self, matrix, rhs, indices, neighbor_indices=None, h=None):
        for idx in indices:
            matrix[idx, idx] = 1
            rhs[idx] = self.value


class NeumannBC(BoundaryCondition):
    """Neumann boundary condition: ∂u/∂n = value (outward normal)."""

    def apply(self, matrix, rhs, indices, neighbor_indices=None, h=None):
        for idx, neighbor_idx in _neighbors(indices, neighbor_indices, h):
            matrix[idx, idx] = 1
            matrix[idx, neighbor_idx] = -1
            rhs[idx] = self.value * h


class RobinBC(BoundaryCondition):
    """
    Robin boundary condition: ∂u/∂n + alpha * u = value (outward normal).

    The capillary condition ∂u/∂n = -u / Λ corresponds to alpha = 1 / Λ and
    value = 0.
    """

    def __init__(self, alpha, value=0.0):
        super().__init__(value)
        self.alpha = alpha

    def apply(self, matrix, rhs, indices, neighbor_indices=None, h=None):
        for idx, neighbor_idx in _neighbors(indices, neighbor_indices, h):
            matrix[idx, idx] = 1 + self.alpha * h
            matrix[idx, neighbor_idx] = -1
            rhs[idx] = self.value * h
