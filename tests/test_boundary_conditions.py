"""
Tests for the boundary condition classes and the boundary handling of the solvers.
"""

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from acinus_diffusion.boundary_conditions import DirichletBC, NeumannBC, RobinBC  # noqa: E402
from acinus_diffusion.geometry import create_deformed_domain  # noqa: E402
from acinus_diffusion.quasistationary import animate_solution  # noqa: E402
from acinus_diffusion.stationary import assemble_system, grid_spacing  # noqa: E402


def test_bc_classes_match_assembled_rows():
    """The Robin and Dirichlet rows written by the classes equal the solver's."""
    N, M, L, lam = 6, 5, 0.01, 0.05
    A, b = assemble_system(N, M, L, 2.0, lam)
    A = A.toarray()
    _, dy = grid_spacing(N, M, L)

    bottom = list(range(N))                      # j = 0
    top = [(M - 1) * N + i for i in range(N)]    # j = M - 1
    matrix, rhs = np.zeros_like(A), np.zeros_like(b)
    RobinBC(alpha=1 / lam).apply(matrix, rhs, bottom, [k + N for k in bottom], dy)
    DirichletBC(2.0).apply(matrix, rhs, top)

    rows = bottom + top
    np.testing.assert_allclose(matrix[rows], A[rows])
    np.testing.assert_allclose(rhs[rows], b[rows])


@pytest.mark.parametrize("bc", [NeumannBC(), RobinBC(alpha=1.0)])
def test_bc_requires_neighbours_and_spacing(bc):
    matrix, rhs = np.zeros((4, 4)), np.zeros(4)
    with pytest.raises(ValueError):
        bc.apply(matrix, rhs, [0, 1])
    with pytest.raises(ValueError):
        bc.apply(matrix, rhs, [0, 1], [2], h=0.1)


def test_full_deformation_keeps_exchange_surfaces():
    """An inscribed ellipse must not cut the Dirichlet or Robin rows."""
    _, _, mask = create_deformed_domain(41, 41, 0.01, 0.01, 1.0)
    assert mask[0, :].all() and mask[-1, :].all()
    assert not mask.all()


def test_animation_restarts_time_series_on_repeat():
    anim = animate_solution(N=11, M=11, duration=1, fps=5)
    for _ in range(2):
        for frame in range(5):
            anim._func(frame)
    times = list(anim._fig.axes[1].lines[0].get_xdata())
    assert times == [0.0, 0.2, 0.4, 0.6, 0.8]
