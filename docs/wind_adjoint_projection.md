# Discretely adjoint variational-wind projection

## Purpose

The production variational-wind correction must remove the same discrete
mass imbalance that HICAR diagnoses.  Calibrating a matrix from a separately
discretized divergence and correction preserves operator equivalence, but it
does not make those two operations adjoints.  On steep, large domains that
composition is strongly non-normal and can have negative Rayleigh quotients;
restarted Krylov and local-preconditioner tuning cannot repair that algebra.

This design replaces the composition by a constrained discrete minimization.
The initial implementation uses a diagonal mass-flux energy.  It reproduces
the non-terrain-cross terms of HICAR's existing correction exactly and yields
a symmetric positive-definite Schur complement by construction.

## Discrete contract

Let `q = (u, v, w_grid)` contain staggered face velocities.  Define `B` as the
volume-integrated form of the existing finite-volume divergence:

```text
B q = cell_volume * calc_divergence(q)
```

For cell volume

```text
V_c = dx^2 dz_k jaco_c / mxy
```

the signed face coefficients in `B` are

```text
b_u = dx dz_k jaco_u rho_u / my_u
b_v = dx dz_k jaco_v rho_v / mx_v
b_w = dx^2 jaco_w rho_w / mxy
```

The correction is the solution of

```text
minimize    1/2 delta_q^T M delta_q
subject to  B (q_0 + delta_q) = 0,
```

with a positive face energy `M`.  The optimality equations give

```text
delta_q = -M^-1 B^T lambda
K lambda = B q_0
K = B M^-1 B^T.
```

Consequently,

```text
lambda^T K lambda
  = (B^T lambda)^T M^-1 (B^T lambda) > 0
```

for every nonzero admissible multiplier.  A converged solve also gives the
identity

```text
B (q_0 + delta_q) = B q_0 - K lambda = 0.
```

These are implementation gates, not continuum arguments: the production
halo, boundary, GPU, and precision paths must satisfy them numerically.

## Initial diagonal energy

The reference implementation uses squared mass flux weighted by physical
dual volume:

```text
M_u = 2 rho_u^2 dx^2 dz_k jaco_u / (mx_u my_u)
M_v = 2 rho_v^2 dx^2 dz_k jaco_v / (mx_v my_v)
M_w = 2 rho_w^2 dx^2 dz_w jaco_w^2 / (mxy alpha_w^2).
```

This gives the correction mobilities

```text
b_u / M_u = mx_u / (2 rho_u dx)
b_v / M_v = my_v / (2 rho_v dx)
b_w / M_w = alpha_w^2 / (2 rho_w dz_w jaco_w),
```

which are exactly HICAR's current flat-grid/non-cross coefficients.  Density
cancels from the Schur face conductances:

```text
c_u = b_u^2 / M_u = dz_k jaco_u mx_u / (2 my_u)
c_v = b_v^2 / M_v = dz_k jaco_v my_v / (2 mx_v)
c_w = b_w^2 / M_w = dx^2 alpha_w^2 / (2 dz_w mxy).
```

The production matrix therefore needs only a seven-point weighted-Laplacian
stencil (or on-the-fly face conductances); it must not retain separate
double-precision `B` and `M^-1` arrays at national scale.

## Boundary contract

All physical cells are active.  HICAR's scalar storage has one lateral ghost
layer and two non-physical vertical ghost planes; their multiplier rows are
identity with zero right-hand side.  A zero exterior multiplier allows the
projection to adjust external lateral face fluxes and removes the constant
nullspace.  Ground and top `w_grid` fluxes remain fixed because only the
`nz-1` interior vertical interfaces are correctable.

This is the current production-compatible policy.  A future policy that
fixes lateral face fluxes must instead define an active-row mask, handle the
resulting global compatibility/nullspace explicitly, and retain symmetry as
`K_A = B_A M^-1 B_A^T`.  Treating an MPI tile edge as a physical boundary
would violate either policy; face geometry and multiplier values must be
exchanged before forming rank-interface corrections.

## Terrain metric

The diagonal energy is deliberately not claimed to be the full physical
kinetic-energy norm over steep terrain.  Terrain cross-coupling may be added
only through an explicitly positive face/cell metric, for example a
factorization `M = T^T W T` based on a documented staggered-to-physical
velocity transform.  Its inverse action and adjoint interpolation must be
defined together.  Adding slope terms to a gradient independently of `B`
would recreate the non-normal operator this design removes.

The diagonal operator is the first production candidate because it is
conservative, SPD, memory-bounded, and has an unambiguous acceptance test.
Physical comparison against the old correction determines whether a fuller
terrain metric is necessary; solver convergence alone does not.

## Required validation sequence

1. Serial algebra: adjoint identity, symmetry, positive energy, and
   manufactured exact projection.
2. Four-rank algebra: the same identities across internal MPI faces.
3. GPU equivalence: host and device applications agree to the selected
   precision tolerance.
4. Domain equivalence: `B q / V_c` matches `calc_divergence` exactly on
   active cells, including density and map factors.
5. Solver gate: true residual at most `1e-5`, status zero, and correction
   residual consistent with the independent `B(q_0 + delta_q)` check.
6. Scale gates: 250 m regional case, a larger bridge domain, then
   Switzerland at 200 m with the validated 80-level SLEVE 2/6 grid.
7. Physical gate: short finite output with bounded winds, plausible mass
   fluxes, and no degradation hidden by terrain or residual smoothing.

The solver may move from FGMRES to CG only after the distributed and GPU
symmetry/positive-energy gates pass.  A horizontal multilevel hierarchy with
vertical line relaxation remains appropriate for the strong vertical
anisotropy, but its transfer and coarse operators must preserve the same
energy.

HICAR stores physical winds in single precision while the Krylov multiplier
is double precision.  `wind_iterations` is therefore used as bounded
iterative refinement on the adjoint path: each pass recomputes `Bq` from the
actual stored winds, resets the multiplier correction, and resolves.  This
does not weaken the matrix residual gate; it closes the independently
recomputed physical constraint to the same standard despite storage
roundoff.

Enable the experimental projection with
`HICAR_WIND_ADJOINT_PROJECTION=1`.  This automatically requires and enables
the exact Galerkin multilevel preconditioner; running the adjoint projection
with only the local line smoother is not a supported production mode.
`HICAR_WIND_MULTILEVEL=1` remains available for testing that hierarchy with
the legacy projection operator.

Every production adjoint update recomputes the volume-weighted constraint
after the final stored-wind halo exchange.  The application aborts before
subsequent physics when its norm has not fallen below `2e-5` of the initial
constraint norm.
