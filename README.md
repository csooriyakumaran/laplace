# laplace

**3D compressible potential-flow solver for wind tunnel contractions.**

`laplace` computes the steady, inviscid, subsonic flow through a rectangular-section wind tunnel contraction. Its output fields are velocity, pressure, temperature and density, which a downstream Lagrangian droplet-trajectory solver consumes.

The method is deliberately simple:

- an algebraic grid built directly from the wall contours;
- a node-centred finite-volume discretization in computational coordinates;
- a line-relaxation solver.

Despite the name, the solver handles compressible flow. Laplace's equation is the incompressible limit ($\rho \equiv 1$) and the first validation case.

---

## 1. Governing equations

### Assumptions

The flow is steady, inviscid, adiabatic and shock-free. It issues from a settling chamber with uniform stagnation pressure $p_0$ and stagnation temperature $T_0$. The governing relations are then:

$$
\begin{aligned}
&\text{continuity:} && \nabla\cdot(\rho\,\mathbf{u}) = 0 \\
&\text{energy:} && h + \tfrac12 q^2 = h_0 \quad (\text{uniform}) \\
&\text{isentropic:} && p/\rho^{\gamma} = p_0/\rho_0^{\gamma} \quad (\text{uniform})
\end{aligned}
$$

where $\mathbf{u} = (u, v, w)$ and $q^2 = u^2 + v^2 + w^2$.

Uniform $h_0$ and uniform entropy make the flow irrotational by Crocco's theorem, $\mathbf{u}\times(\nabla\times\mathbf{u}) = \nabla h_0 - T\nabla s = 0$. The velocity is therefore the gradient of a potential:

$$
\mathbf{u} = \nabla\phi, \qquad u = \phi_x,\quad v = \phi_y,\quad w = \phi_z
$$

Under these assumptions this is exact, not an approximation.

### Full potential equation (conservative form)

Substituting into continuity:

$$
\frac{\partial}{\partial x}\left(\rho\,\phi_x\right) + \frac{\partial}{\partial y}\left(\rho\,\phi_y\right) + \frac{\partial}{\partial z}\left(\rho\,\phi_z\right) = 0
$$

### Density closure

Combining the energy and isentropic relations, with $h = a^2/(\gamma - 1)$, gives the compressible Bernoulli equation. In nondimensional form (§2):

$$
a^2 = 1 - \frac{\gamma - 1}{2}\, q^2, \qquad
\rho = \left(a^2\right)^{\frac{1}{\gamma - 1}} = \left[1 - \frac{\gamma - 1}{2}\left(\phi_x^2 + \phi_y^2 + \phi_z^2\right)\right]^{\frac{1}{\gamma - 1}}
$$

and

$$
T = a^2 = \rho^{\gamma - 1}, \qquad p = \rho^{\gamma}, \qquad M^2 = \frac{q^2}{a^2}
$$

### Quasi-linear (expanded) form

Differentiating the closure gives $\nabla\rho = -\dfrac{\rho}{a^2}\,\nabla\!\left(\tfrac12 q^2\right)$. Expanding the conservative form with this, and dividing by $\rho$:

$$
(a^2 - u^2)\,\phi_{xx} + (a^2 - v^2)\,\phi_{yy} + (a^2 - w^2)\,\phi_{zz}
- 2uv\,\phi_{xy} - 2vw\,\phi_{yz} - 2uw\,\phi_{xz} = 0
$$

This form shows the character of the equation. Aligning $x$ with the local velocity reduces it to:

$$
(1 - M^2)\,\phi_{ss} + \phi_{nn} + \phi_{mm} = 0
$$

So the equation is **elliptic for $M < 1$**. As $M \to 0$ ($a^2 \gg q^2$) it reduces to Laplace's equation, $\phi_{xx} + \phi_{yy} + \phi_{zz} = 0$.

The solver discretizes the **conservative** form, which conserves mass exactly. The quasi-linear form is for understanding only.

The sonic limit, $M = 1 \iff q^2 = 2/(\gamma + 1)$, is checked every iteration. Exceeding it means locally supersonic flow, which this scheme does not support.

---

## 2. Nondimensionalization

| Quantity | Reference |
|---|---|
| length | test-section half-height $h_\text{ts}$ |
| velocity | stagnation speed of sound $a_0$ |
| $\rho,\ p,\ T$ | $\rho_0,\ p_0,\ T_0$ |
| $\phi$ | $a_0\, h_\text{ts}$ |

The mass flux follows from the target test-section Mach number:

$$
m'' \equiv \frac{\rho u}{\rho_0 a_0} = M_\text{ts}\left[1 + \frac{\gamma - 1}{2} M_\text{ts}^2\right]^{-\frac{\gamma + 1}{2(\gamma - 1)}}
$$

The inlet flux is $m''$ scaled by the area ratio, $m''_\text{in} = m''\,A_\text{ts}/A_\text{in}$.

Stagnation conditions are invariant in inviscid flow, so the required inputs are $p_0$, $T_0$, $M_\text{ts}$, $\gamma$ and the wall contours.

---

## 3. Geometry and grid

### Domain

The contraction has a rectangular cross-section described by two wall contours:

- $h(x)$: half-height (walls at $y = \pm h$);
- $b(x)$: half-width (walls at $z = \pm b$).

A square inlet has $h = b$; a rectangular outlet has $h \neq b$. Both are extended upstream into the settling chamber by about one inlet half-height, and downstream into the test section by about two test-section half-heights. This keeps the uniform inlet and outlet conditions valid.

The section is symmetric about $y = 0$ and $z = 0$, so only the **quarter domain** $y \ge 0$, $z \ge 0$ is solved.

### Algebraic grid

The grid is a single structured block with $i$ streamwise ($\xi$), $j$ vertical ($\eta$) and $k$ spanwise ($\zeta$):

$$
x_{ijk} = X_i, \qquad y_{ijk} = h(X_i)\, s_j, \qquad z_{ijk} = b(X_i)\, t_k
$$

- **Streamwise stations $X_i$:** monotone. Use uniform spacing through the contraction and geometric stretching into the extensions.
- **Cross-stream distributions $s_j$ and $t_k$:** run from 0 on the symmetry plane to 1 at the wall, clustered toward the wall with a tanh stretching:

$$
s_j = 1 - \frac{\tanh\big(\beta\,(1 - \eta_j)\big)}{\tanh\beta}, \qquad \eta_j = \frac{j}{n_j - 1}
$$

and likewise for $t_k$ (typically $\beta \approx 1.5$–$2.5$).

Only $h$ and $b$ at the stations $X_i$ are needed. No wall derivatives are used.

Each constant-$i$ plane is flat, $y$ depends only on $(i, j)$, and $z$ depends only on $(i, k)$. The grid is non-orthogonal wherever a wall slopes.

```
Streamwise profile (x-y plane, schematic, not to scale)

  y
  ^
  |                                                    wall:  y = h(x)   (j = n_j-1)
  |   ________________
  |                   \________
  |                           \________
  |                                   \________________
  +------------------------------------------------------> x
  symmetry:  y = 0  (j = 0)
  ^                                                  ^
  i = 0 (inlet)                            i = n_i-1 (outlet)


Cross-section at one streamwise station (y-z plane); only y >= 0, z >= 0 is solved

        z
        ^
        |
  b(x)  +-----------------------+     wall:  z = b(x)   (k = n_k-1)
        |                       |
        |     solved quarter    |
        |        domain         |
        |                       |
      0 +-----------------------+-----> y
        0                     h(x)
        ^
        symmetry: y = 0 (j = 0) is the left edge; z = 0 (k = 0) is the bottom edge
```

The $b(x)$ taper in the cross-section is typically gentler than $h(x)$'s (square inlet, rectangular outlet, per the domain description above) — the two wall contours are independent, only sharing the same streamwise stations $X_i$.

---

## 4. Coordinate transformation

Computational coordinates $(\xi, \eta, \zeta)$ have unit spacing, so $\xi = i$, $\eta = j$, $\zeta = k$. This is the entire point of the transformation: the physical grid is curved and non-orthogonal, but every cell in computational space is an identical unit square, so the same stencil applies everywhere, and all boundary conditions become plain array-index statements instead of geometric ones.

```
  PHYSICAL SPACE  (x, y)                       COMPUTATIONAL SPACE  (xi, eta)

        (i,j+1)________(i+1,j+1)                     eta
              \          \                             ^
               \          \                   (i,j+1)  o-------o  (i+1,j+1)
                \__________\                           |       |
             (i,j)       (i+1,j)                       |       |
                                                 (i,j)  o-------o--> xi
          skewed, non-orthogonal cell                       (i+1,j)
          size and shape vary node to node          always a unit square, identical
                                                      everywhere, however stretched
                                                      the physical grid is

                      mapping:   xi = i,  eta = j   (unit spacing, by construction)
```

### General form

With Jacobian $J = \det\,\partial(x,y,z)/\partial(\xi,\eta,\zeta)$ and contravariant metric tensor $g^{mn} = \nabla\xi^m \cdot \nabla\xi^n$, the conservative equation transforms to:

$$
\frac{\partial}{\partial \xi^m}\left(\rho\, J\, g^{mn}\, \frac{\partial\phi}{\partial \xi^n}\right) = 0
$$

Written out:

$$
\frac{\partial}{\partial\xi}\Big(A^{\xi\xi}\phi_\xi + A^{\xi\eta}\phi_\eta + A^{\xi\zeta}\phi_\zeta\Big)
+ \frac{\partial}{\partial\eta}\Big(A^{\xi\eta}\phi_\xi + A^{\eta\eta}\phi_\eta + A^{\eta\zeta}\phi_\zeta\Big)
+ \frac{\partial}{\partial\zeta}\Big(A^{\xi\zeta}\phi_\xi + A^{\eta\zeta}\phi_\eta + A^{\zeta\zeta}\phi_\zeta\Big) = 0
$$

with $A^{mn} = \rho\, J\, g^{mn}$. The six independent coefficients form a symmetric $3\times 3$ tensor.

Each bracketed term is the mass flux through a face of constant $\xi$, $\eta$ or $\zeta$ respectively.

### Metrics for this grid

The grid's structure makes the Jacobian matrix lower-triangular:

$$
\frac{\partial(x,y,z)}{\partial(\xi,\eta,\zeta)} =
\begin{bmatrix} x_\xi & 0 & 0 \\ y_\xi & y_\eta & 0 \\ z_\xi & 0 & z_\zeta \end{bmatrix},
\qquad J = x_\xi\, y_\eta\, z_\zeta
$$

The inverse metrics are:

$$
\nabla\xi = \left(\frac{1}{x_\xi},\ 0,\ 0\right), \quad
\nabla\eta = \left(-\frac{y_\xi}{x_\xi\, y_\eta},\ \frac{1}{y_\eta},\ 0\right), \quad
\nabla\zeta = \left(-\frac{z_\xi}{x_\xi\, z_\zeta},\ 0,\ \frac{1}{z_\zeta}\right)
$$

so the coefficients are:

$$
\begin{aligned}
A^{\xi\xi} &= \rho\,\frac{y_\eta\, z_\zeta}{x_\xi}, &
A^{\xi\eta} &= -\rho\,\frac{y_\xi\, z_\zeta}{x_\xi}, &
A^{\xi\zeta} &= -\rho\,\frac{z_\xi\, y_\eta}{x_\xi}, \\[4pt]
A^{\eta\eta} &= \rho\,\frac{x_\xi\, z_\zeta}{y_\eta}\left(1 + \frac{y_\xi^2}{x_\xi^2}\right), &
A^{\zeta\zeta} &= \rho\,\frac{x_\xi\, y_\eta}{z_\zeta}\left(1 + \frac{z_\xi^2}{x_\xi^2}\right), &
A^{\eta\zeta} &= \rho\,\frac{y_\xi\, z_\xi}{x_\xi}
\end{aligned}
$$

The off-diagonal terms are driven by the wall slopes $y_\xi$ and $z_\xi$. They vanish in the constant-area extensions. $A^{\eta\zeta}$ is a product of two slopes, so it is small but nonzero.

### Physical velocity

$$
u = \frac{\phi_\xi}{x_\xi} - \frac{y_\xi}{x_\xi\, y_\eta}\,\phi_\eta - \frac{z_\xi}{x_\xi\, z_\zeta}\,\phi_\zeta, \qquad
v = \frac{\phi_\eta}{y_\eta}, \qquad
w = \frac{\phi_\zeta}{z_\zeta}
$$

### Freestream preservation

For uniform flow $\phi = x$, the derivatives are $\phi_\xi = x_\xi$, $\phi_\eta = \phi_\zeta = 0$. The fluxes then reduce to:

- $\rho\, y_\eta z_\zeta$ through $\xi$-faces;
- $-\rho\, y_\xi z_\zeta$ through $\eta$-faces;
- $-\rho\, z_\xi y_\eta$ through $\zeta$-faces.

These are exactly $\rho\,\hat{\mathbf{x}}\cdot\mathbf{S}$ for each face area vector. They telescope to zero over each cell *provided* the metrics are computed from the node coordinates with the same difference stencils as $\phi$.

> **Metric rule:** never use analytic metrics. Compute $x_\xi, y_\xi, y_\eta, z_\xi, z_\zeta$ from node coordinates with the stencils used for $\phi$.

---

## 5. Discretization

The scheme is node-centred finite volume. Node $(i,j,k)$ owns the control volume bounded by the faces $i\pm\tfrac12$, $j\pm\tfrac12$, $k\pm\tfrac12$:

![Node-centred control volume, showing the node, its six bounding faces with their fluxes F/G/H, and the six neighbouring nodes](docs/fv-cell.png)

The face fluxes are:

$$
\begin{aligned}
F_{i+\frac12} &= \big(A^{\xi\xi}\phi_\xi + A^{\xi\eta}\phi_\eta + A^{\xi\zeta}\phi_\zeta\big)_{i+\frac12,j,k} \\
G_{j+\frac12} &= \big(A^{\xi\eta}\phi_\xi + A^{\eta\eta}\phi_\eta + A^{\eta\zeta}\phi_\zeta\big)_{i,j+\frac12,k} \\
H_{k+\frac12} &= \big(A^{\xi\zeta}\phi_\xi + A^{\eta\zeta}\phi_\eta + A^{\zeta\zeta}\phi_\zeta\big)_{i,j,k+\frac12}
\end{aligned}
$$

and the residual is:

$$
\mathcal{R}_{i,j,k} = F_{i+\frac12} - F_{i-\frac12} + G_{j+\frac12} - G_{j-\frac12} + H_{k+\frac12} - H_{k-\frac12} = 0
$$

### Face derivatives

On a $\xi$-face $(i+\tfrac12, j, k)$:

$$
\begin{aligned}
(\cdot)_\xi &= (\cdot)_{i+1,j,k} - (\cdot)_{i,j,k} \\
(\cdot)_\eta &= \tfrac14\big[(\cdot)_{i+1,j+1,k} + (\cdot)_{i,j+1,k} - (\cdot)_{i+1,j-1,k} - (\cdot)_{i,j-1,k}\big] \\
(\cdot)_\zeta &= \tfrac14\big[(\cdot)_{i+1,j,k+1} + (\cdot)_{i,j,k+1} - (\cdot)_{i+1,j,k-1} - (\cdot)_{i,j,k-1}\big]
\end{aligned}
$$

$\eta$- and $\zeta$-faces use the same stencils with indices permuted. Where a 4-point stencil would leave the grid, it falls back to a one-sided 2-point average.

The resulting node stencil has **19 points**: the 7-point core plus 4 neighbours for each of the three cross-derivative pairs.

### Face density

$q^2$ is evaluated at each face from $(u, v, w)$ above. The closure then gives $\rho$, which multiplies the geometric part of $A^{mn}$.

---

## 6. Boundary conditions

Every boundary is a grid plane. Nodes on a boundary face, edge or corner own $\tfrac12$, $\tfrac14$ or $\tfrac18$ control volumes. Fluxes through faces cut by a boundary are scaled by the same fraction.

| Boundary | Location | Treatment |
|---|---|---|
| Symmetry | $j = 0$, $k = 0$ | zero face flux, $G_{-\frac12} = 0$, $H_{-\frac12} = 0$ |
| Walls | $j = n_j - 1$, $k = n_k - 1$ | zero face flux (no penetration) |
| Inlet | $i = 0$ | prescribed flux $F_{-\frac12} = m''_\text{in}\; y_\eta\, z_\zeta$ (uniform mass flux × face area) |
| Outlet | $i = n_i - 1$ | Dirichlet $\phi = 0$ (uniform flow normal to the outlet plane; also fixes the datum) |

---

## 7. Solution: SLOR with density lagging

The cross terms make the operator non-symmetric. The solver uses successive line over-relaxation along $j$-lines.

Each node equation is written as:

$$
a_W\phi_{i-1} + a_E\phi_{i+1} + a_S\phi_{j-1} + a_N\phi_{j+1} + a_B\phi_{k-1} + a_T\phi_{k+1} - a_P\phi_{i,j,k} = -b_{i,j,k}
$$

where unlisted indices are unchanged, and the coefficients are:

$$
a_{W,E} = A^{\xi\xi}_{i\mp\frac12}, \quad
a_{S,N} = A^{\eta\eta}_{j\mp\frac12}, \quad
a_{B,T} = A^{\zeta\zeta}_{k\mp\frac12}, \quad
a_P = \textstyle\sum a_\text{nb}
$$

$b_{i,j,k}$ collects all cross-derivative terms ($A^{\xi\eta}$, $A^{\xi\zeta}$, $A^{\eta\zeta}$), evaluated from the latest $\phi$.

The algorithm:

1. Set $\rho \equiv 1$ and initialize $\phi$ from the quasi-1D solution, or simply $\phi = x$.
2. For each $i$ from inlet to outlet, and each $k$:
   - Solve the tridiagonal system in $j$ ($a_S,\ -a_P,\ a_N$) with the Thomas algorithm. Already-updated neighbours are used where available.
   - Over-relax: $\phi \leftarrow \phi + \omega(\phi^* - \phi)$, with $\omega \approx 1.5$–$1.8$.
3. Every few sweeps, update the face densities with under-relaxation, $\rho \leftarrow \rho + \omega_\rho(\rho_\text{new} - \rho)$, and check the sonic limit.
4. Repeat until $\max|\mathcal{R}|$ and $\max|\Delta\rho|$ fall below tolerance.

Storage is structure-of-arrays, `[i][k][j]` with $j$ contiguous, so each line solve is a unit-stride sweep.

If convergence becomes a bottleneck on large grids, the upgrade path is:

- BiCGSTAB with ILU(0) on the 19-diagonal operator;
- approximate factorization (AF2);
- multigrid.

---

## 8. Output and particle lookup

Nodal velocities come from central differences of $\phi$ (one-sided at boundaries) through §4. Then $\rho$, $p$, $T$ and $M$ follow from §1, and dimensional values from §2. The full duct is recovered from the quarter domain by reflection, with $v \to -v$ across $y = 0$ and $w \to -w$ across $z = 0$.

Particle lookup needs no Newton inversion:

1. Find $i$ by bisection on $X_i$.
2. Compute $s = |y|/h(x)$ and $t = |z|/b(x)$, which give $j$ and $k$ directly.
3. Interpolate trilinearly.

---

## 9. Validation path

1. **Incompressible straight duct** ($\rho \equiv 1$, constant $h$, $b$): uniform flow to machine precision.
2. **Incompressible contraction:** Laplace's equation. Check mass conservation at every station and grid convergence.
3. **Symmetric contraction** ($h = b$): solution symmetric under $y \leftrightarrow z$.
4. **Compressible straight duct:** exact isentropic $\rho$, $p$, $T$ for the prescribed $m''$.
5. **Compressible contraction:** centreline against quasi-1D isentropic theory; stagnation quantities constant everywhere.

---

## Scope

**In scope:** steady, inviscid, shock-free subsonic flow in rectangular-section ducts with independent height and width contours.

**Planned:** wall boundary layer and separation screening by an integral method driven by the wall pressure distributions along both wall pairs. Corners, where the two adverse pressure gradients meet, are the critical region.

**Out of scope:** rotational inflow (screens, non-uniform settling-chamber profiles), transonic flow, and filleted or non-rectangular sections.
