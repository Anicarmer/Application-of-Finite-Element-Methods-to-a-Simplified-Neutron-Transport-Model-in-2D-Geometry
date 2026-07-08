import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as tri

Q         = 1.0
MU        = 1.0 / np.sqrt(2)
NU        = 1.0 / np.sqrt(2)
SIGMA_REF = 1.0


def source_term(x, y):
    """Constant source term f = q."""
    return Q * np.ones_like(x, dtype=float)


def exact_solution(x, y, sigma, mu=MU, nu=NU):
    if np.isclose(mu, 0.0):
        tau = y / nu
    elif np.isclose(nu, 0.0):
        tau = x / mu
    else:
        tau = np.minimum(x / mu, y / nu)

    if np.isclose(sigma, 0.0):
        return Q * tau
    return (Q / sigma) * (1.0 - np.exp(-sigma * tau))


def exact_gradient(xq, yq, sigma, mu=MU, nu=NU):
    if np.isclose(mu, 0.0):
        tau = yq / nu
    elif np.isclose(nu, 0.0):
        tau = xq / mu
    else:
        tau = np.minimum(xq / mu, yq / nu)

    if np.isclose(sigma, 0.0):
        dudt = Q * np.ones_like(tau)
    else:
        dudt = Q * np.exp(-sigma * tau)

    if np.isclose(mu, 0.0):
        grad_x = np.zeros_like(dudt)
        grad_y = dudt / nu
    elif np.isclose(nu, 0.0):
        grad_x = dudt / mu
        grad_y = np.zeros_like(dudt)
    else:
        grad_x = np.where(xq / mu <= yq / nu, dudt / mu, 0.0)
        grad_y = np.where(yq / nu <  xq / mu, dudt / nu, 0.0)

    return grad_x, grad_y


def solve_transport(N, sigma, mu=MU, nu=NU):
    x = np.linspace(0, 1, N + 1)
    y = np.linspace(0, 1, N + 1)
    X, Y = np.meshgrid(x, y)
    nodes = np.column_stack([X.ravel(), Y.ravel()])

    elements = []
    for j in range(N):
        for i in range(N):
            n1 = j * (N + 1) + i
            n2 = n1 + 1
            n3 = n1 + (N + 1)
            n4 = n3 + 1
            elements.append([n1, n2, n4])
            elements.append([n1, n4, n3])
    elements = np.array(elements)

    n_nodes = len(nodes)
    A_mat = np.zeros((n_nodes, n_nodes))
    F_vec = np.zeros(n_nodes)

    quad_points = [(1/6, 1/6, 2/3),
                   (1/6, 2/3, 1/6),
                   (2/3, 1/6, 1/6)]

    for elem in elements:
        coords = nodes[elem]
        x1, y1 = coords[0]; x2, y2 = coords[1]; x3, y3 = coords[2]
        area = 0.5 * abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))

        B    = np.array([[1, x1, y1],
                         [1, x2, y2],
                         [1, x3, y3]])
        invB = np.linalg.inv(B)
        dphidx = invB[1, :]
        dphidy = invB[2, :]

        d = mu * dphidx + nu * dphidy

        Mloc = sigma * (area / 12.0) * np.array([[2, 1, 1],
                                                   [1, 2, 1],
                                                   [1, 1, 2]])
        Cloc = (area / 3.0) * np.outer(np.ones(3), d)

        Aloc = Mloc + Cloc

        Floc = np.zeros(3)
        for l1, l2, l3 in quad_points:
            xq = l1*x1 + l2*x2 + l3*x3
            yq = l1*y1 + l2*y2 + l3*y3
            fq = source_term(xq, yq)
            Floc += area * (1.0/3.0) * fq * np.array([l1, l2, l3])

        for i in range(3):
            I = elem[i]
            F_vec[I] += Floc[i]
            for j in range(3):
                A_mat[I, elem[j]] += Aloc[i, j]

    for k, (xk, yk) in enumerate(nodes):
        is_inflow = False
        if not np.isclose(mu, 0.0) and np.isclose(xk, 0.0):
            is_inflow = True
        if not np.isclose(nu, 0.0) and np.isclose(yk, 0.0):
            is_inflow = True
        if np.isclose(mu, 0.0) and np.isclose(yk, 0.0):
            is_inflow = True
        if np.isclose(nu, 0.0) and np.isclose(xk, 0.0):
            is_inflow = True
        if is_inflow:
            A_mat[k, :] = 0.0; A_mat[:, k] = 0.0
            A_mat[k, k] = 1.0; F_vec[k]    = 0.0

    U = np.linalg.solve(A_mat, F_vec)

    U_exact   = exact_solution(nodes[:, 0], nodes[:, 1], sigma, mu, nu)
    max_error = np.max(np.abs(U_exact - U))

    l2_error_sq = 0.0
    for elem in elements:
        coords = nodes[elem]
        x1, y1 = coords[0]; x2, y2 = coords[1]; x3, y3 = coords[2]
        area = 0.5 * abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))
        Uloc = U[elem]
        for l1, l2, l3 in quad_points:
            xq = l1*x1 + l2*x2 + l3*x3
            yq = l1*y1 + l2*y2 + l3*y3
            u_exact_q = exact_solution(xq, yq, sigma, mu, nu)
            uh = l1*Uloc[0] + l2*Uloc[1] + l3*Uloc[2]
            l2_error_sq += area * (1.0/3.0) * (u_exact_q - uh)**2
    l2_error = np.sqrt(l2_error_sq)

    h1_error_sq = 0.0
    for elem in elements:
        coords = nodes[elem]
        x1, y1 = coords[0]; x2, y2 = coords[1]; x3, y3 = coords[2]
        area = 0.5 * abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))
        B    = np.array([[1, x1, y1],
                         [1, x2, y2],
                         [1, x3, y3]])
        invB = np.linalg.inv(B)
        dphidx = invB[1, :]
        dphidy = invB[2, :]
        Uloc = U[elem]
        grad_uh_x = dphidx @ Uloc
        grad_uh_y = dphidy @ Uloc
        for l1, l2, l3 in quad_points:
            xq = l1*x1 + l2*x2 + l3*x3
            yq = l1*y1 + l2*y2 + l3*y3
            grad_u_x, grad_u_y = exact_gradient(xq, yq, sigma, mu, nu)
            diff_x = grad_u_x - grad_uh_x
            diff_y = grad_u_y - grad_uh_y
            h1_error_sq += area * (1.0/3.0) * (diff_x**2 + diff_y**2)
    h1_error = np.sqrt(h1_error_sq)

    return nodes, elements, U, max_error, l2_error, h1_error


N_REF = 8

nodes, elements, U, max_error, l2_error, h1_error = solve_transport(
    N_REF, SIGMA_REF, MU, NU)

print("=" * 60)
print(f"Single mesh run:  N = {N_REF},  sigma = {SIGMA_REF},  "
      f"b = (1/sqrt(2), 1/sqrt(2))")
print("=" * 60)
print("Number of nodes    =", len(nodes))
print("Number of elements =", len(elements))
print(f"Max error = {max_error:.6e}")
print(f"L2  error = {l2_error:.6e}")
print(f"H1  error = {h1_error:.6e}")

triang = tri.Triangulation(nodes[:, 0], nodes[:, 1], elements)

fig = plt.figure(figsize=(8, 6))
ax  = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(triang, U, cmap='viridis')
ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('$u_h$')
ax.set_title(
    f'FEM Approximation  ($f={Q}$, $\\sigma={SIGMA_REF}$, $N={N_REF}$,  '
    r'$\mathbf{b}=\frac{1}{\sqrt{2}}(1,1)$)')
plt.tight_layout()
plt.show()


mesh_sizes = [4, 8, 16, 32]
h_values   = [1.0 / N for N in mesh_sizes]

errors    = []
l2_errors = []
h1_errors = []

for Nm in mesh_sizes:
    _, _, _, me, l2e, h1e = solve_transport(Nm, SIGMA_REF, MU, NU)
    errors.append(me)
    l2_errors.append(l2e)
    h1_errors.append(h1e)

print(f"\nConvergence study:  sigma = {SIGMA_REF},  "
      f"b = (1/sqrt(2), 1/sqrt(2))")
print(f"{'N':>4}  {'h':>8}  {'Max error':>12}  "
      f"{'L2 error':>12}  {'H1 seminorm':>12}")
print("-" * 56)
for Nm, h, me, l2e, h1e in zip(
        mesh_sizes, h_values, errors, l2_errors, h1_errors):
    print(f"{Nm:>4}  {h:>8.5f}  {me:>12.6e}  "
          f"{l2e:>12.6e}  {h1e:>12.6e}")

print("\nConvergence rates (Max error):")
for i in range(1, len(errors)):
    p = np.log(errors[i-1] / errors[i]) / np.log(2)
    print(f"  N={mesh_sizes[i-1]:>2}→{mesh_sizes[i]:>2}: rate = {p:.4f}")

print("\nConvergence rates (L2 error):")
for i in range(1, len(l2_errors)):
    p = np.log(l2_errors[i-1] / l2_errors[i]) / np.log(2)
    print(f"  N={mesh_sizes[i-1]:>2}→{mesh_sizes[i]:>2}: rate = {p:.4f}")

print("\nConvergence rates (H1 seminorm):")
for i in range(1, len(h1_errors)):
    p = np.log(h1_errors[i-1] / h1_errors[i]) / np.log(2)
    print(f"  N={mesh_sizes[i-1]:>2}→{mesh_sizes[i]:>2}: rate = {p:.4f}")

h_arr = np.array(h_values)

plt.figure(figsize=(6, 5))
plt.loglog(h_arr, errors,    'o-', label='Max error')
plt.loglog(h_arr, l2_errors, 's-', label='$L^2$ error')
plt.loglog(h_arr, h1_errors, '^-', label='$H^1$ seminorm')
plt.loglog(h_arr, 2.0 * h_arr**2, 'k--', label='$O(h^2)$ ref')
plt.loglog(h_arr, 1.0 * h_arr**1, 'k:',  label='$O(h)$ ref')
plt.xlabel('$h$'); plt.ylabel('Error')
plt.title(f'Convergence Study  ($\\sigma = {SIGMA_REF}$,  '
          r'$\mathbf{b}=\frac{1}{\sqrt{2}}(1,1)$)')
plt.grid(True, which='both', ls='--', alpha=0.4)
plt.legend()
plt.tight_layout()
plt.show()

RESONANCE_FILE = "resonance_data_cleaned.csv"

res       = pd.read_csv(RESONANCE_FILE)
Er_arr    = res["Resonance Energy (eV)"].to_numpy()
Gtot_arr  = res["Total Width (eV)"].to_numpy()
Gn_arr    = res["Gamma-N (eV)"].to_numpy()
Ggam_arr  = res["Gamma-Gamma (eV)"].to_numpy()

n_res = len(Er_arr)

h_planck = 6.62607015e-34
m_n      = 1.67492749804e-27
eV_to_J  = 1.602176634e-19
BARN     = 1e-28
G_SPIN   = 0.5


def lambda_r_squared_barn(Er_eV):
    """Squared de Broglie wavelength at resonance energy |Er|, in barns."""
    Er_J = np.abs(Er_eV) * eV_to_J
    lam  = h_planck / np.sqrt(2.0 * m_n * Er_J)
    return (lam ** 2) / BARN


lambda2_arr = lambda_r_squared_barn(Er_arr)

print("\n" + "=" * 60)
print("Multi-resonance Breit-Wigner absorption cross-section")
print(f"  Number of resonances : {n_res}")
print(f"  Energy range         : {Er_arr.min():.2f} to {Er_arr.max():.2f} eV")
print(f"  Statistical factor g : {G_SPIN} (averaged, J=3/J=4 states)")
print("=" * 60)


def sigma_absorption(E):
    """
    Total radiative-capture (absorption) cross-section (barn), summed
    incoherently over all resolved resonances.

    E : scalar or ndarray, incident neutron energy (eV)
    """
    E = np.atleast_1d(np.asarray(E, dtype=float))
    total = np.zeros_like(E)
    for Er_i, Gtot_i, Gn_i, Gg_i, l2_i in zip(
            Er_arr, Gtot_arr, Gn_arr, Ggam_arr, lambda2_arr):
        denom = (E - Er_i) ** 2 + (Gtot_i / 2.0) ** 2
        total += (l2_i * G_SPIN / (4.0 * np.pi)) * Gn_i * Gg_i / denom
    return total if total.size > 1 else total[0]


E_table = np.array([0.05, 0.29, 1.14, 5.6, 8.78, 20.1, 40.5, 60.2, 80.4, 86.9, 100.0])

print(f"{'E (eV)':>10}  {'sigma_a (barn)':>16}")
print("-" * 30)
for E in E_table:
    print(f"{E:>10.2f}  {sigma_absorption(E):>16.6e}")

E_MIN, E_MAX = 0.05, 100.0

E_background = np.linspace(E_MIN, E_MAX, 4000)

refined_blocks = [E_background]
for Er_i, Gtot_i in zip(Er_arr, Gtot_arr):
    if Er_i <= 0:
        continue
    half_width = max(5.0 * Gtot_i, 0.05)
    lo = max(E_MIN, Er_i - half_width)
    hi = min(E_MAX, Er_i + half_width)
    refined_blocks.append(np.linspace(lo, hi, 60))

E_plot = np.unique(np.concatenate(refined_blocks))
E_plot = E_plot[(E_plot >= E_MIN) & (E_plot <= E_MAX)]

sigma_vals = sigma_absorption(E_plot)

fig, ax = plt.subplots(figsize=(9, 5))
ax.semilogy(E_plot, sigma_vals, 'b-', lw=0.8)
ax.scatter(E_table, sigma_absorption(E_table), color='k', zorder=5,
           s=15, label='Table points')
ax.set_xlabel('$E$ (eV)')
ax.set_ylabel('$\\sigma_a(E)$ (barn)')
ax.set_title(f'Full resonance spectrum ({n_res} resonances)\n'
             'Radiative-capture (absorption) cross-section')
ax.legend(fontsize=9)
ax.grid(True, which='both', ls='--', alpha=0.4)
plt.tight_layout()
plt.show()

N_BW = 32
N_ENERGY_FEM = 200

E_fem = np.geomspace(E_MIN, E_MAX, N_ENERGY_FEM)

max_errors_bw = []
l2_errors_bw  = []
h1_errors_bw  = []
max_uh_bw     = []

for E in E_fem:
    sig = sigma_absorption(E)
    _, _, U_bw, me, l2e, h1e = solve_transport(N_BW, sig, MU, NU)
    max_errors_bw.append(me)
    l2_errors_bw.append(l2e)
    h1_errors_bw.append(h1e)
    max_uh_bw.append(np.max(U_bw))

max_errors_bw = np.array(max_errors_bw)
l2_errors_bw  = np.array(l2_errors_bw)
h1_errors_bw  = np.array(h1_errors_bw)
max_uh_bw     = np.array(max_uh_bw)

plt.figure(figsize=(8, 5))
plt.loglog(E_fem, max_errors_bw, '-', label='Max error')
plt.loglog(E_fem, l2_errors_bw,  '-', label='$L^2$ error')
plt.loglog(E_fem, h1_errors_bw,  '-', label='$H^1$ seminorm')
for Er_i in Er_arr[Er_arr > 0]:
    plt.axvline(Er_i, color='r', ls='--', lw=0.3, alpha=0.25)
plt.xlabel('$E$ (eV)')
plt.ylabel('Error')
plt.title(f'Error norms across the full resonance spectrum\n'
          f'($N={N_BW}$, $q={Q}$,  '
          r'$\mathbf{b}=\frac{1}{\sqrt{2}}(1,1)$)')
plt.legend()
plt.grid(True, which='both', ls='--', alpha=0.4)
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 4.5))
plt.semilogx(E_fem, max_uh_bw, 'b-', label='$\\max(u_h)$')
plt.xlabel('$E$ (eV)')
plt.ylabel('$\\max(u_h)$')
plt.title(f'Maximum nodal value across the full resonance spectrum\n'
          f'($N={N_BW}$, $q={Q}$)')
plt.legend()
plt.grid(True, which='both', ls='--', alpha=0.4)
plt.tight_layout()
plt.show()

print(f"\nMulti-resonance Breit-Wigner FEM study:  N = {N_BW},  q = {Q},  "
      f"b = (1/sqrt(2), 1/sqrt(2))")
print(f"{'E (eV)':>8}  {'sigma_a (barn)':>16}  {'Max error':>12}  "
      f"{'L2 error':>12}  {'H1 seminorm':>12}  {'max(u_h)':>10}")
print("-" * 82)
for E in E_table:
    sig = sigma_absorption(E)
    _, _, U_bw, me, l2e, h1e = solve_transport(N_BW, sig, MU, NU)
    print(f"{E:>8.2f}  {sig:>16.6e}  {me:>12.6e}  "
          f"{l2e:>12.6e}  {h1e:>12.6e}  {np.max(U_bw):>10.6f}")


directions = [
    (1.0,          0.0         ),
    (0.0,          1.0         ),
    (1/np.sqrt(2), 1/np.sqrt(2)),
    (0.5,          0.5         ),
    (0.4,          0.3         ),
]
dir_labels = [
    r'$\mathbf{b} = (1,\,0)$',
    r'$\mathbf{b} = (0,\,1)$',
    r'$\mathbf{b} = \frac{1}{\sqrt{2}}(1,\,1)$',
    r'$\mathbf{b} = \frac{1}{2}(1,\,1)$',
    r'$\mathbf{b} = (0.4,\,0.3)$',
]
dir_print_labels = [
    '(1, 0)           |b|=1',
    '(0, 1)           |b|=1',
    '(1,1)/sqrt(2)    |b|=1',
    '(1,1)/2       |b|=0.71',
    '(0.4, 0.3)    |b|=0.50',
]

SIGMA_DIR  = 1.0
N_DIR      = 32
colors_dir = plt.cm.tab10(np.linspace(0, 0.5, len(directions)))

fig, axes = plt.subplots(2, 3, figsize=(14, 8),
                         subplot_kw={'projection': '3d'})
axes = axes.ravel()

print("\n" + "=" * 70)
print(f"Transport direction study:  sigma = {SIGMA_DIR},  N = {N_DIR}")
print("=" * 70)
print(f"{'Direction':>24}  {'Max error':>12}  "
      f"{'L2 error':>12}  {'H1 seminorm':>12}")
print("-" * 66)

max_errors_d = []
l2_errors_d  = []
h1_errors_d  = []

for idx, ((mu, nu), label, plabel) in enumerate(
        zip(directions, dir_labels, dir_print_labels)):
    nodes_d, elems_d, U_d, me, l2e, h1e = solve_transport(
        N_DIR, SIGMA_DIR, mu, nu)
    max_errors_d.append(me)
    l2_errors_d.append(l2e)
    h1_errors_d.append(h1e)
    print(f"{plabel:>24}  {me:>12.6e}  {l2e:>12.6e}  {h1e:>12.6e}")
    triang_d = tri.Triangulation(nodes_d[:, 0], nodes_d[:, 1], elems_d)
    ax = axes[idx]
    ax.plot_trisurf(triang_d, U_d, cmap='viridis', alpha=0.9)
    ax.set_xlabel('x', fontsize=8)
    ax.set_ylabel('y', fontsize=8)
    ax.set_zlabel('$u_h$', fontsize=8)
    ax.set_title(label, fontsize=10)

axes[-1].set_visible(False)
fig.suptitle(
    f'FEM solution for varying transport direction\n'
    f'($\\sigma={SIGMA_DIR}$, $N={N_DIR}$, $q={Q}$)',
    fontsize=12)
plt.tight_layout()
plt.show()

plt.figure(figsize=(7, 5))
x_line = np.linspace(0, 1, 300)
y_line = 0.5 * np.ones_like(x_line)

for (mu, nu), label, col in zip(directions, dir_labels, colors_dir):
    u_line = exact_solution(x_line, y_line, SIGMA_DIR, mu, nu)
    plt.plot(x_line, u_line, color=col, label=label)

plt.xlabel('$x$  (along $y = 0.5$)')
plt.ylabel('$u$')
plt.title('Exact solution profiles along $y = 0.5$\n'
          f'for varying transport direction  ($\\sigma={SIGMA_DIR}$, $q={Q}$)')
plt.legend(fontsize=9)
plt.grid(True, ls='--', alpha=0.4)
plt.tight_layout()
plt.show()

x_pos = np.arange(len(directions))
plt.figure(figsize=(7, 5))
plt.semilogy(x_pos, max_errors_d, 'o-', label='Max error')
plt.semilogy(x_pos, l2_errors_d,  's-', label='$L^2$ error')
plt.semilogy(x_pos, h1_errors_d,  '^-', label='$H^1$ seminorm')
plt.xticks(x_pos, [p.split()[0] for p in dir_print_labels], fontsize=8)
plt.xlabel('Transport direction $\\mathbf{b}$')
plt.ylabel('Error')
plt.title(f'Error norms vs transport direction\n'
          f'($\\sigma={SIGMA_DIR}$, $N={N_DIR}$, $q={Q}$)')
plt.legend()
plt.grid(True, which='both', ls='--', alpha=0.4)
plt.tight_layout()
plt.show()
