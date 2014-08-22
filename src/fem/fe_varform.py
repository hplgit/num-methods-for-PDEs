# purely numerical version

from fe_approx1D_numint import mesh
from fe_approx1D import basis

def finite_element1D(
    vertices, cells, dof_map,
    integrand_lhs, boundary_lhs,
    integrand_rhs, boundary_rhs,
    intrule='GaussLegendre'):
    import sets
    N_n = len(list(set(np.array(dof_map).ravel())))
    N_e = len(cells)
    A = np.zeros((N_n, N_n))
    b = np.zeros(N_n)

    for e in range(N_e):
        Omega_e = [vertices[cells[e][0]], vertices[cells[e][1]]]
        n = len(dof_map[e])
        h = Omega_e[1] - Omega_e[0]
        A_e = np.zeros((n, n))
        b_e = np.zeros(n)

        <Declare element matrix and vector: A_e, b_e>

        # Integrate over the reference cell
        points, weights = <numerical integration rule>
        for X, w in zip(points, weights):
            phi = <basis functions and derivatives at X>
            detJ = h/2
            x = 0.5*(Omega_e[0] + Omega_e[1]) + X*h
            for r in range(n):
                for s in range(n):
                    A_e[r,s] += integrand_lhs(phi, r, s, x)*detJ*w
                b_e[r] += integrand_rhs(phi, r, x)*detJ*w

        # Add boundary terms
        for r in range(n):
            for s in range(n):
                A_e[r,s] += boundary_lhs(phi, r, s, x)*detJ*w
            b_e[r] += boundary_rhs(phi, r, x)*detJ*w

        # Incorporate essential boundary conditions
        for r in range(n):
            global_dof = dof_map[e][r]
            if global_dof in essbc_dofs:
                # dof r is subject to an essential condition
                value = essbc_docs[global_dof]
                # Symmetric modification
                b_e -= value*A_e[:,r]
                A_e[r,:] = 0
                A_e[:,r] = 0
                A_e[r,r] = 1
                b_e[r] = value

        # Assemble
        for r in range(n):
            for s in range(n):
                A[dof_map[e][r], dof_map[e][r]] += A_e[r,s]
            b[dof_map[e][r] += b_e[r]

    return A, b

    h = Omega_e[1] - Omega_e[0]
    detJ = h/2  # dx/dX
    points, weights = integration_rule(d, intrule)
    for X, w in zip(points, weights):
        phi = {0: [phi_r(r, X, d) for r in range(n)]}
        # must perform symbolic diff in a local phi routine here...
        for r in range(n):
            for s in range(r, n):
                A_e[r,s] += integrand_lhs[phi, r, s, X, x)*detJ*w[j]
                A_e[s,r] = A_e[r,s]
    return A_e

def integration_rule(d, tp='GaussLegendre'):
    import numint
    if tp == 'GaussLegendre':
        return numint.GaussLegendre(d+1)
    elif tp == 'NewtonCotes':
        return numint.NewtonCotes(d+1)


def basis(d, X):
    # return phi dict with basis funcs and derivatives
    # create them with symbolic X, diff, then subs X by the value here
    pass


A, b = finite_element1D(...)
c = np.linalg.solve(A, b)
