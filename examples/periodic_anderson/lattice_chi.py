#!/usr/bin/env python3
"""Approximate PAM lattice spin susceptibility from continued cluster data.

Input NPZ arrays (complex values may be stored as complex arrays or real/imag
last-axis pairs):
  omega: (nw,) real frequencies, in the same energy units as t, V, U
  chi_cluster: (nQ, nw) continued physical cluster magnetic susceptibility
  greens_cluster: (nK, nw, 2, 2) retarded cluster Green function in (c,f) basis
  sigma_cluster: (nK, nw, 2, 2) retarded cluster self-energy in (c,f) basis
  q_cluster: (nQ, 2) cluster transfer momenta in radians, in [-pi,pi)
  k_cluster: (nK, 2) cluster momenta in radians, matching sigma_cluster

Usage:
  python lattice_chi.py input.npz output.npz --beta 10 --t 1 --epsilon-f -2 \
      --V 0.5 --V-prime 0 --U 6

Approximation: infer a scalar transfer-frequency kernel from the physical
cluster susceptibility and a cluster bubble, then reuse it in the lattice BSE.
This is not the full fermionic-frequency irreducible vertex. The definition of
the physical susceptibility normalization and the spin factor must match the
continuation. The default convention is chi = chi0 / (1 - V_eff chi0), with
chi0 = -(1/Nk) sum_k integral dnu/pi Im Tr[G(k,nu) G(k+q,nu+w)] weighted by
the Fermi difference. Set --spin-factor to 1 or 2 as appropriate to the input
convention. Verify the static cluster reconstruction before interpreting data.

The Sigma(k,nu) used to build the lattice Green function is assigned to the
nearest cluster momentum (standard piecewise-constant DCA embedding).
"""

import argparse
import numpy as np


def complex_array(x):
    x = np.asarray(x)
    if not np.iscomplexobj(x) and x.ndim and x.shape[-1] == 2:
        return x[..., 0] + 1j * x[..., 1]
    return x.astype(complex)


def wrap_k(k):
    return (np.asarray(k) + np.pi) % (2 * np.pi) - np.pi


def assign_sigma(kpoints, kcluster, sigma):
    delta = wrap_k(kpoints[:, None, :] - kcluster[None, :, :])
    return sigma[np.argmin(np.sum(delta * delta, axis=-1), axis=1)]


def fermi(x, beta):
    y = np.clip(beta * np.asarray(x), -700, 700)
    return 1. / (np.exp(y) + 1.)


def h0_pam(k, t, eps_f, v, vp):
    c = np.cos(k[:, 0]) + np.cos(k[:, 1])
    hk = np.zeros((len(k), 2, 2), dtype=complex)
    hk[:, 0, 0] = -2 * t * c
    hk[:, 1, 1] = eps_f
    hk[:, 0, 1] = hk[:, 1, 0] = v + 2 * vp * c
    return hk


def greens(omega, k, sigma_k, args):
    # G(k,w) = [(w + mu + i0+)I - H0(k) - Sigma_K(w)]^-1
    h = h0_pam(k, args.t, args.epsilon_f, args.V, args.V_prime)
    z = omega[None, :, None, None] + args.mu + 1j * args.eta
    eye = np.eye(2)[None, None, :, :]
    inv = z * eye - h[:, None, :, :] - sigma_k
    return np.linalg.inv(inv)


def bubble_at_q(omega, q, k, g, beta, spin_factor, eta):
    """Retarded bubble from spectral representation, orbital matrix trace.

    chi0(q,w) = -s/Nk sum_k int dnu dnu' A_k(nu) A_k+q(nu')
                  [f(nu)-f(nu')]/[w+nu-nu'+i eta].
    The frequency integrals use trapezoid weights. Runtime scales as
    O(Nk * Nw^2 * Nout_w); use moderate omega and k meshes for initial runs.
    """
    qk = wrap_k(k + q)
    gq = g(qk)
    a = -np.imag(g)
    aq = -np.imag(gq)
    f = fermi(omega, beta)
    dnu = np.gradient(omega)
    out = np.empty(len(omega), dtype=complex)
    for iw, w in enumerate(omega):
        denom = w + omega[:, None] - omega[None, :] + 1j * eta
        thermal = f[:, None] - f[None, :]
        trace = np.einsum("knab,kmba->knm", a, aq, optimize=True)
        out[iw] = -spin_factor * np.einsum(
            "n,m,knm,nm->", dnu, dnu, trace, thermal / denom, optimize=True
        ) / (np.pi**2 * len(k))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", help="NPZ input containing omega, chi_cluster, sigma_cluster, q_cluster, k_cluster")
    p.add_argument("output", help="NPZ output for q, omega, chi_lattice, chi0_lattice, V_eff, chi0_cluster")
    p.add_argument("--beta", type=float, required=True)
    p.add_argument("--t", type=float, required=True)
    p.add_argument("--epsilon-f", type=float, required=True)
    p.add_argument("--V", type=float, required=True)
    p.add_argument("--V-prime", type=float, default=0.)
    p.add_argument("--mu", type=float, default=0.)
    p.add_argument("--eta", type=float, default=0.02, help="positive retarded broadening")
    p.add_argument("--spin-factor", type=float, default=1., help="spin convention multiplier for bubble")
    p.add_argument("--nk", type=int, default=96, help="points per direction in uniform lattice momentum grid")
    p.add_argument("--nq", type=int, default=32, help="points per direction in output q grid")
    args = p.parse_args()

    data = np.load(args.input)
    omega = np.asarray(data["omega"], float)
    chi_c = complex_array(data["chi_cluster"])
    sigma_c = complex_array(data["sigma_cluster"])
    q_c = np.asarray(data["q_cluster"], float)
    k_c = np.asarray(data["k_cluster"], float)
    if omega.ndim != 1 or len(omega) < 3 or np.any(np.diff(omega) <= 0):
        raise ValueError("omega must be a strictly increasing one-dimensional grid")
    if chi_c.shape != (len(q_c), len(omega)):
        raise ValueError("chi_cluster must have shape (nQ, nw)")
    if sigma_c.shape != (len(k_c), len(omega), 2, 2):
        raise ValueError("sigma_cluster must have shape (nK, nw, 2, 2)")
    if args.eta <= 0 or args.beta <= 0:
        raise ValueError("eta and beta must be positive")

    nk = args.nk
    axis = -np.pi + 2 * np.pi * np.arange(nk) / nk
    k = np.stack(np.meshgrid(axis, axis, indexing="ij"), axis=-1).reshape(-1, 2)
    def lattice_g(kpoints):
        return greens(omega, kpoints, assign_sigma(kpoints, k_c, sigma_c), args)
    qaxis = -np.pi + 2 * np.pi * np.arange(args.nq) / args.nq
    q = np.stack(np.meshgrid(qaxis, qaxis, indexing="ij"), axis=-1).reshape(-1, 2)

    chi0_lat = np.zeros((len(q), len(omega)), dtype=complex)
    for iq, qq in enumerate(q):
        chi0_lat[iq] = bubble_at_q(omega, qq, k, lattice_g, args.beta,
                                   args.spin_factor, args.eta)

    # Compute the cluster bubble directly from supplied G_c(K,w), then infer
    # V_eff(Q,w) = chi0_c^-1 - chi_c^-1 using the scalar-kernel convention.
    chi0_c = np.zeros_like(chi_c)
    if "greens_cluster" not in data:
        raise ValueError("input NPZ must also contain greens_cluster with shape (nK,nw,2,2)")
    g_c = complex_array(data["greens_cluster"])
    if g_c.shape != (len(k_c), len(omega), 2, 2):
        raise ValueError("greens_cluster must have shape (nK,nw,2,2)")

    def cluster_g(kpoints):
        idx = np.argmin(np.sum(wrap_k(kpoints[:, None, :] - k_c[None, :, :]) ** 2, axis=-1), axis=1)
        return g_c[idx]

    for iq, qq in enumerate(q_c):
        chi0_c[iq] = bubble_at_q(omega, qq, k_c, cluster_g, args.beta,
                                 args.spin_factor, args.eta)

    tiny = 1e-14
    v_eff = 1. / (chi0_c + tiny) - 1. / (chi_c + tiny)
    chi_lattice = 1. / (1. / (chi0_lat + tiny) - v_eff[np.argmin(
        np.sum(wrap_k(q[:, None, :] - q_c[None, :, :]) ** 2, axis=-1), axis=1)])
    np.savez_compressed(args.output, omega=omega, q=q, chi_lattice=chi_lattice,
                        chi0_lattice=chi0_lat, chi_cluster=chi_c,
                        chi0_cluster=chi0_c, V_eff_cluster=v_eff)


if __name__ == "__main__":
    main()
