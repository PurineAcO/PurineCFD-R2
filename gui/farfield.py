"""来流条件换算: Sutherland 粘度 + (Ma, Re, T, L) -> 静压."""

import math

GAMMA = 1.4
R = 287.05


def sutherland(T, mu0=1.716e-05, t0=273.15, ts=110.4):
    return mu0 * (T / t0) ** 1.5 * (t0 + ts) / (T + ts)


def from_re(Ma, Re, T, L):
    mu = sutherland(T)
    a = math.sqrt(GAMMA * R * T)
    U = Ma * a
    rho = Re * mu / (U * L)
    p = rho * R * T
    return dict(Ma=Ma, T=T, p=p, mu=mu, a=a, U=U, rho=rho, q=0.5 * rho * U * U, Re=Re, L=L)


def from_p(Ma, p, T, L):
    mu = sutherland(T)
    a = math.sqrt(GAMMA * R * T)
    U = Ma * a
    rho = p / (R * T)
    Re = rho * U * L / mu
    return dict(Ma=Ma, T=T, p=p, mu=mu, a=a, U=U, rho=rho, q=0.5 * rho * U * U, Re=Re, L=L)
