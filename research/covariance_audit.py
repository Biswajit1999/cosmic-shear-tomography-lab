"""Single-pair KiDS xi+ covariance sensitivity and nuisance-offset diagnostics.

These are conditional diagnostics, not a combined tomographic likelihood.
"""
import math

def solve(c, b):
    n = len(b)
    if n == 0 or len(c) != n or any(len(row) != n for row in c):
        raise ValueError("dimension mismatch")
    m = [list(map(float, row)) + [float(v)] for row, v in zip(c, b)]
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(m[i][k]))
        if not math.isfinite(m[p][k]) or abs(m[p][k]) < 1e-25:
            raise ValueError("singular matrix")
        m[k], m[p] = m[p], m[k]
        d = m[k][k]
        for j in range(k, n+1):
            m[k][j] /= d
        for i in range(n):
            if i != k:
                f = m[i][k]
                for j in range(k, n+1):
                    m[i][j] -= f*m[k][j]
    return [row[-1] for row in m]

def validate(x, c):
    n = len(x)
    if not n or len(c) != n or any(len(row) != n for row in c):
        raise ValueError("dimension mismatch")
    if not all(math.isfinite(v) for v in x) or not all(math.isfinite(v) for row in c for v in row):
        raise ValueError("nonfinite input")
    l = [[0.0]*n for _ in range(n)]
    for i in range(n):
        if any(abs(c[i][j]-c[j][i]) > 1e-9*max(abs(c[i][j]), abs(c[j][i]), 1e-30) for j in range(n)):
            raise ValueError("asymmetric covariance")
        for j in range(i+1):
            v = c[i][j]-sum(l[i][k]*l[j][k] for k in range(j))
            if i == j:
                if v <= 0:
                    raise ValueError("non-positive covariance")
                l[i][j] = math.sqrt(v)
            else:
                l[i][j] = v/l[j][j]

def covariance_shrinkage_null(x, c, fractions=(0.0, 0.25, 0.5, 1.0)):
    """Zero-shear chi-square, interpolating released C toward diag(C)."""
    validate(x, c)
    n = len(x)
    out = []
    for a in fractions:
        if not math.isfinite(a) or not 0 <= a <= 1:
            raise ValueError("fraction outside [0,1]")
        ca = [[(1-a)*c[i][j] + (a*c[i][i] if i == j else 0.0) for j in range(n)] for i in range(n)]
        out.append((a, sum(v*w for v, w in zip(x, solve(ca, x)))))
    return out

def profile_constant_offset(x, c):
    """GLS constant nuisance offset, sigma and profiled residual chi-square."""
    validate(x, c)
    w1 = solve(c, [1.0]*len(x))
    wx = solve(c, x)
    den = sum(w1)
    a = sum(wx)/den
    residual = [v-a for v in x]
    return {"offset": a, "sigma": math.sqrt(1/den),
            "residual_chi2": max(0.0, sum(v*w for v,w in zip(residual, solve(c, residual)))),
            "dof": len(x)-1}
