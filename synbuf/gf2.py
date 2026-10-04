"""Small GF(2) linear-algebra helpers (dense uint8 matrices)."""
import numpy as np


def rref(M):
    """Row-reduce M over GF(2). Returns (R, pivot_columns)."""
    R = (np.array(M, dtype=np.uint8) & 1).copy()
    rows, cols = R.shape
    pivots = []
    r = 0
    for c in range(cols):
        if r >= rows:
            break
        nz = np.nonzero(R[r:, c])[0]
        if len(nz) == 0:
            continue
        p = r + nz[0]
        if p != r:
            R[[r, p]] = R[[p, r]]
        others = np.nonzero(R[:, c])[0]
        others = others[others != r]
        R[others] ^= R[r]
        pivots.append(c)
        r += 1
    return R[:r], pivots


def rank(M):
    return len(rref(M)[1])


def nullspace(M):
    """Basis (rows) of the right null space of M over GF(2)."""
    M = np.array(M, dtype=np.uint8) & 1
    n = M.shape[1]
    R, piv = rref(M)
    free = [c for c in range(n) if c not in set(piv)]
    basis = []
    for f in free:
        v = np.zeros(n, dtype=np.uint8)
        v[f] = 1
        for i, p in enumerate(piv):
            if R[i, f]:
                v[p] = 1
        basis.append(v)
    return np.array(basis, dtype=np.uint8).reshape(-1, n)


def logical_ops(H_same, H_other):
    """Logical operators of one type for a CSS code.

    Returns rows spanning ker(H_other) modulo rowspace(H_same), e.g. Z-logicals are
    vectors in ker(H_X) that are not in the row space of H_Z.
    """
    K = nullspace(H_other)
    base = np.array(H_same, dtype=np.uint8) & 1
    r0 = rank(base)
    logicals = []
    cur = base.copy()
    for v in K:
        trial = np.vstack([cur, v[None, :]])
        if rank(trial) > r0 + len(logicals):
            logicals.append(v)
            cur = trial
    return np.array(logicals, dtype=np.uint8)
