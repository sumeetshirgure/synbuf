"""Stim circuits for the experiments.

* Rotated surface code memory (Stim's generator, uniform circuit-level depolarizing noise).
* Bivariate bicycle (BB) code memory using the depth-8 syndrome cycle of
  Bravyi et al., Nature 627, 778 (2024), with uniform circuit-level noise
  (the same "uniform_circuit" model as IBM's Relay-BP test circuits: prep/measure
  flips with prob p, DEPOLARIZE2(p) after every CNOT, DEPOLARIZE1(p) on idle data qubits).
"""
from dataclasses import dataclass

import numpy as np
import stim

from .gf2 import logical_ops, rank


def surface_code_circuit(d: int, p: float, rounds: int | None = None, basis: str = "x") -> stim.Circuit:
    return stim.Circuit.generated(
        f"surface_code:rotated_memory_{basis}",
        distance=d,
        rounds=d if rounds is None else rounds,
        after_clifford_depolarization=p,
        before_round_data_depolarization=p,
        before_measure_flip_probability=p,
        after_reset_flip_probability=p,
    )


# ---------------------------------------------------------------- BB codes
@dataclass(frozen=True)
class BBCode:
    name: str
    l: int
    m: int
    A: tuple  # three monomials, each ('x', power) or ('y', power)
    B: tuple
    n: int
    k: int
    d: int


BB_CODES = {
    "72_12_6": BBCode("72_12_6", 6, 6, (("x", 3), ("y", 1), ("y", 2)), (("y", 3), ("x", 1), ("x", 2)), 72, 12, 6),
    "90_8_10": BBCode("90_8_10", 15, 3, (("x", 9), ("y", 1), ("y", 2)), (("x", 0), ("x", 2), ("x", 7)), 90, 8, 10),
    "108_8_10": BBCode("108_8_10", 9, 6, (("x", 3), ("y", 1), ("y", 2)), (("y", 3), ("x", 1), ("x", 2)), 108, 8, 10),
    "144_12_12": BBCode("144_12_12", 12, 6, (("x", 3), ("y", 1), ("y", 2)), (("y", 3), ("x", 1), ("x", 2)), 144, 12, 12),
    "288_12_18": BBCode("288_12_18", 12, 12, (("x", 3), ("y", 2), ("y", 7)), (("y", 3), ("x", 1), ("x", 2)), 288, 12, 18),
}


def _shift(size, power):
    return np.roll(np.eye(size, dtype=np.uint8), power % size, axis=1)


def bb_matrices(code: BBCode):
    """Return ([A1,A2,A3], [B1,B2,B3], HX, HZ) as dense uint8 arrays."""
    l, m = code.l, code.m
    Il, Im = np.eye(l, dtype=np.uint8), np.eye(m, dtype=np.uint8)

    def mono(v, pw):
        if v == "x":
            return np.kron(_shift(l, pw), Im)
        return np.kron(Il, _shift(m, pw))

    As = [mono(*t) for t in code.A]
    Bs = [mono(*t) for t in code.B]
    A = (As[0] + As[1] + As[2]) % 2
    B = (Bs[0] + Bs[1] + Bs[2]) % 2
    HX = np.hstack([A, B]).astype(np.uint8)
    HZ = np.hstack([B.T, A.T]).astype(np.uint8)
    return As, Bs, HX, HZ


def bb_code_parameters(code: BBCode):
    _, _, HX, HZ = bb_matrices(code)
    n = HX.shape[1]
    k = n - rank(HX) - rank(HZ)
    return n, k


def bb_memory_circuit(code: BBCode, p: float, rounds: int | None = None) -> stim.Circuit:
    """Z-basis memory experiment (detects X errors) with Bravyi et al.'s depth-8 cycle.

    Only Z-check detectors are included (standard for CSS memory decoding of one basis).
    Observables: a basis of k logical Z operators.
    """
    rounds = code.d if rounds is None else rounds
    As, Bs, HX, HZ = bb_matrices(code)
    n2 = code.l * code.m
    Xc = list(range(0, n2))
    L = list(range(n2, 2 * n2))
    R = list(range(2 * n2, 3 * n2))
    Zc = list(range(3 * n2, 4 * n2))
    data = L + R
    mats = As + Bs  # X checks: index 0..2 -> A1..A3 on L, 3..5 -> B1..B3 on R

    sX = ["idle", 1, 4, 3, 5, 0, 2]
    sZ = [3, 5, 0, 1, 2, 4, "idle"]

    def x_cnots(k):
        M = mats[k]
        side = L if k < 3 else R
        pairs = []
        for i in range(n2):
            j = int(np.nonzero(M[i])[0][0])
            pairs += [Xc[i], side[j]]
        return pairs

    zmats = Bs + As  # Z checks: index 0..2 -> B1..B3^T on L, 3..5 -> A1..A3^T on R (HZ = [B^T | A^T])

    def z_cnots(k):
        M = zmats[k]
        side = L if k < 3 else R
        pairs = []
        for i in range(n2):
            j = int(np.nonzero(M[:, i])[0][0])  # (M^T)[i, j] = 1
            pairs += [side[j], Zc[i]]
        return pairs

    c = stim.Circuit()
    # initial preparation: data in |0>, Z checks in |0>
    c.append("R", data + Zc)
    c.append("X_ERROR", data + Zc, p)
    c.append("TICK")

    num_meas = 0
    zrec = {}  # cycle -> measurement index of first Z check

    def cycle(ci):
        nonlocal num_meas
        for t in range(7):
            busy = set()
            if sX[t] == "idle":
                c.append("RX", Xc)
                c.append("Z_ERROR", Xc, p)
            else:
                pr = x_cnots(sX[t])
                c.append("CX", pr)
                c.append("DEPOLARIZE2", pr, p)
                busy.update(pr)
            if sZ[t] == "idle":
                c.append("M", Zc, p)
                zrec[ci] = num_meas
                num_meas += n2
            else:
                pr = z_cnots(sZ[t])
                c.append("CX", pr)
                c.append("DEPOLARIZE2", pr, p)
                busy.update(pr)
            idle = [q for q in data if q not in busy]
            if idle:
                c.append("DEPOLARIZE1", idle, p)
            c.append("TICK")
        # round 8: measure X checks, prepare Z checks, data idle
        c.append("MX", Xc, p)
        num_meas += n2
        c.append("R", Zc)
        c.append("X_ERROR", Zc, p)
        c.append("DEPOLARIZE1", data, p)
        c.append("TICK")
        # Z-check detectors
        for i in range(n2):
            cur = zrec[ci] + i
            recs = [stim.target_rec(cur - num_meas)]
            if ci > 0:
                recs.append(stim.target_rec(zrec[ci - 1] + i - num_meas))
            c.append("DETECTOR", recs, [i, ci])

    for ci in range(rounds):
        cycle(ci)

    # final data measurement
    c.append("M", data, p)
    data_base = num_meas
    num_meas += 2 * n2
    for i in range(n2):
        supp = np.nonzero(HZ[i])[0]
        recs = [stim.target_rec(data_base + int(q) - num_meas) for q in supp]
        recs.append(stim.target_rec(zrec[rounds - 1] + i - num_meas))
        c.append("DETECTOR", recs, [i, rounds])
    LZ = logical_ops(HZ, HX)
    assert LZ.shape[0] == code.k, (LZ.shape, code.k)
    for o, row in enumerate(LZ):
        recs = [stim.target_rec(data_base + int(q) - num_meas) for q in np.nonzero(row)[0]]
        c.append("OBSERVABLE_INCLUDE", recs, o)
    return c
