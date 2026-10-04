"""Detector-error-model helpers."""
import numpy as np
import scipy.sparse as sp
import stim


def dem_to_matrices(dem: stim.DetectorErrorModel, keep_detectors=None):
    """Convert a DEM (no decomposition) into (H, L, priors).

    H: detectors x errors (csr, uint8), L: observables x errors (csr, uint8), priors: float64.
    Error mechanisms with identical (detectors, observables) signatures are merged.
    If keep_detectors (sorted index array) is given, other detectors are dropped and H rows are
    renumbered accordingly; mechanisms that become trivial are discarded.
    """
    remap = None
    if keep_detectors is not None:
        remap = {int(d): i for i, d in enumerate(keep_detectors)}
    sig = {}
    for inst in dem.flattened():
        if inst.type != "error":
            continue
        p = inst.args_copy()[0]
        dets, obs = [], []
        for t in inst.targets_copy():
            if t.is_relative_detector_id():
                dets.append(t.val)
            elif t.is_logical_observable_id():
                obs.append(t.val)
        if remap is not None:
            dets = [remap[d] for d in dets if d in remap]
        if not dets and not obs:
            continue
        key = (tuple(sorted(set(d for d in dets if dets.count(d) % 2))),
               tuple(sorted(set(o for o in obs if obs.count(o) % 2))))
        if key in sig:
            q = sig[key]
            sig[key] = q * (1 - p) + p * (1 - q)
        else:
            sig[key] = p
    keys = list(sig)
    priors = np.array([sig[k] for k in keys], dtype=np.float64)
    hr, hc, lr, lc = [], [], [], []
    for j, (dets, obs) in enumerate(keys):
        hr += dets
        hc += [j] * len(dets)
        lr += obs
        lc += [j] * len(obs)
    nd = dem.num_detectors if remap is None else len(remap)
    H = sp.csr_matrix((np.ones(len(hr), dtype=np.uint8), (hr, hc)), shape=(nd, len(keys)))
    L = sp.csr_matrix((np.ones(len(lr), dtype=np.uint8), (lr, lc)), shape=(dem.num_observables, len(keys)))
    return H, L, priors


def poisson_binomial_pmf(probs, nmax=None):
    """Exact PMF of the number of independent Bernoulli(p_i) events that fire.

    Computed by iterated convolution, truncated at nmax (the truncated tail mass is returned
    as the last entry's complement, i.e. sum(pmf) <= 1).
    """
    probs = np.asarray(probs, dtype=np.float64)
    if nmax is None:
        nmax = len(probs)
    pmf = np.zeros(nmax + 1)
    pmf[0] = 1.0
    for p in probs:
        pmf[1:] = pmf[1:] * (1 - p) + pmf[:-1] * p
        pmf[0] *= (1 - p)
    return pmf


def z_type_detectors(circuit: stim.Circuit):
    """Indices of detectors insensitive to every Z error (i.e. Z-check detectors of a CSS circuit)."""
    noiseless = circuit.without_noise().flattened()
    probe = stim.Circuit()
    allq = list(range(circuit.num_qubits))
    for inst in noiseless:
        probe.append(inst)
        if inst.name not in ("DETECTOR", "OBSERVABLE_INCLUDE", "TICK", "QUBIT_COORDS", "SHIFT_COORDS"):
            probe.append("Z_ERROR", allq, 0.001)
    dem = probe.detector_error_model(allow_gauge_detectors=True)
    sensitive = np.zeros(circuit.num_detectors, dtype=bool)
    for inst in dem.flattened():
        if inst.type == "error":
            for t in inst.targets_copy():
                if t.is_relative_detector_id():
                    sensitive[t.val] = True
    return np.nonzero(~sensitive)[0]
