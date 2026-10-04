"""Relay-BP decoding with per-shot work statistics (requires the synbuf-patched relay_bp)."""
import numpy as np
import relay_bp

from .dem import dem_to_matrices

# Columns returned by ObservableDecoderRunner.decode_stats_batch (synbuf patch)
COL_CONVERGED, COL_ITERS, COL_FIRST_ITERS, COL_LEGS, COL_OBS_BEST, COL_OBS_FIRST = range(6)


def build_relay(dem, gamma0=0.125, gamma_dist_interval=(-0.24, 0.66), pre_iter=80, num_sets=100000,
                set_max_iter=60, stop_nconv=1, alpha=None, seed=0, dtype="f32", keep_detectors=None):
    H, L, priors = dem_to_matrices(dem, keep_detectors)
    cls = {"f32": relay_bp.RelayDecoderF32, "f64": relay_bp.RelayDecoderF64}[dtype]
    kwargs = dict(error_priors=priors, gamma0=gamma0, pre_iter=pre_iter, num_sets=num_sets,
                  set_max_iter=set_max_iter, gamma_dist_interval=tuple(gamma_dist_interval),
                  stop_nconv=stop_nconv, seed=seed)
    if alpha is not None:
        kwargs["alpha"] = alpha
    dec = cls(H.tocsr(), **kwargs)
    runner = relay_bp.ObservableDecoderRunner(dec, L.tocsr())
    return runner, H, L, priors


def pack_obs(obs_bits):
    """Pack an (N, k<=64) bool/uint8 observable array into uint64 masks (bit i = observable i)."""
    obs_bits = np.asarray(obs_bits, dtype=np.uint64)
    w = (np.uint64(1) << np.arange(obs_bits.shape[1], dtype=np.uint64))
    return (obs_bits * w).sum(axis=1).astype(np.uint64)
