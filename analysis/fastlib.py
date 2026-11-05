"""
fastlib.py -- vectorized building blocks shared by combined_regression.py, d2min_sensitivity.py and
interphase_all_seeds.py.  Every function here is checked against the project's earlier, validated per-atom
loop implementations (compute_d2min.compute_d2min_frame, bond_order_params.compute_q4_q6,
multiframe_solute_shap.build_coordination) by test_fastlib.py before it is used.

Why it exists: the loop versions take 30-40 s per frame pair at 390k atoms; these take a few seconds, which
makes parameter-sensitivity scans and all-seed analyses affordable.  Conventions are identical to the
earlier tools: orthogonal periodic box with a possibly non-zero origin, neighbor lists built on the
REFERENCE frame, minimum-image displacements, atoms sorted by id.
"""
import itertools, io
import numpy as np
from scipy.spatial import cKDTree
from scipy.special import sph_harm_y

# ------------------------------------------------------------------ dump reading (frame-selective)
def index_dump(path):
    """One pass over a LAMMPS dump: [{timestep, n, box, cols, offset}], offset = byte position of the atom block."""
    frames = []
    with open(path, "rb") as f:
        while True:
            line = f.readline()
            if not line: break
            if not line.startswith(b"ITEM: TIMESTEP"): continue
            ts = int(f.readline()); f.readline(); n = int(f.readline()); f.readline()
            box = np.array([[float(x) for x in f.readline().split()[:2]] for _ in range(3)])
            cols = f.readline().decode().replace("ITEM: ATOMS", "").split()
            off = f.tell()
            for _ in range(n): f.readline()
            frames.append(dict(timestep=ts, n=n, box=box, cols=cols, offset=off))
    return frames

def read_frame(path, info):
    """Parse ONE frame (all columns) -> dict(cols->array), box, timestep; atoms sorted by id."""
    with open(path, "rb") as f:
        f.seek(info["offset"]); data = b"".join(itertools.islice(f, info["n"]))
    try: arr = np.fromstring(data.decode(), dtype=float, sep=" ")
    except Exception: arr = np.loadtxt(io.BytesIO(data))
    arr = arr.reshape(info["n"], len(info["cols"]))
    arr = arr[np.argsort(arr[:, info["cols"].index("id")])]
    out = {c: arr[:, k] for k, c in enumerate(info["cols"])}
    out["pos"] = np.column_stack([out["x"], out["y"], out["z"]]); out["type"] = out["type"].astype(int)
    out["box"] = info["box"]; out["timestep"] = info["timestep"]
    return out

def frame_at_strain(index, strain, erate_per_ps=0.001, t0=20000, dt_ps=0.001):
    """nearest recorded frame to a strain; returns (info, actual_strain)"""
    s = np.array([(fr["timestep"] - t0) * erate_per_ps * dt_ps for fr in index])
    k = int(np.argmin(np.abs(s - strain))); return index[k], s[k], k

# ------------------------------------------------------------------ neighbors
def neighbor_pairs(pos, box, cutoff):
    """DIRECTED neighbor pairs (i, j) within cutoff, with minimum-image displacement d = r_j - r_i and distance r."""
    L = box[:, 1] - box[:, 0]; w = (pos - box[:, 0]) % L
    pr = cKDTree(w, boxsize=L).query_pairs(cutoff, output_type="ndarray")
    i = np.concatenate([pr[:, 0], pr[:, 1]]); j = np.concatenate([pr[:, 1], pr[:, 0]])
    d = w[j] - w[i]; d -= L * np.round(d / L)
    return i, j, d, np.linalg.norm(d, axis=1)

# ------------------------------------------------------------------ Steinhardt q_l
def steinhardt(n_atoms, i, d, r, ls=(4, 6), min_neighbors=4):
    """q_l for every atom that appears as a center in the (already filtered) directed pair list."""
    N = np.bincount(i, minlength=n_atoms)
    theta = np.arccos(np.clip(d[:, 2] / r, -1.0, 1.0)); phi = np.arctan2(d[:, 1], d[:, 0])
    out = {}
    for l in ls:
        s2 = np.zeros(n_atoms)
        for m in range(-l, l + 1):
            Y = sph_harm_y(l, m, theta, phi)
            s2 += np.bincount(i, weights=Y.real, minlength=n_atoms) ** 2 + np.bincount(i, weights=Y.imag, minlength=n_atoms) ** 2
        q = np.sqrt(4 * np.pi / (2 * l + 1) * s2) / np.maximum(N, 1)
        q[N < min_neighbors] = np.nan; out[l] = q
    return out

# ------------------------------------------------------------------ D2_min (Falk & Langer), closed form
def d2min(ref_pos, ref_box, cur_pos, cur_box, cutoff, min_neighbors=4):
    """
    D2_min for every atom, reference-frame neighbor list.  With X = sum dcur (x) dref and Y = sum dref (x) dref,
    the best affine F = X Y^-1 and  D2 = [ sum|dcur|^2 - tr(X Y^-1 X^T) ] / n_neighbors  (algebraically identical to
    the residual of the least-squares fit used by compute_d2min.compute_d2min_frame).
    """
    n = len(ref_pos); i, j, dref, _ = neighbor_pairs(ref_pos, ref_box, cutoff)
    Lc = cur_box[:, 1] - cur_box[:, 0]
    dcur = cur_pos[j] - cur_pos[i]; dcur -= Lc * np.round(dcur / Lc)
    cnt = np.bincount(i, minlength=n)
    X = np.empty((n, 3, 3)); Y = np.empty((n, 3, 3))
    for a in range(3):
        for b in range(3):
            X[:, a, b] = np.bincount(i, weights=dcur[:, a] * dref[:, b], minlength=n)
            Y[:, a, b] = np.bincount(i, weights=dref[:, a] * dref[:, b], minlength=n)
    sc = np.bincount(i, weights=(dcur ** 2).sum(1), minlength=n)
    out = np.full(n, np.nan); ok = np.where(cnt >= min_neighbors)[0]
    det = np.linalg.det(Y[ok]); ok = ok[np.abs(det) > 1e-9]
    Yi = np.linalg.inv(Y[ok]); t = np.einsum("nab,nbc,nac->n", X[ok], Yi, X[ok])
    out[ok] = np.maximum((sc[ok] - t) / cnt[ok], 0.0)
    return out
