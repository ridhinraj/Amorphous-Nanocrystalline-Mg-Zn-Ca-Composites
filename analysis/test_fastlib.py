import numpy as np, os, tempfile
import fastlib as fl
from compute_d2min import compute_d2min_frame, parse_dump_frames
from bond_order_params import compute_q4_q6

rng = np.random.default_rng(5)
origin = np.array([137.4, -22.9, 8.1])

print("TEST 1: vectorized q4/q6 vs the validated per-atom loop (random points, non-zero box origin)")
n, L = 2500, 38.0
pos = rng.uniform(0, L, (n, 3)) + origin; box = np.column_stack([origin, origin + L])
for cutoff in (3.6, 4.0):
    q4_ref, q6_ref = compute_q4_q6(pos, box, cutoff)
    i, j, d, r = fl.neighbor_pairs(pos, box, cutoff); q = fl.steinhardt(n, i, d, r)
    same_nan = np.array_equal(np.isnan(q4_ref), np.isnan(q[4]))
    diff = max(np.nanmax(np.abs(q[4] - q4_ref)), np.nanmax(np.abs(q[6] - q6_ref)))
    print(f"   cutoff {cutoff}: NaN pattern identical: {same_nan}; max |difference| over {np.sum(~np.isnan(q4_ref))} atoms = {diff:.2e}")
    assert same_nan and diff < 1e-9
a = 4.5; reps = 7; basis = np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5]])
cells = np.array(np.meshgrid(*[range(reps)]*3)).T.reshape(-1, 3); fcc = (cells[:, None, :] + basis[None]).reshape(-1, 3) * a
fbox = np.array([[0, reps*a]]*3); i, j, d, r = fl.neighbor_pairs(fcc, fbox, a*0.8); q = fl.steinhardt(len(fcc), i, d, r)
print(f"   FCC lattice: q4 = {np.nanmean(q[4]):.3f}, q6 = {np.nanmean(q[6]):.3f}  (literature 0.191 / 0.575)")
assert abs(np.nanmean(q[4]) - 0.191) < 0.002 and abs(np.nanmean(q[6]) - 0.575) < 0.002

print("\nTEST 2: vectorized D2min vs compute_d2min_frame (loop), non-zero origin, strained + one defect + thermal noise")
n, L = 3000, 40.0; ids = np.arange(1, n + 1)
ref = rng.uniform(0, L, (n, 3)) + origin; rbox = np.column_stack([origin, origin + L])
eps, nu = 0.04, 0.3; F = np.diag([1 + eps, 1 - nu*eps, 1 - nu*eps])
cur = (ref - origin) @ F.T + origin + rng.normal(scale=0.03, size=(n, 3))
cbox = np.column_stack([origin, origin + L * np.diag(F)])
cur[700] += rng.normal(scale=2.0, size=3)
_, ref_d2 = compute_d2min_frame(ids, ref, rbox, ids, cur, cbox, cutoff=4.5)
new_d2 = fl.d2min(ref, rbox, cur, cbox, 4.5)
same_nan = np.array_equal(np.isnan(ref_d2), np.isnan(new_d2)); ok = ~np.isnan(ref_d2)
print(f"   NaN pattern identical: {same_nan}; max |difference| = {np.max(np.abs(ref_d2[ok]-new_d2[ok])):.2e}; max D2min = {np.nanmax(new_d2):.3f}")
assert same_nan and np.max(np.abs(ref_d2[ok] - new_d2[ok])) < 1e-8
aff = (ref - origin) @ F.T + origin
z = fl.d2min(ref, rbox, aff, cbox, 4.5); print(f"   pure affine strain: max D2min = {np.nanmax(z):.2e} (should be ~0)"); assert np.nanmax(z) < 1e-8

print("\nTEST 3: frame-selective reader vs parse_dump_frames (atoms deliberately written in shuffled order)")
tmp = tempfile.mkdtemp(); path = os.path.join(tmp, "t.dump"); n = 800
with open(path, "w") as f:
    for ts in (20000, 25000, 30000):
        p = rng.uniform(0, 30, (n, 3)) + origin; perm = rng.permutation(n)
        f.write(f"ITEM: TIMESTEP\n{ts}\nITEM: NUMBER OF ATOMS\n{n}\nITEM: BOX BOUNDS pp pp pp\n")
        f.write(f"{origin[0]} {origin[0]+30}\n{origin[1]} {origin[1]+30}\n{origin[2]} {origin[2]+30}\n")
        f.write("ITEM: ATOMS id type x y z c_peratom_stress[1] c_peratom_stress[2] c_peratom_stress[3] c_peratom_stress[4] c_peratom_stress[5] c_peratom_stress[6]\n")
        for k in perm: f.write(f"{k+1} {1+k%3} {p[k,0]:.5f} {p[k,1]:.5f} {p[k,2]:.5f} 1 2 3 4 5 {ts/1000:.3f}\n")
idx = fl.index_dump(path); old = parse_dump_frames(path)
print("   frames indexed:", [fr["timestep"] for fr in idx])
for k in range(3):
    fr = fl.read_frame(path, idx[k]); assert fr["timestep"] == old[k][0] and np.allclose(fr["pos"], old[k][2]) and np.array_equal(fr["id"].astype(int), old[k][1])
    assert np.allclose(fr["c_peratom_stress[6]"], idx[k]["timestep"]/1000)
info, s, k = fl.frame_at_strain(idx, 0.0053, erate_per_ps=0.001); print(f"   frame_at_strain(0.0053) -> timestep {info['timestep']} (strain {s:.4f}) as expected: {info['timestep']==25000}")
assert info["timestep"] == 25000
print("\nALL FASTLIB TESTS PASSED")
