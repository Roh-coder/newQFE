import json
import glob
import os

results_dir = "K_from_continuum/results/production_iso111_qtr_tests_L128_grid7x7_20260515/test_data/grid"
meta_files = glob.glob(os.path.join(results_dir, "**/*.meta.json"), recursive=True)

times_by_L = {}
for f in meta_files:
    with open(f, 'r') as jfile:
        data = json.load(jfile)
        L = data['L']
        wall_time = data.get('wall_seconds', 0)
        times_by_L.setdefault(L, []).append(wall_time)

# Expected grid: 7x7 = 49 points.
# L values are 16, 32, 48, 64, 80, 96, 112, 128 (8 sizes based on meta.json observed)
# Total intended points = 49 * 8 = 392 points? 
# Or 49 points total across different L?
# The prompt says "grid7x7" and "L128". 
# Usually it's 49 points PER L, but let's check directory structure.

L_dirs = glob.glob(os.path.join(results_dir, "L*"))
L_values = sorted([int(os.path.basename(d)[1:]) for d in L_dirs])

print(f"L values found: {L_values}")
for L in L_values:
    avg = sum(times_by_L.get(L, [0])) / max(len(times_by_L.get(L, [1])), 1)
    print(f"L={L}: {len(times_by_L.get(L, []))} completed, Avg wall_time: {avg:.2f}s")

# Total points in 7x7 grid is 49. Assuming 49 points for each of the 8 L's.
# Total points = 8 * 49 = 392.
# Completed = 16.
# Remaining = 376.

# Let's refine the estimation based on scaling.
