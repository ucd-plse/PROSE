import numpy as np
import sys
import os

log_path = sys.argv[1]
log_path_for_baseline = sys.argv[2]

cost = 0
try:               
    with open(os.path.join(log_path, "outlog.txt"), "r") as f:
        for line in f.readlines():
            if line.startswith("Total runtime"):
                cost += float(line.split()[4])
except:
    pass

file_name = "ocean.stats"
with open("ocean.stats", "r") as f:
    lines = f.readlines()

column_names = [
    "Energy/Mass",
    "Maximum CFL",
    "Mean Sea Level",
    "Total Mass",
    "Mean Salin",
    "Frac Mass Err",
    "Salin Err",
    "Temp Err",
]

rows = []
i = 1
while i + 1 < len(lines):
    i += 1

    row = []
    fields = lines[i].split(",")
    for j in range(4,12):
        _, value = tuple([val.strip() for val in fields[j].split()])
        row.append(float(value))
    rows.append(row)

X = np.array(rows)
with open(os.path.join(log_path, "output.npy"), "wb") as f:
    np.save(f, X)

if log_path.endswith("prose_logs/0000"):
    X_baseline = X
else:
    with open(os.path.join(log_path_for_baseline, "output.npy"), "rb") as f:
        X_baseline = np.load(f)

X_errs = np.where(X_baseline == 0, 0, np.abs((X_baseline - X)/X_baseline))
X_errs_L2 = np.linalg.norm(X_errs, axis=0)
X_errs_inf = np.linalg.norm(X_errs, ord=np.inf, axis=0,)

with open(os.path.join(log_path, "errors_L2.npy"), "wb") as f:
    np.save(f, X_errs_L2)
with open(os.path.join(log_path, "errors_inf.npy"), "wb") as f:
    np.save(f, X_errs_inf)

# check maximum CFL
if X_errs_L2[1] > 0.25:
    cost = -1 * abs(cost)

print(cost)