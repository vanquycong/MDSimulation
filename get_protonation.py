import glob
import numpy as np
import matplotlib.pyplot as plt
import re
import csv

def moving_average_with_std(x, y, window_size):
    x_avg = []
    y_avg = []
    y_std = []
    half_window = window_size // 2
    for i in range(half_window, len(y) - half_window):
        window_y = y[i - half_window:i + half_window]
        window_x = x[i - half_window:i + half_window]
        x_avg.append(np.mean(window_x))
        y_avg.append(np.mean(window_y))
        y_std.append(np.std(window_y))
    return np.array(x_avg), np.array(y_avg), np.array(y_std)
    
def parse_ires_from_file(fname):
    with open(fname, 'r') as f:
        for line in f:
            if line.startswith('#  ires'):
                return list(map(int, line.split()[2:]))
    raise ValueError(f"Cannot find # ires line in {fname}")

def read_lambda_data(files, target_residues_dict, time_per_step_ns=0.000002):
    target_residues = list(target_residues_dict.keys())
    # Parse residues in columns from first file
    ires = parse_ires_from_file(files[0])

    # Find column indices for each residue's two lambda columns
    residue_indices = {}
    for res in target_residues:
        idxs = [i for i, r in enumerate(ires) if r == res]
        if len(idxs) != 2:
            raise ValueError(f"Residue {res} does not have exactly two lambda columns!")
        residue_indices[res] = idxs

    # Storage for all data
    
    data = {res: {'time': [], 'prot': [], 'side': []} for res in target_residues}
    cumulative_time = 0.0
    for fname in files:
        with open(fname, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split()
                step = int(fields[0])
                time = cumulative_time + step * time_per_step_ns

                for res in target_residues:
                    idx1, idx2 = residue_indices[res]
                    lam1 = float(fields[1 + idx1])
                    lam2 = float(fields[1 + idx2])

                    data[res]['time'].append(time)
                    data[res]['prot'].append(lam1)
                    data[res]['side'].append(lam2)

        cumulative_time += 5_000_000 * time_per_step_ns  # 10 ns per file

    return data

def parse_itauto_from_file(fname):
    with open(fname, 'r') as f:
        for line in f:
            if line.startswith('# itauto'):
                return list(map(int, line.split()[2:]))
    raise ValueError(f"Cannot find '# itauto' line in {fname}")
    
def get_protonation_state_from_data(data, residues_dict, target_time_ns):
    result = {}
    intermediates = []

    for resid, vals in data.items():
        times = vals['time']
        prots = vals['prot']
        
        # Find closest frame to target_time_ns
        idx = min(range(len(times)), key=lambda i: abs(times[i] - target_time_ns))
        lam = prots[idx]
        
        if lam < 0.2:
            result[resid] = 'protonated'
        elif lam > 0.8:
            result[resid] = 'deprotonated'
        else:
            result[resid] = 'intermediate'
            intermediates.append(resid)

    return result, intermediates

def get_protonation_table(data, residues_dict, itauto, target_time_ns, output_csv="protonation_state.csv"):
    table = []
    intermediates = []
    
    # Flatten the mapping from residue -> column indices
    # This duplicates logic from read_lambda_data
    ires = [r for r in data.keys()]
    residue_indices = {}
    all_residues = list(data.keys())

    flat_index = 0
    # Build mapping from residue to its two column indices
    for res in residues_dict:
        # Find 2 columns for that residue
        count = 0
        indices = []
        for i, r in enumerate(parse_ires_from_file(lambda_files[0])):
            if r == res:
                indices.append(i)
                count += 1
        if count != 2:
            raise ValueError(f"Residue {res} does not have exactly 2 λ columns.")
        residue_indices[res] = indices

    for resid, vals in data.items():
        times = vals['time']
        prots = vals['prot']
        idx = min(range(len(times)), key=lambda i: abs(times[i] - target_time_ns))
        lam = prots[idx]

        idx1, idx2 = residue_indices[resid]
        lam1 = data[resid]['prot'][idx]
        lam2 = data[resid]['side'][idx]
        site1 = itauto[idx1]
        site2 = itauto[idx2]

        print(resid, residues_dict[resid], lam1, lam2, site1, site2)
        if lam1 < 0.2 and lam2 > 0.8:
            state = "protonated"
            site = site1
        if lam1 < 0.2 and lam2 < 0.2:
            state = "protonated"
            site = site2
        if lam1 < 0.2 and 0.2 <= lam2 <=0.8:
            state = "protonated"
            site = site1 
        #elif lam1 > 0.8:
        #    state = "deprotonated"
        #    site = "-"
        elif lam1 > 0.8 and lam2 < 0.2:
            state = "deprotonated"
            site = site1 
        elif lam1 > 0.8 and lam2 > 0.8:
            state = "deprotonated"
            site = site2 
        elif lam1 > 0.8 and 0.2 <= lam2 <=0.8:
            state = "deprotonated"
            site = site1 
        #else:
        elif 0.2 <= lam1 <= 0.8:
            state = "intermediate"
            site = "-"
            intermediates.append(resid)

        table.append({
            "Residue": resid,
            "Name": residues_dict[resid],
            "State": state,
            "Protonation Site": site
        })

    # Save to CSV
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["Residue", "Name", "State", "Protonation Site"])
        writer.writeheader()
        writer.writerows(table)

    return table, intermediates
    

# 1. Load data and metadata
lambda_files = sorted(glob.glob("../../protonation/pH5.0_*.lambda"), key=lambda x: int(re.findall(r"_(\d+)\.lambda", x)[0]))
residues_to_plot = {31: "Asp", 39: "His", 47: "Glu", 54: "Asp", 56: "Glu", 62: "His", 110: "His", 125: "Glu", 133: "His", 
                    144: "Glu", 200: "Glu", 205: "His", 210: "Asp", 211: "Glu", 212: "Glu", 217: "Glu", 218: "Glu", 222: "Glu", 
                    231: "Glu", 272: "Glu", 273: "Glu", 278: "Glu", 290: "His", 305: "Asp", 330: "Glu", 341: "His", 368: "His", 
                    370: "Asp", 408: "Glu", 409: "Glu", 421: "His", 435: "Glu", 441: "Asp"}

data = read_lambda_data(lambda_files, residues_to_plot)
itauto = parse_itauto_from_file(lambda_files[0])

# 2. Analyze and save
frame_time = 470.0
table, intermediates = get_protonation_table(data, residues_to_plot, itauto, frame_time)

# 3. Optional print
print(f"Saved protonation states at {frame_time} ns to 'protonation_state.csv'.")
if intermediates:
    print("Residues with intermediate protonation:")
    for r in intermediates:
        print(f"  Residue {r} ({residues_to_plot[r]})")