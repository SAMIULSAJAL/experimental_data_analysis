import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import os

# Base directory
base_dir = "/Users/ssajal/Desktop/academic/Atomic_and_molecular_lab/M_1.8_Mikrowellen_Spektroskopie/our_lab/code"
data_path = os.path.join(base_dir, "3,3.csv")

# Load data
data = np.loadtxt(data_path, delimiter='\t')
x = data[:, 0]
y = data[:, 1]

# Normalize
y_norm = y - np.min(y)

# Find peaks
peaks, properties = find_peaks(y_norm, distance=20, height=np.max(y_norm)*0.05)

# Filter peaks by cutoff
y_cutoff = 3
filtered_peaks = [p for p in peaks if y[p] >= y_cutoff]

if len(filtered_peaks) < 5:
    print(f"⚠️ Warning: fewer than 5 peaks found after filtering with cutoff {y_cutoff}.")

# Sort and select first 5
sorted_peaks = np.array(filtered_peaks)[np.argsort(x[filtered_peaks])]
selected_peaks = sorted_peaks[:5]

peak_labels = ["Left peak", "Inner left peak", "Centre peak", "Inner right peak", "Right peak"]
colors = ['red', 'red', 'blue', 'green', 'green']
window_width = 20

plt.figure(figsize=(12, 7))

for i, peak in enumerate(selected_peaks):
    peak_x = x[peak]
    mask = (x >= peak_x - window_width) & (x <= peak_x + window_width)
    if i == 2:  # center peak: no vertical line, only scatter + colored region line
        plt.plot(x[mask], y[mask], color=colors[i], alpha=0.7, linewidth=1)
        plt.scatter(peak_x, y[peak], color=colors[i], zorder=5, label=f"Peak {i+1}")
    else:
        plt.plot(x[mask], y[mask], color=colors[i], alpha=0.7, linewidth=1, label=f"Peak {i+1}")
        plt.axvline(peak_x, linestyle='--', color=colors[i], alpha=0.5)
        plt.scatter(peak_x, y[peak], color=colors[i], zorder=5)

plt.plot(x, y, color='black', alpha=0.9, linewidth=2, label='_nolegend_', zorder=10)

plt.xlabel("Frequency (MHz)", fontsize=20)
plt.ylabel("Intensity", fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(fontsize=20)
plt.grid(True)
plt.tight_layout()

plot_path = os.path.join(base_dir, "peaks_plot.png")
plt.savefig(plot_path, dpi=300)
plt.show()


# Calculate ΔE = |ν0 − νi| where ν0 is center peak frequency
v0 = x[selected_peaks[2]]

table_header = r"""\begin{table}[h]
    \centering
    \begin{tabular}{|c|c|c|}
        \hline
        Peak & Frequency (MHz) & $\Delta E$ (MHz) \\
        \hline
"""

table_footer = r"""    \hline
    \end{tabular}
    \caption{Peak centres and energy differences in the 3,3 ammonia transition}
    \label{tab:peak_centres}
\end{table}
"""

table_rows = ""
for label, peak in zip(peak_labels, selected_peaks):
    freq = x[peak]
    delta_e = abs(v0 - freq)
    table_rows += f"        {label} & {freq:.3f} & {delta_e:.3f} \\\\\n"

table_text = table_header + table_rows + table_footer

table_path = os.path.join(base_dir, "peaks_table.tex")
with open(table_path, "w") as f:
    f.write(table_text)

print(f"Plot saved as '{plot_path}'")
print(f"LaTeX table saved as '{table_path}'")
