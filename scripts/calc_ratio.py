import pandas as pd

# Load the count file (adjust the filename if needed)
df = pd.read_csv("Phase_out.pe.count", sep="\t", header=None, names=["inv_id", "direction", "count"])

# Pivot to reshape: one row per inv_id, columns for F and R
pivot = df.pivot_table(index="inv_id", columns="direction", values="count", aggfunc="sum", fill_value=0)

# Calculate ratio: R / (F + R)
pivot["R_to_total"] = pivot["R"] / (pivot["F"] + pivot["R"])

# Output the result
print(pivot[["R_to_total"]])

# Optional: Save to a new file
pivot[["R_to_total"]].to_csv("ratios.tsv", sep="\t")
