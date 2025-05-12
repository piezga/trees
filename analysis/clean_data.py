import pandas as pd
import matplotlib.pyplot as plt
from variables import *

census = 4

file = wanang_path + wanang_template.format(census = census) 

# Load the CSV, drop empty rows, and remove missing/zero coordinates
df = pd.read_csv(file, skip_blank_lines=True)
df = df.dropna(how="all")  # Drop fully empty rows
df = df.dropna(subset=["px", "py"])  # Drop rows with missing coordinates
df = df[(df["px"] != 0) & (df["py"] != 0)]  # Remove zero coordinates
df = df.reset_index(drop=True)                # Reset index to avoid gaps

# Set reference tree (e.g., first entry)
ref_x, ref_y = df["px"].min(), df["py"].min()

# Calculate relative coordinates
df["x"] = df["px"] - ref_x
df["y"] = df["py"] - ref_y



output_df = df[["x", "y", "mnemonic"]].rename(columns={"mnemonic": "name"})

print(output_df)

output_df.to_csv(wanang_path + f"wanang{census}.csv", index=False)