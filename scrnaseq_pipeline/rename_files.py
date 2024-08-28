import os
import glob
import pandas as pd

# Define the data directory
data_dir = "/Users/lzy/Desktop/Final/data"

# Get all files in the directory (excluding those already in TSV format)
files = sorted(glob.glob(os.path.join(data_dir, "*")))
tsv_files = [f for f in files if f.endswith('.tsv')]

# Iterate over each file and convert to TSV format if necessary
for file_path in files:
    if file_path in tsv_files:
        print(f"File {file_path} is already in TSV format. Skipping conversion.")
        continue

    # Read the file
    try:
        # Assume the file is in CSV format or another format separated by commas
        df = pd.read_csv(file_path, sep=None, engine='python')
        
        # Get the base name of the file (without the extension)
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        
        # Generate the new TSV file path
        new_name = os.path.join(data_dir, f"{base_name}.tsv")
        
        # Save the data as a TSV file
        df.to_csv(new_name, sep='\t', index=False)
        print(f"Converted {file_path} to {new_name}.")
        
        # Delete the original file
        os.remove(file_path)
        print(f"Deleted original file {file_path}.")
    except Exception as e:
        print(f"Error processing {file_path}: {e}")

print("All files have been successfully converted to TSV format and original files have been deleted.")
