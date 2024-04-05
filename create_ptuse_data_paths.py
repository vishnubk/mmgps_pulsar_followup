import pandas as pd
import sys

import argparse
import pandas as pd

# Set up argument parser
parser = argparse.ArgumentParser(description="Process grouped files and generate ptuse paths.")
parser.add_argument("-f", "--filename", type=str, required=True, help="Path to the input file (grouped_files.txt).")

# Parse arguments
args = parser.parse_args()

# Read the arguments
file_path = args.filename
ptuse_path1 = "/beegfs/DATA/MeerTIME/SCI-20200703-MK-01/search"
ptuse_path2 = "/beegfs/DATA/PTUSE/SCI-20200703-MK-01/search"
ptuse_path3 = "/beegfs/DATA/PTUSE/SCI-20200703-MK-02/search"

# Read the file, assuming tab-separated values and skipping the header
df = pd.read_csv(file_path, sep='\t', header=0)

# Deduplicate based on 'utc_start', keeping the first occurrence
deduplicated_df = df.drop_duplicates(subset=['utc_start'])

ptuse_data_root = []
# For each deduplicated entry, generate ptuse paths
for index, row in deduplicated_df.iterrows():
    date_with_hour = row['utc_start'].split(' ')[0] + '-' + row['utc_start'].split(' ')[1].split(':')[0]
    psrfits_path1 = ptuse_path1 + '/' + date_with_hour + '*/' + row['target'] + '/**/*.sf' 
    psrfits_path2 = ptuse_path2 + '/**/' + date_with_hour + '*/' + row['target'] + '/**/*.sf'  
    psrfits_path3 = ptuse_path3 + '/**/' + date_with_hour + '*/' + row['target'] + '/**/*.sf'
    ptuse_data_root.append([psrfits_path1, psrfits_path2, psrfits_path3, row['target'], "ptuse", date_with_hour])
    

# Print or save the output
df = pd.DataFrame(ptuse_data_root, columns=['PTUSE_PATH1', 'PTUSE_PATH2', 'PTUSE_PATH3', 'target', 'beam', 'utc' ])

df.to_csv('ptuse_data_paths.txt', index=False, header=True, sep='\t')




