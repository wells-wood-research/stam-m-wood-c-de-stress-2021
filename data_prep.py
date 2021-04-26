# This script reads in csv files of DE-STRESS output data, combines them and removes fields that
# are not needed for analysis

# 0. Loading the releavnt packages-----------------------------------

# Loading relevant packages
import os
import numpy as np
import pandas as pd

# 1. Loading csv files and appending them all together---------------

# Defining the data path
DESTRESS_OUTPUT_PATH = os.getenv("DESTRESS_OUTPUT_PATH")

# Listing all the files in this folder
file_name_list = [
    f
    for f in os.listdir(DESTRESS_OUTPUT_PATH)
    if os.path.isfile(os.path.join(DESTRESS_OUTPUT_PATH, f))
]

# Looping through the de-stress output files
# and appending them all together
for i, file_name in enumerate(file_name_list):
    # Reading in the first file and naming it
    # combined_data
    if i == 0:
        # Reading the csv file using pandas
        combined_data = pd.read_csv(DESTRESS_OUTPUT_PATH + file_name)

    # For the other files we read them in and then append
    # them to combined_data
    else:
        # Reading the csv file using pandas
        data = pd.read_csv(DESTRESS_OUTPUT_PATH + file_name)

        # Appending the data together
        combined_data = combined_data.append(data, ignore_index=True)

# 2. Creating some new columns and removing columns
# that are not needed-------------------------------------------------------

# Dropping the file_name column
combined_data.drop(["file name", "tags"], axis=1, inplace=True)

# Creating a column to indicate decoy or native
combined_data["decoy or native"] = np.where(
    combined_data["design name"].str.contains("decoy"),
    "decoy",
    np.where(combined_data["design name"].str.contains("native"), "native", ""),
)

# Creating a column to indicate the pdb id
combined_data["pdb id"] = combined_data["design name"].str.split(pat="_").str[0]

# print(combined_data["number of residues"].unique())
# print(combined_data["mass (da)"].unique())

# 3. Saving data--------------------------------------------------------------

# Outputting the pandas data frame as a csv file
combined_data.to_csv(DESTRESS_OUTPUT_PATH + "combined_data.csv", index=False)
