# This script reads in csv files of DE-STRESS output data,
# combines them and removes fields that are not needed for analysis.

# 0. Loading the releavnt packages-----------------------------------

# Loading relevant packages
import os
import numpy as np
import pandas as pd

# 1. Loading csv files and appending them all together---------------

# Setting the file paths
destress_output = "de-stress_output/"

# Listing all the files in this folder
file_name_list = [
    f
    for f in os.listdir(destress_output)
    if os.path.isfile(os.path.join(destress_output, f))
    and f.endswith(".csv")
    and f != "combined_data.csv"
]

# Looping through the de-stress output files
# and appending them all together
for i, file_name in enumerate(file_name_list):
    # Reading in the first file and naming it
    # combined_data
    if i == 0:
        # Reading the csv file using pandas
        combined_data = pd.read_csv(destress_output + file_name)

    # For the other files we read them in and then append
    # them to combined_data
    else:
        # Reading the csv file using pandas
        data = pd.read_csv(destress_output + file_name)

        # Appending the data together
        combined_data = combined_data.append(data, ignore_index=True)

# 2. Creating some new columns and removing columns
# that are not needed-------------------------------------------------------

# Dropping the file_name column
combined_data.drop(["file name", "tags"], axis=1, inplace=True)

# Creating a column to indicate the pdb id
combined_data["pdb id"] = combined_data["design name"].str.split(pat="_").str[0]

# Creating a column to indicate decoy or native
combined_data["decoy or native"] = np.where(
    combined_data["design name"].str.contains("decoy"),
    "decoy",
    np.where(
        combined_data["pdb id"].isin(
            [
                "1ZIC",
                "1ZIX",
                "1ZI9",
                "2HS2",
                "3S53",
                "1PZKD",
                "1PZJD",
                "3NJHA",
                "3NJMA",
                "3WDEA",
                "3WDDA",
                "1N8UA",
                "3AB5A",
                "2X2P",
                "3OUS",
                "3R65",
                "3LDD",
            ]
        ),
        "additional native",
        np.where(combined_data["design name"].str.contains("native"), "native", ""),
    ),
)


# Creating a column to indicate structure group
combined_data["structure group"] = np.where(
    combined_data["pdb id"].isin(["1ZIC", "1ZIX", "1ZI9"]),
    "1ZI8A",
    np.where(
        combined_data["pdb id"].isin(["2HS2", "3S53"]),
        "2HS1A",
        np.where(
            combined_data["pdb id"].isin(["1PZKD", "1PZJD"]),
            "3CHBD",
            np.where(
                combined_data["pdb id"].isin(["3NJHA", "3NJMA"]),
                "3NJN",
                np.where(
                    combined_data["pdb id"].isin(["3WDEA", "3WDDA"]),
                    "3WDC",
                    np.where(
                        combined_data["pdb id"].isin(["1N8UA"]),
                        "1N8V",
                        np.where(
                            combined_data["pdb id"].isin(["3AB5A"]),
                            "3WCQA",
                            np.where(
                                combined_data["pdb id"].isin(["2X2P"]),
                                "2XODA",
                                np.where(
                                    combined_data["pdb id"].isin(
                                        ["3OUS", "3R65", "3LDD"]
                                    ),
                                    "3LDCA",
                                    combined_data["design name"]
                                    .str.split(pat="_")
                                    .str[0],
                                ),
                            ),
                        ),
                    ),
                ),
            ),
        ),
    ),
)


# 3. Saving data--------------------------------------------------------------

# Outputting the pandas data frame as a csv file
combined_data.to_csv(destress_output + "combined_data.csv", index=False)
