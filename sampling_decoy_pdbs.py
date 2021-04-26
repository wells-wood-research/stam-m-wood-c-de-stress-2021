# This script copies over the native protein structure for each pdb id and
# samples decoys for each one

# 0. Loading the relevant packages-----------------------------------

import os
import numpy as np
import pandas as pd
import subprocess
import random
from random import sample

# 1. Setting the file paths and parameters-------------------------------------------

# Setting the DECOY_DATA_FOLDER_PATH and DECOY_SUBSET_FOLDER_PATH from the .env file
DECOY_DATA_FOLDER_PATH = os.getenv("DECOY_DATA_FOLDER_PATH")
DECOY_SUBSET_FOLDER_PATH = os.getenv("DECOY_SUBSET_FOLDER_PATH")

# Setting the seed
random.seed(42)

# Setting the number of structures
num_structures = 10

# Setting the number of decoys
num_decoys = 20

# 2. Creating the subset data-------------------------------------------------------------------

# Listing all the sub dirs (structures) in the 3DRobot Decoy Data
subdirs = [x[0] for x in os.walk(DECOY_DATA_FOLDER_PATH)]

# Removing the DECOY_DATA_FOLDER_PATH from this list
subdirs[:] = [d for d in subdirs if d not in [DECOY_DATA_FOLDER_PATH]]

# Sampling a subset of these structures
subdirs_subset = sample(subdirs, num_structures)

# Printing how many structures the script will run for
print(
    "Sampling "
    + str(num_decoys)
    + " decoys for "
    + str(len(subdirs_subset))
    + " structures"
)

# Looping through all the sub directories
for subdir in subdirs_subset:

    # Extracting the name of the subdir
    subdir_name = os.path.basename(subdir)

    # Defining the native structure file path
    native_file_path = DECOY_DATA_FOLDER_PATH + subdir_name + "/native.pdb"

    # Creating the subfolder in the subset folder
    subprocess.run(["mkdir", DECOY_SUBSET_FOLDER_PATH + subdir_name])

    # Copying this structure across to the subset folder
    subprocess.run(
        [
            "cp",
            native_file_path,
            DECOY_SUBSET_FOLDER_PATH + subdir_name + "/" + subdir_name + "_native.pdb",
        ]
    )

    # Extracting all the file names of the decoy structures
    decoy_structures = [
        f
        for f in os.listdir(subdir)
        if os.path.isfile(os.path.join(subdir, f)) and "decoy" in f and subdir_name in f
    ]

    # Sampling decoy files
    sample_decoys = sample(decoy_structures, num_decoys)

    # Copying over these decoy files
    for decoy_structure in sample_decoys:

        # Copying this structure across to the subset folder
        subprocess.run(
            [
                "cp",
                DECOY_DATA_FOLDER_PATH + subdir_name + "/" + decoy_structure,
                DECOY_SUBSET_FOLDER_PATH + subdir_name + "/" + decoy_structure,
            ]
        )
