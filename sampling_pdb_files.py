# This script copies over the native protein structure for each pdb id and
# randomly samples decoys for each native structure. These pdb files are saved
# in another folder called 3DRobot_subset.

# 0. Loading the relevant packages--------------------------------------------

import os
import numpy as np
import pandas as pd
import subprocess
import pathlib
import random
from random import sample

# 1. Setting different parameters----------------------------------------------

# Setting the file paths
data_3drobot_set = "data/3DRobot_set/"
data_3drobot_subset = "data/3DRobot_subset/"

# Setting the seed
random.seed(42)

# Setting the number of structures
num_structures = 10

# Setting the number of decoys
num_decoys = 40

# 2. Creating the subset data----------------------------------------------------

# Defining the subset of data
subdirs_subset = [
    "data/3DRobot_set/1N8VA",
    "data/3DRobot_set/1ZI8A",
    "data/3DRobot_set/2HS1A",
    "data/3DRobot_set/3CHBD",
    "data/3DRobot_set/3NJNA",
    "data/3DRobot_set/3WCQA",
    "data/3DRobot_set/3WDCA",
    "data/3DRobot_set/3LDCA",
    "data/3DRobot_set/2XODA",
]

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
    native_file_path = (
        data_3drobot_set + subdir_name + "/" + subdir_name + "_native.pdb"
    )

    # Creating the subfolder in the subset folder
    subprocess.run(["mkdir", data_3drobot_subset + subdir_name])

    # Copying this structure across to the subset folder
    subprocess.run(
        [
            "cp",
            native_file_path,
            data_3drobot_subset + subdir_name + "/" + subdir_name + "_native.pdb",
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
                data_3drobot_set + subdir_name + "/" + decoy_structure,
                data_3drobot_subset + subdir_name + "/" + decoy_structure,
            ]
        )
