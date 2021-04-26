# The 3DRobot Decoy set is used for the DE-STRESS decoy analysis. However, all of the PDB files for the decoy files
# are not valid pdb files. This script fixes the pdb files so that we can use them in DE-STRESS.

# Loading packages
import ampal
import os

# Setting the DECOY_DATA_FOLDER_PATH from the .env file
DECOY_DATA_FOLDER_PATH = os.getenv("DECOY_DATA_FOLDER_PATH")

# Listing all the sub dirs in the 3DRobot Decoy Data
subdirs = [x[0] for x in os.walk(DECOY_DATA_FOLDER_PATH)]

# Removing the DECOY_DATA_FOLDER_PATH from this list
subdirs[:] = [d for d in subdirs if d not in [DECOY_DATA_FOLDER_PATH]]
print(subdirs)

# Printing how many directories the script will run for
print("Fixing the pdb files for " + str(len(subdirs)) + " directories")

# Looping through all the sub directories
for subdir in subdirs:

    # Extracting the name of the subdir
    subdir_name = os.path.basename(subdir)
    print("Subdir: " + subdir_name)

    # Loading the native structure
    native_structure = ampal.load_pdb(subdir + "/" + "native.pdb")

    # Renaming the native pdb file to include the pdb id and saving it
    with open(subdir + "/" + subdir_name + "_native.pdb", "w") as outf:
        outf.write(native_structure.pdb)

    # Extracting all the file names of the decoy structures
    decoy_structures = [
        f
        for f in os.listdir(subdir)
        if os.path.isfile(os.path.join(subdir, f)) and "decoy" in f
    ]

    # Looping through the decoy structure files
    for decoy_structure in decoy_structures:

        # Load the structure
        decoy_pdb = ampal.load_pdb(subdir + "/" + decoy_structure)

        # Add the missing atom label e.g. res_label = "CA", it assigns element as "C"
        for atom in decoy_pdb.get_atoms():
            atom.element = atom.res_label[0]

        # Delete the "H" atoms
        for residue in decoy_pdb.get_monomers():
            if "H" in residue.atoms:
                del residue.atoms["H"]

        # Write the new file
        with open(subdir + "/" + subdir_name + "_" + decoy_structure, "w") as outf:
            outf.write(decoy_pdb.pdb)
