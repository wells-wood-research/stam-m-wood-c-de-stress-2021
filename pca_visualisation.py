# This script reads in the csv file of the combined data from the DE-STRESS output,
# reduces the dimensionality with PCA and plots this data

# 0. Loading the relevant packages-------------------------------------------

# Loading relevant packages
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn import decomposition


# 1. Reading in the data and preprocessing-----------------------------------

# Defining the data path
DESTRESS_OUTPUT_PATH = os.getenv("DESTRESS_OUTPUT_PATH")

# Loading in the combined_data
combined_data = pd.read_csv(DESTRESS_OUTPUT_PATH + "combined_data.csv")

# Extracting the labels
decoy_or_native = combined_data["decoy or native"]
pdb_id = combined_data["pdb id"]

# Dropping some columns
combined_data.drop(
    [
        "design name",
        "number of residues",
        "mass (da)",
        "decoy or native",
        "evoef2 - interD total",
        "rosetta - yhh_planarity",
        "pdb id",
    ],
    axis=1,
    inplace=True,
)

# Standardising the data
# scaled_data = StandardScaler().fit_transform(combined_data)
# scaled_data_df = pd.DataFrame(scaled_data, columns=combined_data.columns)

# Normalising the data with min max transform
scaled_data = MinMaxScaler(feature_range=(0, 1)).fit_transform(combined_data)
scaled_data_df = pd.DataFrame(scaled_data, columns=combined_data.columns)

# Outputting the pandas data frame as a csv file
scaled_data_df.to_csv(DESTRESS_OUTPUT_PATH + "scaled_data.csv", index=False)

# 2. Performing PCA and plotting---------------------------------------------------
pca = decomposition.PCA(n_components=6)
pca.fit(scaled_data_df)
print(np.sum(pca.explained_variance_ratio_))
transformed_data = pca.transform(scaled_data_df)

# # Converting to data frame and renaming columns
# transformed_data = pd.DataFrame(transformed_data).rename(
#     columns={0: "pca_dim0", 1: "pca_dim1"}
# )

# # Adding back on labels
# transformed_data = pd.concat([transformed_data, decoy_or_native, pdb_id], axis=1)

# # Plotting
# sns.scatterplot(
#     x="pca_dim0",
#     y="pca_dim1",
#     data=transformed_data,
#     hue="decoy or native",
# )
# plt.show()