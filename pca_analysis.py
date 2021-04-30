# This script reads in the csv file of the combined data from the DE-STRESS output,
# normalises the data with min max scaling and reduces the dimensionality with PCA.

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

# Setting the file paths
destress_output = "de-stress_output/"
analysis_output = "analysis/"

# Loading in the combined_data
combined_data = pd.read_csv(destress_output + "combined_data.csv")

# Extracting the labels
decoy_or_native = combined_data["decoy or native"]
pdb_id = combined_data["pdb id"].str.slice(start=0, stop=4)

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

# 2. Histograms of the features------------------------------------------------------------------------

combined_data["dfire2 - total"].plot.hist(
    grid=True, bins=50, rwidth=0.9, color="#87ceeb"
)
plt.title("Histogram of the dfire2 - total metric values")
plt.xlabel("dfire2 - total")
plt.ylabel("counts")
plt.savefig(
    analysis_output + "hist_dfire2.png",
    bbox_inches="tight",
)
plt.close()


# 2. Performing PCA and plotting components
# against variance with two different scaling methods---------------------------------------------------
n_components = [2, 3, 4, 5, 6, 7, 8, 9, 10]
scaling_method = ["stand", "minmax"]

# Creating a data frame to collect results
var_explained_df = pd.DataFrame(
    columns=["n_components", "scaling_method", "var_explained"]
)

# Looping through different scaling methods and components
for i in n_components:
    for j in scaling_method:

        # Performing the scaling method
        if j == "stand":
            scaled_data = StandardScaler().fit_transform(combined_data)

        elif j == "minmax":
            scaled_data = MinMaxScaler(feature_range=(0, 1)).fit_transform(
                combined_data
            )

        # Scaling the data
        scaled_data_df = pd.DataFrame(scaled_data, columns=combined_data.columns)

        # Performing PCA with the specified components
        pca = decomposition.PCA(n_components=i)
        pca.fit(scaled_data_df)

        # Calculating the variance explained
        var_explained = np.sum(pca.explained_variance_ratio_)

        # Appending to the data frame
        var_explained_df = var_explained_df.append(
            {"n_components": i, "scaling_method": j, "var_explained": var_explained},
            ignore_index=True,
        )

# Saving as a csv file
var_explained_df.to_csv(analysis_output + "var_explained.csv", index=False)

# Plotting the data and saving
plot = sns.lineplot(
    x="n_components",
    y="var_explained",
    data=var_explained_df,
    hue="scaling_method",
)
plt.title("""Variance explained by number of pca components and scaling method.""")
plt.legend(loc="upper right", ncol=2, handletextpad=0.1)
plt.savefig(analysis_output + "var_explained.png")
plt.close()


# 3. Performing PCA with 2 components and plotting against different labels

# Normalising the data with min max transform as not all the
# features will have a gaussian distribution
scaled_data = MinMaxScaler(feature_range=(0, 1)).fit_transform(combined_data)
scaled_data_df = pd.DataFrame(scaled_data, columns=combined_data.columns)

# Outputting the pandas data frame as a csv file
scaled_data_df.to_csv(analysis_output + "scaled_data.csv", index=False)

# Performing PCA
pca = decomposition.PCA(n_components=2)
pca.fit(scaled_data_df)

# Saving contributions of the features to the principal components
feat_contr_to_cmpts = pd.DataFrame(abs(pca.components_), columns=combined_data.columns)
feat_contr_to_cmpts.to_csv(analysis_output + "feat_contr_to_cmpts.csv", index=False)

# Selecting the 10 largest contributers to pca component 1
comp_1_contr = feat_contr_to_cmpts.iloc[0].nlargest(10, keep="first")
comp_1_contr.to_csv(analysis_output + "comp_1_contr.csv", index=True)

# Selecting the 10 largest contributers to pca component 2
comp_2_contr = feat_contr_to_cmpts.iloc[1].nlargest(10, keep="first")
comp_2_contr.to_csv(analysis_output + "comp_2_contr.csv", index=True)

# Transforming the data
transformed_data = pca.transform(scaled_data_df)

# Converting to data frame and renaming columns
transformed_data = pd.DataFrame(transformed_data).rename(
    columns={0: "pca_dim0", 1: "pca_dim1"}
)

# Adding the labels back
transformed_data = pd.concat([transformed_data, decoy_or_native, pdb_id], axis=1)

# Plotting and saving the figure
plot = sns.scatterplot(
    x="pca_dim0",
    y="pca_dim1",
    data=transformed_data,
    hue="pdb id",
    hue_order=[
        "1N8V",
        "1ZI8",
        "2HS1",
        "2YH5",
        "3CHB",
        "3D32",
        "3MMH",
        "3NJN",
        "3WCQ",
        "3WDC",
    ],
    style="decoy or native",
    size="decoy or native",
    markers=["o", "*"],
    sizes=[50, 150],
)
plt.xlabel("PCA Component 1")
plt.ylabel("PCA Component 2")
h, l = plot.get_legend_handles_labels()
plt.legend(h[1:11], l[1:11], loc="best", ncol=3, handletextpad=0.1)
plt.savefig(
    analysis_output + "pca_2dproj.svg",
    bbox_inches="tight",
)
