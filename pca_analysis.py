# This script reads in the csv file of the combined data from the DE-STRESS output,
# normalises the data with min max scaling and reduces the dimensionality with PCA.

# 0. Loading the relevant packages-------------------------------------------

# Loading relevant packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn import decomposition
from scipy.stats import multivariate_normal
from helper_functions import *
from scipy import stats


# 1. Reading in the data and preprocessing-----------------------------------

# Setting the file paths
destress_output = "de-stress_output/"
analysis_output = "analysis/"

# Loading in the combined_data
combined_data = pd.read_csv(destress_output + "combined_data.csv")

# Extracting the labels
decoy_or_native = combined_data["decoy or native"]
pdb_id = combined_data["pdb id"].str.slice(start=0, stop=4)
structure_group = combined_data["structure group"].str.slice(start=0, stop=4)
# extra_native_flag = combined_data["extra native flag"]

# Dividing energy values by number of residues
combined_data.loc[
    :,
    [
        "hydrophobic fitness",
        "budeff: total",
        "budeff: steric",
        "budeff: desolvation",
        "budeff: charge",
        "evoef2: total",
        "evoef2: ref total",
        "evoef2: intraR total",
        "evoef2: interS total",
        "evoef2 - interD total",
        "dfire2 - total",
        "rosetta - total",
        "rosetta - fa_atr",
        "rosetta - fa_rep",
        "rosetta - fa_intra_rep",
        "rosetta - fa_elec",
        "rosetta - fa_sol",
        "rosetta - lk_ball_wtd",
        "rosetta - fa_intra_sol_xover4",
        "rosetta - hbond_lr_bb",
        "rosetta - hbond_sr_bb",
        "rosetta - hbond_bb_sc",
        "rosetta - hbond_sc",
        "rosetta - dslf_fa13",
        "rosetta - rama_prepro",
        "rosetta - p_aa_pp",
        "rosetta - fa_dun",
        "rosetta - omega",
        "rosetta - pro_close",
        "rosetta - yhh_planarity",
    ],
] = combined_data.loc[
    :,
    [
        "hydrophobic fitness",
        "budeff: total",
        "budeff: steric",
        "budeff: desolvation",
        "budeff: charge",
        "evoef2: total",
        "evoef2: ref total",
        "evoef2: intraR total",
        "evoef2: interS total",
        "evoef2 - interD total",
        "dfire2 - total",
        "rosetta - total",
        "rosetta - fa_atr",
        "rosetta - fa_rep",
        "rosetta - fa_intra_rep",
        "rosetta - fa_elec",
        "rosetta - fa_sol",
        "rosetta - lk_ball_wtd",
        "rosetta - fa_intra_sol_xover4",
        "rosetta - hbond_lr_bb",
        "rosetta - hbond_sr_bb",
        "rosetta - hbond_bb_sc",
        "rosetta - hbond_sc",
        "rosetta - dslf_fa13",
        "rosetta - rama_prepro",
        "rosetta - p_aa_pp",
        "rosetta - fa_dun",
        "rosetta - omega",
        "rosetta - pro_close",
        "rosetta - yhh_planarity",
    ],
].div(
    combined_data["number of residues"], axis=0
)


# Dropping some columns
combined_data.drop(
    [
        "design name",
        "number of residues",
        "mass (da)",
        "decoy or native",
        "evoef2 - interD total",
        "rosetta - yhh_planarity",
        "aggrescan3d: total_value",
        "isoelectric point (pH)",
        "pdb id",
        "structure group",
        "evoef2: ref total",
        "rosetta - dslf_fa13",
    ],
    axis=1,
    inplace=True,
)


# Dropping composition metrics as these only change across different pdb ids
# and not natives vs decoy structures. This is important as we are interested
# in looking at metrics that separate out decoys vs native structures. If these
# values are included they dominate the variance in the data set as each structure
# has a different amino acid composition.
combined_data = combined_data.loc[
    :, ~combined_data.columns.str.startswith("composition")
]

# Saving a list of the features used
print(combined_data.columns)


# 2. Histograms of the features------------------------------------------------------------------------
combined_data["rosetta - pro_close"].plot.hist(
    grid=True, bins=50, rwidth=0.9, color="#87ceeb"
)
plt.title("Histogram of the rosetta - pro_close metric values")
plt.xlabel("rosetta - pro_close")
plt.ylabel("counts")
plt.savefig(
    analysis_output + "hist_rosetta_pro_close.png",
    bbox_inches="tight",
)
plt.close()


# 3. Performing PCA and plotting components
# against variance with two different scaling methods---------------------------------------------------
n_components = [2, 3, 4, 5, 6, 7, 8, 9, 10]
scaling_method = ["minmax"]

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
var_plot = sns.lineplot(
    x="n_components",
    y="var_explained",
    data=var_explained_df,
    hue="scaling_method",
)
plt.title("""Variance explained by number of pca components and scaling method.""")
plt.xlabel("Number of components")
plt.legend(loc="upper right", ncol=2, handletextpad=0.1)
plt.savefig(analysis_output + "var_explained.png")
plt.close()


fig, (ax1, ax2) = plt.subplots(2, 1)
fig.set_size_inches(6, 7)

combined_data["dfire2 - total"].plot.hist(
    grid=True, bins=50, rwidth=0.9, color="#87ceeb", ax=ax1
)
ax1.set_xlabel("DFIRE2 - Total")
ax1.set_ylabel("Counts")

var_plot = sns.lineplot(
    x="n_components",
    y="var_explained",
    data=var_explained_df,
    hue="scaling_method",
    ax=ax2,
    legend=False,
)
ax2.set_xlabel("Number of Components")
ax2.set_ylabel("Variance Explained")
plt.savefig(
    analysis_output + "hist_dfire2_var_explained.png",
    dpi=600,
)
plt.close()


# 4. Performing PCA with 2 components and plotting against different labels

# Normalising the data with min max transform as not all the
# features will have a gaussian distribution
scaled_data = MinMaxScaler(feature_range=(0, 1)).fit_transform(combined_data)
scaled_data_df = pd.DataFrame(scaled_data, columns=combined_data.columns)

# Checking normality of the features
for col in scaled_data_df.columns:
    scaled_data_df[col].plot.hist(grid=True, bins=50, rwidth=0.9, color="#87ceeb")
    plt.title(col)
    plt.xlabel("col")
    plt.ylabel("counts")
    plt.savefig(
        analysis_output + "hist_transf_" + col + ".png",
        bbox_inches="tight",
    )
    plt.close()

# Outputting the pandas data frame as a csv file
scaled_data_df.to_csv(analysis_output + "scaled_data.csv", index=False)

# Performing PCA
pca = decomposition.PCA(n_components=2)
pca.fit(scaled_data_df)

# Saving contributions of the features to the principal components
feat_contr_to_cmpts = pd.DataFrame(
    np.round(abs(pca.components_), 4), columns=combined_data.columns
)
feat_contr_to_cmpts.to_csv(analysis_output + "feat_contr_to_cmpts.csv", index=True)

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
transformed_data = pd.concat(
    [transformed_data, decoy_or_native, pdb_id, structure_group],
    axis=1,
)

# 5. Fitting 2d Gaussians to each of the sets of decoy structures by pdb id----------------------------

# 6. Various different plots for paper and report---------------------------------------

# Scatter plot of all PCA component 1 and PCA component 2 for all structures
plot = sns.scatterplot(
    x="pca_dim0",
    y="pca_dim1",
    data=transformed_data,
    hue="structure group",
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
    markers=["o", "*", "s"],
    alpha=0.6,
    edgecolor="black",
    sizes=[50, 150, 50],
)
plt.xlabel("Principal Component 1")
plt.ylabel("Principal Component 2")
h, l = plot.get_legend_handles_labels()
h[0] = "PDB ID"
plt.legend(
    # h[0:11],
    # l[0:11],
    bbox_to_anchor=(1.05, 1),
    loc="upper left",
)
# plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.savefig(analysis_output + "pca_2dproj.png", bbox_inches="tight", dpi=600)
plt.savefig(analysis_output + "pca_2dproj.svg", bbox_inches="tight", dpi=600)
plt.close()


# Individual scatter plots for each structure
sns.set(font_scale=1.5)
# Creating facet grid
g = sns.FacetGrid(
    data=transformed_data,
    col="structure group",
    col_wrap=3,
    hue="decoy or native",
    height=5,
    aspect=1.5,
    sharex=True,
    sharey=True,
)
# Scatter plot split out by decoy or native
g.map(
    sns.scatterplot,
    "pca_dim0",
    "pca_dim1",
    s=200,
    legend=False,
)

# Adding labels
axes = g.fig.axes
for ax in axes:
    ax.set_xlabel(
        "Principal Component 1",
    )
    ax.set_ylabel("Principal Component 2")

plt.savefig(analysis_output + "pca_2dproj_subplots.svg", bbox_inches="tight", dpi=600)
plt.savefig(analysis_output + "pca_2dproj_subplots.png", bbox_inches="tight", dpi=600)
plt.close()

# Strip plot of PCA component 2 for all structures
# split out by decoy or native structure
sns.set(style="white")
strip_plot1 = sns.stripplot(
    x="structure group",
    y="pca_dim1",
    data=transformed_data[transformed_data["decoy or native"] == "decoy"],
    hue=decoy_or_native,
    edgecolor="black",
    alpha=0.7,
    jitter=True,
    linewidth=1,
    marker="o",
    s=6,
)
strip_plot1.legend_.remove()
strip_plot2 = sns.stripplot(
    x="structure group",
    y="pca_dim1",
    data=transformed_data[transformed_data["decoy or native"] == "native"],
    hue=decoy_or_native,
    edgecolor="black",
    alpha=0.7,
    jitter=True,
    linewidth=1,
    marker="*",
    s=10,
)
strip_plot2.legend_.remove()
plt.xlabel("PDB ID")
plt.ylabel("Principal Component 2")
plt.savefig(analysis_output + "strip_plot.svg", bbox_inches="tight", dpi=600)
plt.savefig(analysis_output + "strip_plot.png", bbox_inches="tight", dpi=600)
plt.close()


fig, (ax1, ax2) = plt.subplots(2, 1)
fig.set_size_inches(6, 8)

plot = sns.scatterplot(
    x="pca_dim0",
    y="pca_dim1",
    data=transformed_data,
    hue="structure group",
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
    markers=["o", "*", "s"],
    alpha=0.7,
    edgecolor="black",
    sizes=[50, 150, 50],
    ax=ax1,
)
ax1.set_xlabel("Principal Component 1")
ax1.set_ylabel("Principal Component 2")
h, l = plot.get_legend_handles_labels()
ax1.legend(h[1:11], l[1:11], bbox_to_anchor=(1.05, 1), loc="upper left")


# Strip plot of PCA component 2 for all structures
# split out by decoy or native structure
strip_plot1 = sns.stripplot(
    x="structure group",
    y="pca_dim1",
    data=transformed_data[transformed_data["decoy or native"] == "decoy"],
    hue=decoy_or_native,
    edgecolor="black",
    alpha=0.7,
    jitter=True,
    linewidth=1,
    marker="o",
    s=6,
    ax=ax2,
)
strip_plot1.legend_.remove()
strip_plot2 = sns.stripplot(
    x="structure group",
    y="pca_dim1",
    data=transformed_data[transformed_data["decoy or native"] == "native"],
    hue=decoy_or_native,
    edgecolor="black",
    alpha=0.7,
    jitter=True,
    linewidth=1,
    marker="*",
    s=10,
    ax=ax2,
)
strip_plot2.legend_.remove()
ax2.set_xlabel("PDB ID")
ax2.set_ylabel("Principal Component 2")
plt.savefig(analysis_output + "strip_plot_combined.png", bbox_inches="tight", dpi=600)
plt.close()
