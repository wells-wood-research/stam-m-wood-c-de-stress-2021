# DE-STRESS: A user-friendly web application for the evaluation of protein designs - Decoy Analysis

## Overview

This repository contains the python code and the output for the decoy structure analysis included in 
the DE-STRESS paper. This decoy analysis gives one example of how the metrics in DE-STRESS could be 
used. A set of 10 experimentally-determined structures along with 20 folding decoys were randomly 
sampled from the 3DRobot_set. Using the DE-STRESS web server, we generated and exported metrics for 
all of these structures. After this, the metrics were normalised with min-max scaling, so that the 
features ranged between 0 and 1. Finally, principal component analysis (pca) was performed on the 
normalised data set. 

By plotting the first two principal components, the experimentally-determined structures and 
their decoys formed clusters. Furthermore, the experimentally-determined structures were close to, 
but distinct from, the main cluster, indicating that the metrics included in DE-STRESS could be used
to automatically identify high-quality structure models using machine learning. 

## Data

### Data Set used for Decoy Analysis

The data set used for the decoy analysis was created using the 3DRobot_set generated by 
3DRobot [Deng et al., 2016](https://doi.org/10.1093/bioinformatics/btv601). A set of 10 
experimentally-determined structures along with 20 folding decoys were randomly sampled from 
the 3DRobot_set. After some pre-processing, the DE-STRESS webserver was used to generate metrics 
for this data and the final data set is saved as `de-stress_output/combined_data.csv`.


### Creating the Data Set

In order to recreate the full analysis, the 3DRobot_set can be downloaded from 
[here](https://zhanglab.dcmb.med.umich.edu/3DRobot/decoys/). The `3DRobot_set.tar.bz2` folder
was extracted in the `data/` folder in this repo. After all the subfolders were extracted,
the `fixing_pdb_files.py` script was ran. This script is needed as the decoy pdb files 
generated by 3DRobot have the element column missing. Also, some of the pdb files contained 
hydrogen atoms which are not accepted by a lot of the metrics ran in DE-STRESS. Both of these 
issues needed to be addressed before running these pdb files through the DE-STRESS web server. 
**(Note: This script takes quite a bit of time to run)**.

After the pdb files had been fixed, the `sampling_pdb_files.py` script was ran to sample 10 
experimentally-determined structures and 20 decoys were randomly sampled for each one. 
The decoy structures in the 3DRobot set have RMSD ranging from 0 to 12 Å, so the random sampling 
was done to ensure a spread across this range. 

Once the structures had been sampled, they were 
ran through the DE-STRESS web server in 10 batches and the data was downloaded as `.csv` files. 
These `.csv` files are saved as;

1. `de-stress_output/1N8VA.csv`
2. `de-stress_output/1ZI8A.csv`
3. `de-stress_output/2HS1A.csv`
4. `de-stress_output/2YHSA.csv`
5. `de-stress_output/3CHBD.csv`
6. `de-stress_output/3D32A.csv`
7. `de-stress_output/3MMHA.csv`
8. `de-stress_output/3NJNA.csv`
9. `de-stress_output/3WCQA.csv`
10. `de-stress_output/3WDCA.csv`.
                                 
Finally the `data_prep.py` script was used to combine all the `.csv` files together, to drop 
some columns that were not needed and to create the `decoy or native` and `pdb id` columns.
This data set was then saved as `de-stress_output/combined_data.csv`. 

## Principal Component Analysis (PCA)

### Feature Selection and Scaling

Firstly, before performing PCA we removed the amino acid composition features such as `composition_ALA`. 
This is because these metrics are constant for the decoys and native structures of one `pdb_id` and only 
change across the different structures. In this analysis, we are interested in features that can distinguish 
between the decoys and native structures, so that's why these metrics were excluded. Other metrics that were
excluded were `design name` and `pdb_id` as these are categorical features, `number of residues` and `mass (da)` 
as similarly to the composition metrics, these don't distinguish between decoys and native structures, 
`evoef2 - interD total` and `rosetta - yhh_planarity` as these were constant across the data set, and 
obviously `decoy or native` as we are interested in seeing if the DE-STRESS features can separate decoys and
native structures out themselves.
 

After this, we scaled the remaining features so that the values were between 0 and 1. 
This is because PCA is extremely sensitive to the magnitude of the values and higher magnitude features 
(e.g energy function values) could skew the analysis. We chose to normalise the features with min-max scaling 
rather than using standardisation, as some of the features didn't appear to have a Gaussian distribution. 
An example of this is shown in the plot below.

[<img src=./analysis/hist_dfire2.png>]()

### Performing PCA

Firstly, we checked the variance explained by different numbers of PCA components. The chart below
shows that the first 2 principal components explain roughly 67% of the variance in the data set, and
roughly 90% of the variance is explained with 6 principal components.

[<img src=./analysis/var_explained.png>]()

After this, we plotted the first two principal components for all 10 of the 
experimentally-determined structures along with their decoys, in separate plots. 

[<img src=./analysis/pca_2dproj_subplots.png>]()

From this plot, we see that for each `pdb id`, the native structure consistently has the maximum value
for principal component 2 compared to the decoy structures for the same `pdb id`. Based on this observation, 
we plotted principal component 2 for each `pdb id` on the same chart, with the circles representing the 
decoy structures and the stars representing the native structures. 

[<img src=./analysis/strip_plot.png>]()

The plot above shows that the native structures have the maximum value for principal component 2 across all 
the different structures, and shows that the DE-STRESS metrics are able to consistently separate out the native 
structures from the decoys. This plot is included in the updated DE-STRESS paper.  

Furthermore, we investigated what features contributed to principal component 2 and the top 10 features, ordered
by relative contribution are shown below. 

1. `rosetta - hbond_lr_bb`
2. `aggrescan3d: total_value`
3. `aggrescan3d: avg_value`
4. `rosetta - rama_prepro`
5. `rosetta - hbond_sr_bb`
6. `rosetta - omega`
7. `aggrescan3d: max_value`
8. `rosetta - total`
9. `evoef2: intraR total`
10. `rosetta - pro_close`.























