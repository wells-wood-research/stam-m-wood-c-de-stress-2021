# DE-STRESS: A user-friendly web application for the evaluation of protein designs - Decoy Analysis

## Overview

This repository contains the python code and the output for the decoy structure analysis included 
in the DE-STRESS paper. The decoy analysis gives one example of how the metrics in DE-STRESS could
be used. A set of 10 experimentally-determined structures along with 40 folding decoys were 
randomly sampled from the 3DRobot_set. Using the DE-STRESS web server, we generated and exported 
metrics for all of these structures. For some of the experimentally-determined structures, other 
crystallographic structures were available from the PDB, that were not part of the 3DRobot_set. 
These were included in the analysis, as we would expect these to have similar metric values to the 
experimentally-determined structure in the 3DRobot_set.  After this, the metrics were scaled with 
min-max scaling, so that the features ranged between 0 and 1. Finally, principal component analysis
(pca) was performed on the scaled data set. 

By plotting the first two principal components, the experimentally-determined structures were shown
to have the largest values for both principal components compared to their decoys. The extra 
crystallographic structures from the PDB were very close to the experimentally-determined structure
from the 3DRobot_set. This shows that these results are robust to different crystallographic structures
for the same protein. Also, the 2d PCA projection of the experimentally-determined structures and 
decoys are shown to be linerarly separable, which means it would be easy to train a classifier 
to predict the experimentally-determined structures from the decoys.

## Data

### Data Set used for Decoy Analysis

The data set used for the decoy analysis was created using the 3DRobot_set generated by 
3DRobot [Deng et al., 2016](https://doi.org/10.1093/bioinformatics/btv601). A set of 10 
experimentally-determined structures along with 40 folding decoys were randomly sampled from 
the 3DRobot_set. The extra crystallographic structures were downloaded from the PDB and added 
to the data set. After some pre-processing, the DE-STRESS webserver was used to generate metrics 
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
experimentally-determined structures and 40 decoys were randomly sampled for each one. 
The decoy structures in the 3DRobot set have RMSD ranging from 0 to 12 Å, so the random sampling 
was done to ensure a spread across this range. 

The extra crystallographic structures that were included in the analysis were found using the 
"Find similar assemblies" search functionality on the PDB, and structures with a similarity 
score greater than 70% were selected. Not all structures had similar structures that met this 
criteria. The table below shows pdb ids and chains for the structures that were selected from 
the PDB, and the experimentally determined structures from the 3DRobot_set. 

| Experimentally-determined structure <br />  pdb id and chain  |  Similar structures pdb ids <br /> and chains | 
|---|---|
| 1N8V A | 1N8U A |
| 1ZI8 A | 1ZIC A, 1ZIX A, 1ZI9 A |
| 2HS1 A | 2HS2 A, 3S53 A |
| 2YH5 A |  None |
| 3CHB D | 1PZK D, 1PZJ D |
| 3D32 A | None |
| 3MMH A | None  |
| 3NJN A | 3NJH A, 3NJM A |
| 3WCQ A | None |
| 3WDC A | 3WDE A, 3WDD A |

All of these structures, along with their decoys, were ran through the DE-STRESS web server in batches, and the data was 
downloaded as `.csv` files. These `.csv` files are saved as;

1. `de-stress_output/1N8VA_40decoys.csv`
2. `de-stress_output/1N8VA_extra_natives.csv`
3. `de-stress_output/1ZI8A_40decoys.csv`
4. `de-stress_output/1ZI8A_extra_natives.csv`
5. `de-stress_output/2HS1A_40decoys.csv`
6. `de-stress_output/2HS1A_extra_natives.csv`
7. `de-stress_output/2YHSA_40decoys.csv`
8. `de-stress_output/3CHBD_40decoys.csv`
9. `de-stress_output/3CHBD_extra_natives.csv`
10. `de-stress_output/3D32A_40decoys.csv`
11. `de-stress_output/3MMHA_40decoys.csv`
12. `de-stress_output/3NJNA_40decoys.csv`
13. `de-stress_output/3NJNA_extra_natives.csv`
14. `de-stress_output/3WCQA_40decoys.csv`
15. `de-stress_output/3WDCA_40decoys.csv`
15. `de-stress_output/3WDCA_extra_natives.csv`

Finally the `data_prep.py` script was used to combine all the `.csv` files together, to drop 
some columns that were not needed and to create the `decoy or native` and `pdb id` columns.
This data set was then saved as `de-stress_output/combined_data.csv`. 

## Principal Component Analysis (PCA)

### Feature Selection and Scaling

Firstly, before performing PCA, a number of data pre processing steps were peformed. As the energy
values are dependent on the size of the protein, we divided all of these values by the number
of residues in the structure. After this, we removed a number of variables from the data set.
Features such as `composition - ALA`, `mass (da)`, `number of residues` and `isoelectric point (pH)`
were constant for structures from the same pdb id. In this analysis, we are interested in features
that can distinguish between the decoys and experimentally-determined structures, so that's why 
these metrics were excluded. Other metrics that were excluded were `design name`, `pdb_id` and 
`structure group` as these are categorical features, `evoef2 - interD total` 
and `rosetta - yhh_planarity` as these were constant across the data set, and obviously 
`decoy or native` as we are interested in seeing if the DE-STRESS features can separate decoys and
experimentally-determined structures out themselves. Two other metrics `evoef2: ref total` and
`rosetta - dslf_fa13` were excluded because they were discrete variables and `aggrescan3d: total_value`
was excluded as we already have `aggrescan3d: avg_value` in the data set which is normalised for 
the size of the structure. The two sections below show the full list of included and excluded metrics.

Included Metrics

    'hydrophobic fitness', 'packing density', 'budeff: total',
    'budeff: steric', 'budeff: desolvation', 'budeff: charge',
    'evoef2: total', 'evoef2: intraR total', 'evoef2: interS total', 
    'dfire2 - total', 'rosetta - total','rosetta - fa_atr', 
    'rosetta - fa_rep', 'rosetta - fa_intra_rep','rosetta - fa_elec', 
    'rosetta - fa_sol', 'rosetta - lk_ball_wtd', 'rosetta - fa_intra_sol_xover4', 
    'rosetta - hbond_lr_bb', 'rosetta - hbond_sr_bb', 'rosetta - hbond_bb_sc', 
    'rosetta - hbond_sc', 'rosetta - rama_prepro', 'rosetta - p_aa_pp',
    'rosetta - fa_dun', 'rosetta - omega', 'rosetta - pro_close',
    'aggrescan3d: avg_value', 'aggrescan3d: min_value', 'aggrescan3d: max_value'

Excluded Metrics

    'design name', 'number of residues', 'mass (da)',
    'decoy or native', 'evoef2 - interD total', 'rosetta - yhh_planarity',
    'aggrescan3d: total_value', 'isoelectric point (pH)',
    'pdb id', 'structure group', 'evoef2: ref total',
    'rosetta - dslf_fa13', 'composition: ALA', 'composition: CYS',
    'composition: ASP',	'composition: GLU',	'composition: PHE',
    'composition: GLY',	'composition: HIS',	'composition: ILE',
    'composition: LYS',	'composition: LEU',	'composition: MET',	
    'composition: ASN',	'composition: PRO',	'composition: GLN',
    'composition: ARG',	'composition: SER',	'composition: THR',
    'composition: VAL',	'composition: TRP',	'composition: UNK'



After excluding these metrics, the remaining features were scaled so that the values were 
between 0 and 1. This is because PCA is extremely sensitive to the magnitude of the values 
and higher magnitude features (e.g energy function values) could skew the analysis. 
We chose to scale the features with min-max scaling rather than using standardisation, 
as some of the features didn't appear to have a Gaussian distribution. 
An example of this is shown in the plot below.

[<img src=./analysis/hist_transf_rosetta-hbond_sc.png>]()

### Performing PCA

Firstly, we checked the variance explained by different numbers of PCA components. The chart below
shows that the first 2 principal components explain roughly 64% of the variance in the data set, and
roughly 90% of the variance is explained with 8 principal components.

[<img src=./analysis/var_explained.png>]()

After this, we plotted the first two principal components for all 10 of the 
experimentally-determined structures, along with their decoys and the other 
crystallgraphic structures.

[<img src=./analysis/pca_2dproj_subplots.png>]()

From this plot, we see that for each `pdb id`, the experimentally determined structures, shown as
stars, consistently have the maximum values for principal components 1 and 2. It is clear from this
plot that they are distinct from the decoy structures shown as circles. In addition to this, the
other crystallographic structures, shown as squares, are close to the experimentally determined structures from the 
3DRobot_set, which shows that this analysis is robust to small variations in the experimentally 
determined structure.

This can also be seen by plotting all of the structures on one chart. There is a clear cluster
of experimentally determined structures, shown as stars, and their other crystallographic 
structures, shown as squares, in the top right section of the chart. A dotted line has been 
added to the plot to show that the classes of "native" structures and "decoy structures are 
linearly separable. This means that it would be trivial to train a classifier to predict "native"
vs "decoy" for these structures. 

[<img src=./analysis/pca_2dproj_dottedline.png>]()

Furthermore, we investigated what features contributed to both principal components. The table 
below shows the top 6 contributer to principal component 1 and 2.

| Top 6 contributers to PC1  |  Top 6 contributers to PC2 | 
|---|---|
| budeff: total | rosetta - total |
| rosetta - hbond_sr_bb | rosetta - pro_close |
| rosetta - hbond_lr_bb | budeff: steric |
| evoef2: intraR total | rosetta - fa_dun  |
| dfire2 - total | rosetta - hbond_lr_bb |
| budeff: charge | rosetta - fa_sol |











































