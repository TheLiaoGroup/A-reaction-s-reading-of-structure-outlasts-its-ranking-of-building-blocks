# A reaction's reading of structure outlasts its ranking of building blocks

This repository contains the experimental data and the analysis code supporting the Article "A reaction's reading of structure outlasts its ranking of building blocks" (Lin Tan, Yaxian Kong, Chan Zhu, Kai Xue and Kuangbiao Liao).

## Contents

| Folder | Contents |
|---|---|
| `data/` | The 165,131 high-throughput experimentation wells analysed in the Article, and the metadata that define the analysis sets |
| `code/` | The analysis code |
| `source_data/` | `Source_Data.xlsx`: the data underlying Figures 1–4, Extended Data Figs 2–8 and Extended Data Table 1, one sheet per panel |

## Data

`data/experimental_wells.csv` has one row per well:

| Column | Content |
|---|---|
| `well_id` | Well identifier |
| `component` | `existing_release`: the multi-reaction core and the second condition set; `additional_campaign`: the second Buchwald–Hartwig campaign |
| `rxn_class` | Reaction |
| `condition_id` | Reaction condition; identifiers are specific to each reaction |
| `electrophile_smiles`, `nucleophile_smiles` | The two substrates |
| `product_smiles` | The product |
| `plate` | Plate, where recorded |
| `yield` | Measured yield, as a fraction from 0 to 1. Every well was first screened for the product by DESI–TOF–MS; wells with product, and one to three wells without as a check, were quantified by UPLC–UV/MS. A well in which screening found no product and that was not chromatographed is recorded as 0 (63,375 of the 165,131 wells are exactly zero) |

The 165,131 wells comprise:
- 75,870 wells of the multi-reaction core (fifteen reactions, twenty conditions each, 3,982 substrate pairs);
- 71,687 wells of the second condition set (twenty further conditions per reaction, 3,866 substrate pairs);
- 17,574 wells of the second Buchwald–Hartwig campaign (285 conditions, 213 substrate pairs; 212 of them form the connected set from which the molecular contributions are fitted).

Every well enters the analysis: the membership file assigns each of the 165,131 wells to one of these three sets.

`data/metadata/` assigns the wells to the analyses:

| File | Content |
|---|---|
| `analysis_well_membership.csv` | The analysis set of every well |
| `base_condition_roles.csv` | The twenty conditions of each reaction in each condition set: coverage rank, reagents, solvent, temperature and time |
| `base_pair_roles.csv` | The substrate pairs of each reaction in each condition set |
| `extended_conditions.csv` | The 285 conditions of the second Buchwald–Hartwig campaign: reagents, amounts, solvent, temperature and time |
| `calibration_cells.json` | The simulation settings of the interval calibration |

## Reproduce

With Python 3.13.3:

```
pip install -r requirements.txt
python run_analysis.py
```

This runs the analyses of the Article and the Supplementary Information, including the sensitivity analyses, and writes the Supplementary Tables, Extended Data Table 1, the figures and the source data (about 7 minutes).

```
python run_analysis.py --all
```

also runs the two longer analyses: the heterogeneity test of the cross-reaction agreement (about 15 minutes) and the calibration of the intervals on simulated data (about 80 minutes with the default 8 processes; `--workers N` sets the number).

Random draws use fixed seeds, so a rerun reproduces the same outputs. The repository includes an independent implementation of the principal core estimates (`code/independent`).

## Outputs

`output/` is generated and is not included in the repository. Each step writes to its own folder:

| Folder | Contents |
|---|---|
| `output/core`, `output/core_dense`, `output/core_intervals` | Molecule contributions, agreement across condition halves, partner halves and reactions, held-out pair prediction, the second condition set, and their intervals |
| `output/extended`, `output/extended_intervals` | The second Buchwald–Hartwig campaign and the structural profiles |
| `output/partner_profiles` | Structural profiles across partner halves |
| `output/independent` | The independent implementation |
| `output/sensitivity` | The sensitivity analyses of Supplementary Notes 8, 10, 12, 14, 17, 18 and 19 and further values of the main text |
| `output/presentation` | Supplementary Tables 1–11 (`tables`), Extended Data Table 1, the figures (`figures`), the source data (`source_data` and `Source_Data.xlsx`) |
| `output/heterogeneity`, `output/calibration` | The two longer analyses run by `--all` |

Extended Data Fig. 1 shows the reaction schemes and is not computed.

## Figure and table numbers

The code and `Source_Data.xlsx` use the numbering of the Nature-format manuscript. In the Science-format manuscript Figures 1–4 and the Supplementary Notes keep their numbers (panel letters are capitals), and the other items are numbered as follows:

| This repository | Science format |
|---|---|
| Extended Data Fig. *k* (*k* = 1–8) | fig. S*k* |
| Extended Data Table 1 | table S1 |
| Supplementary Tables 1, 2, 3, 4, 5, 6 | tables S3, S2, S9, S11, S7, S8 |
| Supplementary Tables 7, 8, 9, 10, 11 | tables S10, S5, S4, S12, S6 |

## Citation

This release is archived at https://doi.org/10.5281/zenodo.23310573 and is also available at https://github.com/TheLiaoGroup/A-reaction-s-reading-of-structure-outlasts-its-ranking-of-building-blocks. When using the data or code, please cite the Article and this archive.

## Data note

A subset of the amide-coupling measurements was previously reported by Zhang et al. (Chem. Sci. 16, 11809–11822, 2025).
