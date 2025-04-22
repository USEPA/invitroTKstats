# Smeltz-Clint file descriptions 

Descriptions of files found within invitroTKstats/data-raw/Smeltz-Clint. 

The unprocessed data file \insertCite{smeltz2023plasma}{invitroTKstats} is used to generate the readily available example Clint dataset which can be accessed with `data("Clint-example")`. 

The processed data files are outputs from generating the data during level-4 processing. They are also used within the "Clint" vignette to demonstrate example level-4 files for users. 

### Unprocessed data 
The raw data file contains mass spectrometry measurements of intrinsic hepatic clearance (Cl~int~) for cryopreserved pooled human hepatocytes. Chemicals were per- and poly-fluorinated alkyl substance (PFAS) samples. The experiments were led by Dr.s Marci Smeltz and Barbara Wetmore. 
  
  * Hep12 Data for Uncertainty Feb2022.xlsx - Raw data file containing chemical ID mappings and all experimental data 

### Processed data 
The raw data file was pipelined through `invitroTKstats` to generate the "Clint-example" dataset. The following files are intermediate and output files from the Level 4 processing `calc_clint()`.
  
  * Example-Clint-Level4Analysis-2025-04-17.RData - Complete Level 4 results (output file from `calc_clint()`)
  * Example-Clint-Level4.tsv - Level 4 TSV written to one chemical at a time (output file from `calc_clint()`)
  * Example-Clint-PREJAGS.RData - Arguments given to JAGS (intermediate file from `calc_clint()`)
  * Example-Clint-Level2-heldout.tsv - Unverified Level 2 samples (intermediate file from `calc_clint()`)
  
CAUTION: We do not anticipate any changes in the unprocessed data files. If there are any new updates to these files, you must include a note log containing the date and explanation of changes at the end of this README. 