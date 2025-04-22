# Smeltz-RED file descriptions 

Descriptions of files found within invitroTKstats/data-raw/Smeltz-RED. 

The unprocessed data file \insertCite{smeltz2023plasma}{invitroTKstats} is used to generate the readily available example Fup RED datasets which can be accessed with `data("fup-red-example")`. The example Fup RED datasets are also used in the function examples for various functions such as `format_fup_red` and `calc_fup_red_point`. 

The processed data files are outputs from generating the data during level-4 processing. They are also used within the "Fup RED" vignette to demonstrate example level-4 files for users. 

### Unprocessed data 
The raw data file contains mass spectrometry measurements of plasma protein binding (PPB) via rapid equilibrium dialysis (RED) for per- and poly-fluorinated alkyl substance (PFAS) samples. Experiments were led by Dr.s Marci Smeltz and Barbara Wetmore.

  * PFAS LC-MS RED Summary 20220709.xlsx - Raw data file containing chemical ID mappings and all experimental data 

### Processed data 
The raw data file was pipelined through `invitroTKstats` to generate the "Fup-RED-example" dataset. The following files are intermediate and output files from the Level 4 processing `calc_fup_red()`.
  
  * Example-fup-RED-Level4Analysis-2025-04-17.RData - Complete Level 4 results (output file from `calc_fup_red()`)
  * Example-fup-RED-Level4.tsv - Level 4 TSV written to one chemical at a time (output file from `calc_fup_red()`)
  * Example-fup-RED-PREJAGS.RData - Arguments given to JAGS (intermediate file from `calc_fup_red()`)
  * Example-fup-RED-Level2-heldout.tsv - Unverified Level 2 samples (intermediate file from `calc_fup_red()`)
  
**CAUTION: We do not anticipate any changes in the unprocessed data files. If there are any new updates to these files, you must include a note log containing the date and explanation of changes at the end of this README.**