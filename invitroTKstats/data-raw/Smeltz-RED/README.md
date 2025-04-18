# Smeltz-RED file descriptions 

Descriptions of files found within invitroTKstats/data-raw/Smeltz-RED. 

### Unprocessed data 
The raw data file contains mass spectrometry measurements of plasma protein binding (PPB) via rapid equilibrium dialysis (RED) for per- and poly-fluorinated alkyl substance (PFAS) samples. Experiments were led by Dr.s Marci Smeltz and Barbara Wetmore.

  * PFAS LC-MS RED Summary 20220709.xlsx - Raw data file containing chemical ID mappings and all experimental data 

### Processed data 
The raw data file was pipelined through `invitroTKstats` to generate the "Fup-RED-example" dataset. The following files are intermediate and output files from the Level 4 processing `calc_fup_red()`.
  
  * Example-fup-RED-Level4Analysis-2025-04-17.RData - Complete Level 4 results (output file from `calc_fup_red()`)
  * Example-fup-RED-Level4.tsv - Level 4 TSV written to one chemical at a time (output file from `calc_fup_red()`)
  * Example-fup-RED-PREJAGS.RData - Arguments given to JAGS (intermediate file from `calc_fup_red()`)
  * Example-fup-RED-Level2-heldout.tsv - Unverified Level 2 samples (intermediate file from `calc_fup_red()`)