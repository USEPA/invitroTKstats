# Honda-Caco2 file descriptions 

Descriptions of files found within invitroTKstats/data-raw/Honda-Caco2. 

### Unprocessed data 
The raw data files contains mass spectrometry measurements of membrane permeability (P~app~) from Caco2 cells. The experiments were led by Cyprotex.

  * Caco2.xlsx - Raw data file containing chemical ID mappings 
  * SupTable1-AnalyticalMethods.xlsx - Raw data file containing chemical ID mappings 
  * EPA_Task 10_13_Caco-2 Compiled_LCMSGC_10032017_Data Summary_GZ.xlsm - Original raw data file containing all experimental data 
  * Edited_EPA_Task 10_13_Caco-2 Compiled_LCMSGC_10032017_Data Summary_GZ.xlsm - Edited raw data file containing all experimental data 
  
    * Added headers (column names) to the beginning of each compound chunk
    * Added a 'Test Concentration' column
    * Added a 'Type' column

### Processed data 
The raw data file was pipelined through `invitroTKstats` to generate the "Caco2-example" dataset. The following file is an optional output file from the Level 2 processing `sample_verification()`.
  
  * Examples-Caco-2-Level2.tsv - Level 2 TSV (optional output file from `sample_verification`)