#This file is used by roxygen2 to generate man files (documentation) for data
#sets included in the package.

#' Caco-2 Chemical Information Example Data set
#'
#' The chemical ID mapping information from tandem mass spectrometry (MS/MS) measurements 
#' of Caco-2 assay-specific data. 
#' This data set contains 520 unique compounds/chemicals.
#' 
#' @name caco2_cheminfo
#' @aliases caco2_cheminfo
#' @docType data
#' @format A chemical info data.frame with 554 rows and 7 variables: \describe{
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard)}
#' \item{\code{PREFERRED_NAME}}{Preferred compound name from the CompTox Chemicals Dashboard (CCD)}
#' \item{\code{CASRN}}{CAS Registry Number of the test compound}
#' \item{\code{MOLECULAR_FORMULA}}{Molecular formula of the test compound}
#' \item{\code{AVERAGE_MASS}}{Molecular weight of the compound in daltons}
#' \item{\code{QSAR_READY_SMILES}}{SMILES (Simplified molecular-input line-entry system) chemical structure description.}
#' \item{\code{test_article}}{Compound ID used in the laboratory}
#' }
#'
"caco2_cheminfo"

#' Caco-2 Level-0 Example Data set
#' 
#' A subset of tandem mass spectrometry (MS/MS) measurements of Caco-2 assay-specific
#' data. This subset contains samples for 3 test analytes/compounds. 
#' 
#' @name caco2_L0
#' @aliases caco2_L0
#' @docType data
#' @format A level-0 data.frame with 48 rows and 17 variables: \describe{
#' \item{\code{Compound}}{Compound name}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard)}
#' \item{\code{Lab.Compound.ID}}{Compound ID used in the laboratory}
#' \item{\code{Date}}{Date MS/MS assay data acquired from instrument}
#' \item{\code{Sample}}{Sample Name}
#' \item{\code{Type}}{Type of Caco-2 sample}
#' \item{\code{Compound.Conc}}{Expected (or nominal) concentration of analyte (for calibration curve)}
#' \item{\code{Peak.Area}}{Peak area of analyte (target compound)}
#' \item{\code{ISTD.Peak.Area}}{Peak area of internal standard (pixels)}
#' \item{\code{ISTD.Name}}{Name of compound used as internal standard (ISTD)}
#' \item{\code{Analysis.Params}}{General description of chemical analysis method}
#' \item{\code{Level0.File}}{Name of data file from laboratory that was used to compile level-0 data.frame}
#' \item{\code{Level0.Sheet}}{Name of "sheet" (for Excel workbooks) from which the laboratory data were read}
#' \item{\code{Direction}}{Direction of the Caco-2 permeability experiment}
#' \item{\code{Vol.Donor}}{The media volume (in cm^3) of the donor portion of the Caco-2 experimental well}
#' \item{\code{Vol.Receiver}}{The media volume (in cm^3) of the receiver portion of the Caco-2 experimental well}
#' \item{\code{Dilution.Factor}}{Number of times the sample was diluted}
#' }
#'
"caco2_L0"

#' Caco-2 Level-1 Example Data set
#'
#' A subset of tandem mass spectrometry (MS/MS) measurements of Caco-2 assay-specific
#' data. This subset contains samples for 3 test analytes/compounds.
#' 
#' @name caco2_L1
#' @aliases caco2_L1
#' @docType data
#' @format A level-1 data.frame with 48 rows and 28 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample name as described in the laboratory}
#' \item{\code{Date}}{Date MS/MS assay data acquired from instrument}
#' \item{\code{Compound.Name}}{Compound name}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of Caco-2 sample}
#' \item{\code{Direction}}{Direction of the Caco-2 permeability experiment}
#' \item{\code{Dilution.Factor}}{Number of times the sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{Biological.Replicates}}{Identifier for measurements of multiple samples with the same analyte}
#' \item{\code{Technical.Replicates}}{Identifier for repeated measurements of one sample of a compound}
#' \item{\code{Test.Compound.Conc}}{Measured concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{Test.Nominal.Conc}}{Expected initial concentration of chemical added to donor side (uM)}
#' \item{\code{Time}}{Time when sample was measured (h)}
#' \item{\code{ISTD.Name}}{Name of compound used as internal standard (ISTD)}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (pixels)}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Membrane.Area}}{The area of the Caco-2 monolayer.}
#' \item{\code{Vol.Donor}}{The media volume (in cm^3) of the donor portion of the Caco-2 experimental well}
#' \item{\code{Vol.Receiver}}{The media volume (in cm^3) of the receiver portion of the Caco-2 experimental well}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Additional information}
#' \item{\code{Level0.File}}{Name of data file from laboratory that was used to compile level-0 data.frame)}
#' \item{\code{Level0.Sheet}}{Name of "sheet" (for Excel workbooks) from which the laboratory data were read}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' }
#'
"caco2_L1"

#' Caco-2 Level-2 Example Data set
#'
#' A subset of tandem mass spectrometry (MS/MS) measurements of Caco-2 assay-specific
#' data. This subset contains samples for 3 test analytes/compounds.
#' 
#' @name caco2_L2
#' @aliases caco2_L2
#' @docType data
#' @format A level-2 data.frame with 48 rows and 29 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample name as described in the laboratory}
#' \item{\code{Date}}{Date MS/MS assay data acquired from instrument}
#' \item{\code{Compound.Name}}{Compound name}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of Caco-2 sample}
#' \item{\code{Direction}}{Direction of the Caco-2 permeability experiment}
#' \item{\code{Dilution.Factor}}{Number of times the sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{Biological.Replicates}}{Identifier for measurements of multiple samples with the same analyte}
#' \item{\code{Technical.Replicates}}{Identifier for repeated measurements of one sample of a compound}
#' \item{\code{Test.Compound.Conc}}{Measured concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{Test.Nominal.Conc}}{Expected initial concentration of chemical added to donor side (uM)}
#' \item{\code{Time}}{Time when sample was measured (h)}
#' \item{\code{ISTD.Name}}{Name of compound used as internal standard (ISTD)}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (pixels)}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Membrane.Area}}{The area of the Caco-2 monolayer.}
#' \item{\code{Vol.Donor}}{The media volume (in cm^3) of the donor portion of the Caco-2 experimental well}
#' \item{\code{Vol.Receiver}}{The media volume (in cm^3) of the receiver portion of the Caco-2 experimental well}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Additional information}
#' \item{\code{Level0.File}}{Name of data file from laboratory that was used to compile level-0 data.frame)}
#' \item{\code{Level0.Sheet}}{Name of "sheet" (for Excel workbooks) from which the laboratory data were read}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' \item{\code{Verified}}{If "Y", then sample is included in the analysis. (Any other causes the data to be ignored.)}
#' }
#'
"caco2_L2"

#' Caco-2 Level-3 Example Data set
#'
#' A subset of tandem mass spectrometry (MS/MS) measurements of Caco-2 assay-specific
#' data. This subset contains samples for 3 test analytes/compounds.
#' 
#' @name caco2_L3
#' @aliases caco2_L3
#' @docType data
#' @format A level-3 data.frame with 3 rows and 20 variables: \describe{
#' \item{\code{Compound.Name}}{Compound name}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard)}
#' \item{\code{Time}}{Time when sample was measured (h)}
#' \item{\code{Membrane.Area}}{The area of the Caco-2 monolayer}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{C0_A2B}}{Initial concentration in the apical side}
#' \item{\code{dQdt_A2B}}{Rate of permeation from the apical to the basolateral side}
#' \item{\code{Papp_A2B}}{Apparent membrane permeability from the apical to the basolateral side}
#' \item{\code{Frec_A2B.vec}}{Fraction of the initial compound in the apical side recovered in the basolateral side (collapsed numeric vector, values for replicates separated by a "|")}
#' \item{\code{Frec_A2B.mean}}{Mean of fraction recovered values in the apical to basolateral direction}
#' \item{\code{Recovery_Class_A2B.vec}}{Recovery classification of fraction recovered values in the apical to basolateral direction (collapsed character vector, values for replicates separated by a "|")}
#' \item{\code{Recovery_Class_A2B.mean}}{Recovery classification of mean fraction recovered in the apical to basolateral direction}
#' \item{\code{C0_B2A}}{Initial concentration in the basolateral side}
#' \item{\code{dQdt_B2A}}{Rate of permeation from the basolateral to the apical side}
#' \item{\code{Papp_B2A}}{Apparent membrane permeability from the basolateral to the apical side}
#' \item{\code{Frec_B2A.vec}}{Fraction of the initial compound in the basolateral side recovered in the apical side (collapsed numeric vector, values for replicates separated by a "|")}
#' \item{\code{Frec_B2A.mean}}{Mean of fraction recovered values in the basolateral to apical direction}
#' \item{\code{Recovery_Class_B2A.vec}}{Recovery classification of fraction recovered values in the basolateral to apical direction (collapsed character vector, values for replicates separated by a "|")}
#' \item{\code{Recovery_Class_B2A.mean}}{Recovery classification of mean fraction recovered in the basolateral to apical direction}
#' \item{\code{Refflux}}{Efflux ratio}
#' }
#'
"caco2_L3"

#' Fup UC Chemical Information Example Data set
#' 
#' The chemical ID mapping information from mass spectrometry measurements of plasma protein binding (PPB) via 
#' ultracentrifugation (UC) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples. Experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. 
#' This data set contains 75 unique compounds/chemicals.
#'
#' @name fup_uc_cheminfo
#' @aliases fup_uc_cheminfo
#' @docType data
#' @format A chemical info data.frame with 75 rows and 4 variables: \describe{
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Chemical Name (Common Abbreviation)}}{Name of the test analyte/compound and abbreviation used by the lab as the compound ID}
#' \item{\code{Compound}}{Name of the test analyte/compound}
#' \item{\code{Chem.Lab.ID}}{Common abbreviation of the test analyte/compound as described in the laboratory}
#' }
#'
#' @references
#' \insertRef{howard2010plasma}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#' @keywords data
"fup_uc_cheminfo"

#' Fup UC Level-0 Example Data set
#' 
#' Mass Spectrometry measurements of plasma protein binding (PPB) via 
#' ultracentrifugation (UC) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples. Experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. 
#' This data set is a subset of experimental data containing samples for 
#' 3 test analytes/compounds.
#'
#' @name fup_uc_L0
#' @aliases fup_uc_L0
#' @docType data
#' @format A level-0 data.frame with 240 rows and 17 variables: \describe{
#' \item{\code{Compound}}{Name of the test analyte/compound}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Lab.Compound.ID}}{Compound as described in the laboratory}
#' \item{\code{Date}}{Date the sample was added to the MS analyzer}
#' \item{\code{Sample}}{Sample description used in the laboratory}
#' \item{\code{Type}}{Type of UC sample, annotated by the laboratory}
#' \item{\code{Compound.Conc}}{Expected (or nominal) concentration of analyte (for calibration curve)}
#' \item{\code{Peak.Area}}{Peak area of analyte (target compound)}
#' \item{\code{ISTD.Peak.Area}}{Peak area of internal standard (ISTD) compound (pixels)}
#' \item{\code{ISTD.Name}}{Name of the internal standard (ISTD) analyte/compound}
#' \item{\code{Analysis.Params}}{Column contains the retention time} 
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level-0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level-0 sample data was extracted}
#' \item{\code{Sample.Text}}{Additional notes on the sample}
#' \item{\code{Sample.Type}}{Type of UC sample in \code{invitroTKstats} package annotations}
#' \item{\code{Dilution.Factor}}{Number of times the sample was diluted}
#' \item{\code{Replicate}}{Identifier for repeated measurements of one sample of a compound}
#' }
#'
#' @references
#' \insertRef{howard2010plasma}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#' @keywords data
"fup_uc_L0"

#' Fup UC Level-1 Example Data set
#' 
#' Mass Spectrometry measurements of plasma protein binding (PPB) via 
#' ultracentrifugation (UC) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples. Experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. 
#' This data set is a subset of experimental data containing samples for 
#' 3 test analytes/compounds.
#'
#' @name fup_uc_L1
#' @aliases fup_uc_L1
#' @docType data
#' @format A level-1 data.frame with 240 rows and 23 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample description used in the laboratory}
#' \item{\code{Date}}{Date the sample was added to the MS analyzer}
#' \item{\code{Compound.Name}}{Name of the test analyte/compound}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of UC sample}
#' \item{\code{Dilution.Factor}}{Number of times the sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{ISTD.Name}}{Name of the internal standard (ISTD) analyte/compound}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (ISTD) compound (pixels)}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level-0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level-0 sample data was extracted}
#' \item{\code{Test.Compound.Conc}}{Measured concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{Test.Nominal.Conc}}{Expected initial concentration of chemical added to T1 sample (uM)}
#' \item{\code{Biological.Replicates}}{Identifier for measurements of multiple samples with the same analyte}
#' \item{\code{Technical.Replicates}}{Identifier for repeated measurements of one sample of a compound}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' }
#'
#' @references
#' \insertRef{howard2010plasma}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#' @keywords data
"fup_uc_L1"

#' Fup UC Level-2 Example Data set
#' 
#' Mass Spectrometry measurements of plasma protein binding (PPB) via 
#' ultracentrifugation (UC) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples. Experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. 
#' This data set is a subset of experimental data containing samples for 
#' 3 test analytes/compounds.
#'
#' @name fup_uc_L2
#' @aliases fup_uc_L2
#' @docType data
#' @format A level-2 data.frame with 240 rows and 24 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample description used in the laboratory}
#' \item{\code{Date}}{Date the sample was added to the MS analyzer}
#' \item{\code{Compound.Name}}{Name of the test analyte/compound}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of UC sample}
#' \item{\code{Dilution.Factor}}{Number of times the sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{ISTD.Name}}{Name of the internal standard (ISTD) analyte/compound}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (ISTD) compound (pixels)}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level-0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level-0 sample data was extracted}
#' \item{\code{Test.Compound.Conc}}{Measured concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{Test.Nominal.Conc}}{Expected initial concentration of chemical added to T1 sample (uM)}
#' \item{\code{Biological.Replicates}}{Identifier for measurements of multiple samples with the same analyte}
#' \item{\code{Technical.Replicates}}{Identifier for repeated measurements of one sample of a compound}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' \item{\code{Verified}}{If "Y", then sample is included in the analysis. (Any other value causes the data to be ignored.)}
#' }
#'
#' @references
#' \insertRef{howard2010plasma}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#' @keywords data
"fup_uc_L2"

#' Fup UC Level-3 Example Data set
#' 
#' Mass Spectrometry measurements of plasma protein binding (PPB) via 
#' ultracentrifugation (UC) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples. Experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. 
#' This data set is a subset of experimental data containing samples for 
#' 3 test analytes/compounds.
#'
#' @name fup_uc_L3
#' @aliases fup_uc_L3
#' @docType data
#' @format A level-3 data.frame with 3 rows and 5 variables: \describe{
#' \item{\code{Compound.Name}}{Name of the test analyte/compound}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{Fup}}{Fraction unbound in plasma}
#' }
#'
#' @references
#' \insertRef{howard2010plasma}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#' @keywords data
"fup_uc_L3"

#' Fup UC Level-4 Example Data set
#' 
#' Mass Spectrometry measurements of plasma protein binding (PPB) via 
#' ultracentrifugation (UC) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples. Experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. 
#' This data set is a subset of experimental data containing samples for 
#' 3 test analytes/compounds.
#'
#' @name fup_uc_L4
#' @aliases fup_uc_L4
#' @docType data
#' @format A level-4 data.frame with 3 rows and 10 variables: \describe{
#' \item{\code{Compound}}{Name of the test analyte/compound}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Fstable.Med}}{Median stability fraction}
#' \item{\code{Fstable.Low}}{2.5th quantile of stability fraction}
#' \item{\code{Fstable.High}}{97.5th quantile of stability fraction}
#' \item{\code{Fup.Med}}{Median fraction unbound in plasma}
#' \item{\code{Fup.Low}}{2.5th quantile of fraction unbound in plasma}
#' \item{\code{Fup.High}}{97.5th quantile of fraction unbound in plasma}
#' \item{\code{Fup.point}}{Point estimate of fraction unbound in plasma}
#' }
#'
#' @references
#' \insertRef{howard2010plasma}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#' @keywords data
"fup_uc_L4"

#' Fup UC Level-4 PREJAGS arguments
#' 
#' The arguments given to JAGS for the tested compound during level-4 processing of mass spectrometry measurements of plasma protein binding (PPB) via 
#' ultracentrifugation (UC) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples. Experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. This list is overwritten for each tested 
#' compound. Therefore, only contains arguments given to JAGS for the last tested compound. 
#' 
#' @name fup_uc_PREJAGS
#' @aliases fup_uc_PREJAGS
#' @docType data
#' @format A named list with 10 elements: \describe{
#' \item{\code{Num.cal}}{Unique number of \code{Calibration} values for the tested compound}
#' \item{\code{Num.obs}}{Total number of observations for the tested compound}
#' \item{\code{Response.obs}}{\code{Response} of all samples for the tested compound}
#' \item{\code{obs.conc}}{Indices of the \code{Test.Compound.Conc} values that corresponds to all samples' \code{Test.Compound.Conc} for the tested compound.}
#' \item{\code{obs.cal}}{Indices of the unique \code{Calibration} values that corresponds to all samples' \code{Calibration} for the tested compound.}
#' \item{\code{Conc}}{\code{Test.Compound.Conc} of the "CC" sample types + three placeholder concentrations ("T1", "T5", "AF") per \code{Biological.Replicates} series}
#' \item{\code{Num.cc.obs}}{Number of "CC" sample types for the tested compound}
#' \item{\code{Num.series}}{Unique number of \code{Biological.Replicates} series}
#' \item{\code{Dilution.Factor}}{\code{Dilution.Factor} of all samples for the tested compound (number of times the sample was diluted)}
#' \item{\code{Test.Nominal.Conc}}{Unique \code{Test.Nominal.Conc} values (expected initial concentration) of all samples for the tested compound}
#' }
#' 
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#'
"fup_uc_PREJAGS"

#' Fup UC Level-2 Heldout Example Data set
#' 
#' The unverified level-2 samples from mass spectrometry measurements of plasma protein binding (PPB) via 
#' ultracentrifugation (UC) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples. Experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. 
#' This data set is a subset of experimental data containing samples for 
#' 0 test analytes/compounds. No data samples are unverified.
#'
#' @name fup_uc_L2_heldout
#' @aliases fup_uc_L2_heldout
#' @docType data
#' @format A level-2 data.frame with 0 rows and 24 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample description used in the laboratory}
#' \item{\code{Date}}{Date the sample was added to the MS analyzer}
#' \item{\code{Compound.Name}}{Name of the test analyte/compound}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of UC sample}
#' \item{\code{Dilution.Factor}}{Number of times the sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{ISTD.Name}}{Name of the internal standard (ISTD) analyte/compound}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (ISTD) compound (pixels)}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level-0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level-0 sample data was extracted}
#' \item{\code{Test.Compound.Conc}}{Measured concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{Test.Nominal.Conc}}{Expected initial concentration of chemical added to T1 sample (uM)}
#' \item{\code{Biological.Replicates}}{Identifier for measurements of multiple samples with the same analyte}
#' \item{\code{Technical.Replicates}}{Identifier for repeated measurements of one sample of a compound}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' \item{\code{Verified}}{If "Y", then sample is included in the analysis. (Any other value causes the data to be ignored.)}
#' }
#'
#' @references
#' \insertRef{howard2010plasma}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#' @keywords data
"fup_uc_L2_heldout"

#' Fup RED Chemical Information Example Data set
#'
#' The chemical ID mapping information from mass spectrometry measurements of plasma protein binding (PPB) via rapid 
#' equilibrium dialysis (RED) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples. Experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}.
#' This data set contains 26 unique compounds/chemicals.
#' 
#' @name fup_red_cheminfo
#' @aliases fup_red_cheminfo
#' @docType data
#' @format A chemical info data.frame with 26 rows and 4 variables: \describe{
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{NAME (Abbreviation)}}{Name of the test analyte/compound and abbreviation used by the lab as the compound ID}
#' \item{\code{Compound}}{Name of the test analyte/compound}
#' \item{\code{Chem.Lab.ID}}{Abbreviation of the test analyte/compound as described in the laboratory}
#' }
#' 
#' @references
#' \insertRef{waters2008validation}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#'
"fup_red_cheminfo"

#' Fup RED Level-0 Example Data set
#'
#' Mass Spectrometry measurements of plasma protein binding (PPB) via rapid 
#' equilibrium dialysis (RED) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples. Experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}.
#' This data set is a subset of experimental data containing samples for 
#' 3 test analytes/compounds.
#' 
#' @name fup_red_L0
#' @aliases fup_red_L0
#' @docType data
#' @format A level-0 data.frame with 660 rows and 18 variables: \describe{
#' \item{\code{Compound}}{Name of the test analyte/compound}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Lab.Compound.ID}}{Compound as described in the laboratory}
#' \item{\code{Date}}{Date the sample was added to the MS analyzer}
#' \item{\code{Sample}}{Sample description used in the laboratory}
#' \item{\code{Type}}{Type of RED sample, annotated by the laboratory}
#' \item{\code{Compound.Conc}}{Expected (or nominal) concentration of analyte (for calibration curve)}
#' \item{\code{Peak.Area}}{Peak area of analyte (target compound)}
#' \item{\code{ISTD.Peak.Area}}{Peak area of internal standard (ISTD) compound (pixels)}
#' \item{\code{ISTD.Name}}{Name of the internal standard (ISTD) analyte/compound}
#' \item{\code{Analysis.Params}}{Column contains the retention time}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level-0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level-0 sample data was extracted}
#' \item{\code{Sample Text}}{Additional notes on the sample}
#' \item{\code{Sample.Type}}{Type of RED sample in \code{invitroTKstats} package annotations}
#' \item{\code{Replicate}}{Identifier for repeated measurements of one sample of a compound}
#' \item{\code{Time}}{Time when the sample was measured - in hours (h)}
#' \item{\code{Dilution.Factor}}{Number of times the sample was diluted}
#' }
#' 
#' @references
#' \insertRef{waters2008validation}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#'
"fup_red_L0"

#' Fup RED Level-1 Example Data set
#'
#' Mass Spectrometry measurements of plasma protein binding (PPB) via rapid 
#' equilibrium dialysis (RED) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples. Experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}.
#' This data set is a subset of experimental data containing samples for 
#' 3 test analytes/compounds.
#' 
#' @name fup_red_L1
#' @aliases fup_red_L1
#' @docType data
#' @format A level-1 data.frame with 636 rows and 25 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample description used in the laboratory}
#' \item{\code{Date}}{Date the sample was added to the MS analyzer}
#' \item{\code{Compound.Name}}{Name of the test analyte/compound}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of RED sample}
#' \item{\code{Dilution.Factor}}{Number of times the sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{ISTD.Name}}{Name of the internal standard (ISTD) analyte/compound}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (ISTD) compound (pixels)}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level-0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level-0 sample data was extracted}
#' \item{\code{Time}}{Time when the sample was measured - in hours (h)}
#' \item{\code{Test.Compound.Conc}}{Measured concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{Test.Nominal.Conc}}{Expected initial concentration of chemical added to RED plate (uM)}
#' \item{\code{Percent.Physiologic.Plasma}}{Percent of physiological plasma concentration in RED plate (in percent)}
#' \item{\code{Technical.Replicates}}{Identifier for repeated measurements of a sample of a compound}
#' \item{\code{Biological.Replicates}}{Identifier for measurements of multiple samples with the same analyte}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' }
#'
#' @references
#' \insertRef{howard2010plasma}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#' @keywords data
"fup_red_L1"

#' Fup RED Level-2 Example Data set
#'
#' Mass Spectrometry measurements of plasma protein binding (PPB) via rapid 
#' equilibrium dialysis (RED) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples. Experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}.
#' This data set is a subset of experimental data containing samples for 
#' 3 test analytes/compounds.
#' 
#' @name fup_red_L2
#' @aliases fup_red_L2
#' @docType data
#' @format A level-2 data.frame with 636 rows and 26 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample description used in the laboratory}
#' \item{\code{Date}}{Date the sample was added to the MS analyzer}
#' \item{\code{Compound.Name}}{Name of the test analyte/compound}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of RED sample}
#' \item{\code{Dilution.Factor}}{Number of times the sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{ISTD.Name}}{Name of the internal standard (ISTD) analyte/compound}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (ISTD) compound (pixels)}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level-0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level-0 sample data was extracted}
#' \item{\code{Time}}{Time when the sample was measured - in hours (h)}
#' \item{\code{Test.Compound.Conc}}{Measured concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{Test.Nominal.Conc}}{Expected initial concentration of chemical added to RED plate (uM)}
#' \item{\code{Percent.Physiologic.Plasma}}{Percent of physiological plasma concentration in RED plate (in percent)}
#' \item{\code{Technical.Replicates}}{Identifier for repeated measurements of one sample of a compound}
#' \item{\code{Biological.Replicates}}{Identifier for measurements of multiple samples with the same analyte}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' \item{\code{Verified}}{If, "Y" then sample is included in the analysis. (Any other value causes the data to be ignored.)}
#' }
#'
#' @references
#' 
#' \insertRef{waters2008validation}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
"fup_red_L2"

#' Fup RED Level-3 Example Data set
#'
#' Mass Spectrometry measurements of plasma protein binding (PPB) via rapid 
#' equilibrium dialysis (RED) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples. Experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}.
#' This data set is a subset of experimental data containing samples for 
#' 3 test analytes/compounds.
#' 
#' @name fup_red_L3
#' @aliases fup_red_L3
#' @docType data
#' @format A level-3 data.frame with 3 rows and 4 variables: \describe{
#' \item{\code{Compound.Name}}{Name of the test analyte/compound}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{Fup}}{Fraction unbound in plasma}
#' }
#'
#' @references
#' 
#' \insertRef{waters2008validation}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
"fup_red_L3"

#' Fup RED Level-4 Example Data set
#'
#' Mass Spectrometry measurements of plasma protein binding (PPB) via rapid 
#' equilibrium dialysis (RED) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples. Experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}.
#' This data set is a subset of experimental data containing samples for 
#' 3 test analytes/compounds.
#' 
#' @name fup_red_L4
#' @aliases fup_red_L4
#' @docType data
#' @format A level-4 data.frame with 3 rows and 7 variables: \describe{
#' \item{\code{Compound.Name}}{Name of the test analyte/compound}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Fup.point}}{Point estimate of fraction unbound in plasma}
#' \item{\code{Fup.Med}}{Median fraction unbound in plasma}
#' \item{\code{Fup.Low}}{2.5th quantile of fraction unbound in plasma}
#' \item{\code{Fup.High}}{97.5th quantile of fraction unbound in plasma}
#' }
#'
#' @references
#' 
#' \insertRef{waters2008validation}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
"fup_red_L4"

#' Fup RED Level-4 PREJAGS arguments
#' 
#' The arguments given to JAGS for the tested compound during level-4 processing of mass spectrometry measurements of plasma protein binding (PPB) via rapid 
#' equilibrium dialysis (RED) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples. Experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. This list is overwritten for each tested
#' compound. Therefore, only contains arguments given to JAGS for the last tested compound. 
#' 
#' @name fup_red_PREJAGS
#' @aliases fup_red_PREJAGS
#' @docType data
#' @format A named list with 33 elements: \describe{
#' \item{\code{Test.Nominal.Conc}}{Unique \code{Test.Nominal.Conc} values (expected initial concentration) for the tested compound}
#' \item{\code{Num.cal}}{Unique number of \code{Calibration} values for the tested compound}
#' \item{\code{Physiological.Protein.Conc}}{The assumed physiological protein concentration 
#' for plasma protein binding calculations. (Defaults to 70/(66.5*1000)*1000000.
#' According to \insertCite{berg2011pathology;textual}{invitroTKstats}: 60-80 mg/mL, albumin is 66.5 kDa,
#' assume all protein is albumin to estimate default in uM.)}
#' \item{\code{Assay.Protein.Perecent}}{\code{Percent.Physiologic.Plasma} values for each "Plasma" sample type replicate group}
#' \item{\code{Num.Plasma.Blank.obs}}{Number of "Plasma.Blank" sample types for the tested compound}
#' \item{\code{Plasma.Blank.obs}}{\code{Response} of the "Plasma.Blank" sample types for the tested compound}
#' \item{\code{Plasma.Blank.cal}}{Indices of the unique \code{Calibration} values that corresponds to the "Plasma.Blank" sample types' \code{Calibration} for the tested compound}
#' \item{\code{Plasma.Blank.df}}{Unique \code{Dilution Factor} of the "Plasma.Blank" sample types for the tested compound}
#' \item{\code{Plasma.Blank.rep}}{Integer representing "Plasma.Blank" replicate group for the tested compound}
#' \item{\code{Num.NoPlasma.Blank.obs}}{Number of "NoPlasma.Blank" sample types for the tested compound}
#' \item{\code{NoPlasma.Blank.obs}}{\code{Response} of the "NoPlasma.Blank" sample types for the tested compound}
#' \item{\code{NoPlasma.Blank.cal}}{Indices of the unique \code{Calibration} values that corresponds to the "NoPlasma.Blank" sample types' \code{Calibration} for the tested compound}
#' \item{\code{NoPlasma.Blank.df}}{Unique \code{Dilution Factor} of the "NoPlasma.Blank" sample types for the tested compound}
#' \item{\code{Num.CC.obs}}{Number of "CC" sample types with non-NA \code{Test.Compound.Conc} values for the tested compound}
#' \item{\code{CC.conc}}{\code{Test.Compound.Conc} (non-NA) of the "CC" sample types for the tested compound}
#' \item{\code{CC.obs}}{\code{Response} of the "CC" sample types with non-NA \code{Test.Compound.Conc} for the tested compound}
#' \item{\code{CC.cal}}{Indices of the unique \code{Calibration} values that corresponds to the "CC" sample types' \code{Calibration} for the tested compound}
#' \item{\code{CC.df}}{Unique \code{Dilution Factor} of the "NoPlasma.Blank" sample types for the tested compound}
#' \item{\code{Num.T0.obs}}{Number of "T0" sample types for the tested compound}
#' \item{\code{T0.obs}}{\code{Response} of the "T0" sample types for the tested compound}
#' \item{\code{T0.cal}}{Indices of the unique \code{Calibration} values that corresponds to the "T0" sample types' \code{Calibration} for the tested compound}
#' \item{\code{T0.df}}{Unique \code{Dilution Factor} of the "T0" sample types for the tested compound}
#' \item{\code{Num.rep}}{Unique number of (\code{Calibration} + \code{Technical.Replicates}) combinations for "PBS" and "Plasma" sample types for the tested compound}
#' \item{\code{Num.PBS.obs}}{Number of "PBS" sample types for the tested compound}
#' \item{\code{PBS.obs}}{\code{Response} of the "PBS" sample types for the tested compound}
#' \item{\code{PBS.cal}}{Indices of the unique \code{Calibration} values that corresponds to the "PBS" sample types' \code{Calibration} for the tested compound}
#' \item{\code{PBS.df}}{Unique \code{Dilution Factor} of the "PBS" sample types for the tested compound}
#' \item{\code{PBS.rep}}{Integer representing "PBS" replicate group for the tested compound}
#' \item{\code{Num.Plasma.obs}}{Number of "Plasma" sample types for the tested compound}
#' \item{\code{Plasma.obs}}{\code{Response} of the "Plasma" sample types for the tested compound}
#' \item{\code{Plasma.cal}}{Indices of the unique \code{Calibration} values that corresponds to the "Plasma" sample types' \code{Calibration} for the tested compound}
#' \item{\code{Plasma.df}}{Unique \code{Dilution Factor} of the "Plasma" sample types for the tested compound}
#' \item{\code{Plasma.rep}}{Integer representing "Plasma" replicate group for the tested compound}
#' }
#' 
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#'
"fup_red_PREJAGS"

#' Fup RED Level-2 Heldout Example Data set
#'
#' The unverified level-2 samples from mass spectrometry measurements of plasma protein binding (PPB) via rapid 
#' equilibrium dialysis (RED) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples. Experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}.
#' This data set is a subset of experimental data containing samples for 
#' 0 test analytes/compounds. No data samples are unverified. 
#' 
#' @name fup_red_L2_heldout
#' @aliases fup_red_L2_heldout
#' @docType data
#' @format A level-2 data.frame with 0 rows and 26 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample description used in the laboratory}
#' \item{\code{Date}}{Date the sample was added to the MS analyzer}
#' \item{\code{Compound.Name}}{Name of the test analyte/compound}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of RED sample}
#' \item{\code{Dilution.Factor}}{Number of times the sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{ISTD.Name}}{Name of the internal standard (ISTD) analyte/compound}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (ISTD) compound (pixels)}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level-0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level-0 sample data was extracted}
#' \item{\code{Time}}{Time when the sample was measured - in hours (h)}
#' \item{\code{Test.Compound.Conc}}{Measured concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{Test.Nominal.Conc}}{Expected initial concentration of chemical added to RED plate (uM)}
#' \item{\code{Percent.Physiologic.Plasma}}{Percent of physiological plasma concentration in RED plate (in percent)}
#' \item{\code{Technical.Replicates}}{Identifier for repeated measurements of one sample of a compound}
#' \item{\code{Biological.Replicates}}{Identifier for measurements of multiple samples with the same analyte}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' \item{\code{Verified}}{If "Y", then sample is included in the analysis. (Any other value causes the data to be ignored.)}
#' }
#'
#' @references
#' 
#' \insertRef{waters2008validation}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
"fup_red_L2_heldout"

#' Clint Chemical Information Example Data set
#'
#' The chemical ID mapping information from mass spectrometry measurements of intrinsic
#' hepatic clearance (Clint) for cryopreserved pooled human hepatocytes.
#' Chemicals were per- and poly-fluorinated alkyl substance
#' (PFAS) samples. The experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. This data set contains 7 unique compounds/chemicals. 
#' 
#' @name clint_cheminfo
#' @aliases clint_cheminfo
#' @docType data
#' @format A chemical info data.frame with 7 rows and 6 variables: \describe{
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Analyte Name}}{Name of the test analyte/compound and the name used by the laboratory}
#' \item{\code{Internal Standard}}{Name of the internal standard (ISTD)}
#' \item{\code{Mix}}{Mix used for the sample}
#' \item{\code{Compound}}{Name of the test analyte/compound}
#' \item{\code{Chem.Lab.ID}}{Compound as described in the chemistry laboratory}
#' }
#' 
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#'
"clint_cheminfo"

#' Clint Level-0 Example Data set
#'
#' Mass Spectrometry measurements of intrinsic hepatic clearance (Clint) for cryopreserved 
#' pooled human hepatocytes. Chemicals were per- and poly-fluorinated alkyl substance
#' (PFAS) samples. The experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. This data set is a subset of 
#' experimental data containing samples for 3 test analytes/compounds.
#' 
#' @name clint_L0
#' @aliases clint_L0
#' @docType data
#' @format A level-0 data.frame with 247 rows and 16 variables: \describe{
#' \item{\code{Compound}}{Name of the test analyte/compound}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Lab.Compound.ID}}{Compound as described in the laboratory}
#' \item{\code{Date}}{Date the sample was added to the MS analyzer}
#' \item{\code{Sample}}{Sample description used in the laboratory}
#' \item{\code{Type}}{Type of Clint sample}
#' \item{\code{Compound.Conc}}{Expected (or nominal) concentration of analyte (for calibration curve)}
#' \item{\code{Peak.Area}}{Peak area of analyte (target compound)}
#' \item{\code{ISTD.Peak.Area}}{Peak area of internal standard (ISTD) compound (pixels)}
#' \item{\code{ISTD.Name}}{Name of the internal standard (ISTD) analyte/compound}
#' \item{\code{Analysis.Params}}{Column contains the retention time}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level-0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level-0 sample data was extracted}
#' \item{\code{Sample.Text}}{Additional notes on the sample}
#' \item{\code{Time}}{Time when the sample was measured - in hours (h)}
#' \item{\code{Dilution.Factor}}{Number of times the sample was diluted}
#' }
#' 
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#'
"clint_L0"

#' Clint Level-1 Example Data set
#'
#' Mass Spectrometry measurements of intrinsic hepatic clearance (Clint) for cryopreserved 
#' pooled human hepatocytes. Chemicals were per- and poly-fluorinated alkyl substance
#' (PFAS) samples. The experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. This data set is a subset of 
#' experimental data containing samples for 3 test analytes/compounds.
#' 
#' @name clint_L1
#' @aliases clint_L1
#' @docType data
#' @format A level-1 data.frame with 229 rows and 24 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample description used in the laboratory}
#' \item{\code{Date}}{Date the sample was added to the MS analyzer}
#' \item{\code{Compound.Name}}{Name of the test analyte/compound}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of Clint sample}
#' \item{\code{Dilution.Factor}}{Number of times the sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{ISTD.Name}}{Name of the internal standard (ISTD) analyte/compound}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (pixels)}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level-0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level-0 sample data was extracted}
#' \item{\code{Time}}{Time when the sample was measured - in hours (h)}
#' \item{\code{Test.Compound.Conc}}{Measured concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{Test.Nominal.Conc}}{Expected initial concentration of chemical added to well (uM)}
#' \item{\code{Hep.Density}}{The density (units of millions of hepatocytes per mL) hepatocytes in the \emph{in vitro} incubation}
#' \item{\code{Biological.Replicates}}{Identifier for measurements of multiple samples with the same analyte}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' }
#' 
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#'
"clint_L1"

#' Clint Level-2 Example Data set
#'
#' Mass Spectrometry measurements of intrinsic hepatic clearance (Clint) for cryopreserved 
#' pooled human hepatocytes. Chemicals were per- and poly-fluorinated alkyl substance
#' (PFAS) samples. The experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. This data set is a subset of 
#' experimental data containing samples for 3 test analytes/compounds.
#' 
#' @name clint_L2
#' @aliases clint_L2
#' @docType data
#' @format A level-2 data.frame with 229 rows and 25 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample description used in the laboratory}
#' \item{\code{Date}}{Date the sample was added to the MS analyzer}
#' \item{\code{Compound.Name}}{Name of the test analyte/compound}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of Clint sample}
#' \item{\code{Dilution.Factor}}{Number of times the sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{ISTD.Name}}{Name of the internal standard (ISTD) analyte/compound}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (pixels)}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level-0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level-0 sample data was extracted}
#' \item{\code{Time}}{Time when the sample was measured - in hours (h)}
#' \item{\code{Test.Compound.Conc}}{Measured concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{Test.Nominal.Conc}}{Expected initial concentration of chemical added to well (uM)}
#' \item{\code{Hep.Density}}{The density (units of millions of hepatocytes per mL) hepatocytes in the \emph{in vitro} incubation}
#' \item{\code{Biological.Replicates}}{Identifier for measurements of multiple samples with the same analyte}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' \item{\code{Verified}}{If "Y", then sample is included in the analysis. (Any other value causes the data to be ignored.)}
#' }
#' 
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#'
"clint_L2"

#' Clint Level-3 Example Data set
#'
#' Mass Spectrometry measurements of intrinsic hepatic clearance (Clint) for cryopreserved 
#' pooled human hepatocytes. Chemicals were per- and poly-fluorinated alkyl substance
#' (PFAS) samples. The experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. This data set is a subset of 
#' experimental data containing samples for 3 test analytes/compounds.
#' 
#' @name clint_L3
#' @aliases clint_L3
#' @docType data
#' @format A level-3 data.frame with 3 rows and 13 variables: \describe{
#' \item{\code{Compound.Name}}{Name of the test analyte/compound}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{Clint}}{Intrinsic hepatic clearance}
#' \item{\code{Clint.pValue}}{p-value of estimated \code{Clint} value}
#' \item{\code{Fit}}{Test nominal concentrations}
#' \item{\code{AIC}}{Akaike Information Criterion of the linear regression fit}
#' \item{\code{AIC.Null}}{Akaike Information Criterion of the exponential decay assuming a constant rate of decay}
#' \item{\code{Clint.1}}{Intrinsic hepatic clearance at 1 uM}
#' \item{\code{Clint.10}}{Intrinsinc hepatic clearance at 10 uM}
#' \item{\code{AIC.Sat}}{Akaike Information Criterion of the exponential decay with a saturation probability}
#' \item{\code{Sat.pValue}}{p-value of exponential decay with a saturation probability}
#' }
#' 
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#'
"clint_L3"

#' Clint Level-4 Example Data set
#' 
#' Mass Spectrometry measurements of intrinsic hepatic clearance (Clint) for cryopreserved 
#' pooled human hepatocytes. Chemicals were per- and poly-fluorinated alkyl substance
#' (PFAS) samples. The experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. This data set is a subset of 
#' experimental data containing samples for 3 test analytes/compounds.
#' 
#' @name clint_L4
#' @aliases clint_L4
#' @docType data
#' @format A level-4 data.frame with 3 rows and 12 variables: \describe{
#' \item{\code{Compound.Name}}{Name of the test analyte/compound}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Clint.1.Med}}{Median intrinsic hepatic clearance at 1 uM}
#' \item{\code{Clint.1.Low}}{2.5th quantile of intrinsic hepatic clearance at 1 uM}
#' \item{\code{Clint.1.High}}{97.5th quantile of intrinsic hepatic clearance at 1 uM}
#' \item{\code{Clint.10.Med}}{Median of intrinsic hepatic clearance at 10 uM}
#' \item{\code{Clint.10.Low}}{2.5th quantile of intrinsic hepatic clearance at 10 uM}
#' \item{\code{Clint.10.High}}{97.5th quantile of intrinsic hepatic clearance at 1 uM}
#' \item{\code{Clint.pValue}}{Probability that a decrease is observed}
#' \item{\code{Sat.pValue}}{Saturation probability that a lower \code{Clint} is observed at a higher concentration}
#' \item{\code{degrades.pValue}}{Probability of abiotic degradation}
#' }
#' 
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#'
"clint_L4"

#' Clint Level-4 PREJAGS arguments
#' 
#' The arguments given to JAGS for the tested compound during level-4 processing of mass spectrometry
#' measurements of intrinsic hepatic clearance (Clint) for cryopreserved 
#' pooled human hepatocytes. Chemicals were per- and poly-fluorinated alkyl substance
#' (PFAS) samples. The experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. This list is overwritten for each tested
#' compound. Therefore, only contains arguments given to JAGS for the last tested compound. 
#' 
#' @name clint_PREJAGS
#' @aliases clint_PREJAGS
#' @docType data
#' @format A named list with 26 elements: \describe{
#' \item{\code{obs}}{\code{Response} of the "Cvst" sample types for the tested compound}
#' \item{\code{Test.Nominal.Conc}}{Unique \code{Test.Nominal.Conc} values (expected initial concentration) of "Cvst" sample types}
#' \item{\code{Num.cal}}{Unique number of \code{Calibration} values}
#' \item{\code{Num.obs}}{Number of \code{Response} of the "Cvst" sample types for the tested compound}
#' \item{\code{obs.conc}}{Indices of the \code{Test.Nominal.Conc} values that corresponds to the "Cvst" sample types' \code{Test.Nominal.Conc}}
#' \item{\code{obs.time}}{\code{Time} of the "Cvst" sample types for the tested compound}
#' \item{\code{obs.cal}}{Indices of the unique "Cvst" \code{Calibration} values that corresponds to the "Cvst" sample types' \code{Calibration}}
#' \item{\code{obs.Dilution.Factor}}{\code{Dilution Factor} of the "Cvst" sample types for the tested compound (number of times the sample was diluted)}
#' \item{\code{Num.blank.obs}}{Number of "Blank" sample types for the tested compound}
#' \item{\code{Blank.obs}}{\code{Response} of the "Blank" sample types for the tested compound}
#' \item{\code{Blank.cal}}{Indices of the unique "Blank" \code{Calibration} values that corresponds to the "Blank" sample types' \code{Calibration}}
#' \item{\code{Blank.Dilution.Factor}}{\code{Dilution Factor} of the "Blank" sample types for the tested compound (number of times the sample was diluted)}
#' \item{\code{Num.cc}}{Number of "CC" sample types with non-NA \code{Test.Compound.Conc} values for the tested compound}
#' \item{\code{cc.obs.conc}}{\code{Test.Compound.Conc} (non-NA) of the "CC" sample types for the tested compound}
#' \item{\code{cc.obs}}{\code{Response} of the "CC" sample types with non-NA \code{Test.Compound.Conc} for the tested compound}
#' \item{\code{cc.obs.cal}}{Indices of the unique "CC" \code{Calibration} values that corresponds to the "CC" sample types' \code{Calibration}}
#' \item{\code{cc.obs.Dilution.Factor}}{\code{Dilution Factor} of the "CC" sample types (number of times the sample was diluted) with non-NA \code{Test.Compound.Conc} for the tested compound}
#' \item{\code{Num.abio.obs}}{Number of "Inactive" samples types for the tested compound}
#' \item{\code{abio.obs}}{\code{Response} of the "Inactive" sample types for the tested compound}
#' \item{\code{abio.obs.conc}}{Indices of the \code{Test.Nominal.Conc} values that corresponds to the "Inactive" sample types' \code{Test.Nominal.Conc}}
#' \item{\code{abio.obs.time}}{\code{Time} of the "Inactive" sample types for the tested compound}
#' \item{\code{abio.obs.cal}}{Indices of the unique "Inactive" \code{Calibration} values that corresponds to the "Inactive" sample types' \code{Calibration}}
#' \item{\code{abio.obs.Dilution.Factor}}{\code{Dilution Factor} of the "Inactive" sample types for the tested compound (number of times the sample was diluted)}
#' \item{\code{DECREASE.PROB}}{Prior probability that a chemical will decrease in the assay. (Defaults to 0.5.)}
#' \item{\code{SATURATE.PROB}}{Prior probability that a chemicals rate of metabolism will decrease between 1 and 10 uM. (Defaults to 0.25.)}
#' \item{\code{DEGRADE.PROB}}{Prior probability that a chemical will be unstable (degrade abiotically) in the assay. (Defaults to 0.05.)}
#' }
#' 
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#'
"clint_PREJAGS"

#' Clint Level-2 Heldout Example Data set
#' 
#' The unverified level-2 samples from mass spectrometry measurements of intrinsic hepatic clearance (Clint) for cryopreserved 
#' pooled human hepatocytes. Chemicals were per- and poly-fluorinated alkyl substance
#' (PFAS) samples. The experiments were led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. This data set is a subset of 
#' experimental data containing samples for 2 test analytes/compounds.
#' 
#' @name clint_L2_heldout
#' @aliases clint_L2_heldout
#' @docType data
#' @format A level-2 data.frame with 10 rows and 25 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample description used in the laboratory}
#' \item{\code{Date}}{Date the sample was added to the MS analyzer}
#' \item{\code{Compound.Name}}{Name of the test analyte/compound}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard - CCD)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of Clint sample}
#' \item{\code{Dilution.Factor}}{Number of times the sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{ISTD.Name}}{Name of the internal standard (ISTD) analyte/compound}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (pixels)}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level-0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level-0 sample data was extracted}
#' \item{\code{Time}}{Time when the sample was measured - in hours (h)}
#' \item{\code{Test.Compound.Conc}}{Measured concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{Test.Nominal.Conc}}{Expected initial concentration of chemical added to well (uM)}
#' \item{\code{Hep.Density}}{The density (units of millions of hepatocytes per mL) hepatocytes in the \emph{in vitro} incubation}
#' \item{\code{Biological.Replicates}}{Identifier for measurements of multiple samples with the same analyte}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' \item{\code{Verified}}{If "Y", then sample is included in the analysis. (Any other value causes the data to be ignored.)}
#' }
#' 
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#'
"clint_L2_heldout"

#' Common Columns in Level-1
#' 
#' Common column names across the various \emph{in vitro} assays used for collecting
#' \emph{in vitro} toxicokinetic parameters.
#' 
#' @name L1.common.col
#' @aliases L1.common.col
#' @docType data
#' @format A named character vector containing the default/standard column names
#' across HTTK assays, where the element names are the corresponding L1 arguments.
"L1.common.cols"

#' Standard Data Catalog (Data Guide) Columns
#' 
#' Standardized column names for data catalogs (i.e. data guides) used for
#' collecting the minimum information to merge level-0 data files.
#' 
#' @name std.catcols
#' @aliases std.catcols
#' @docType data
#' @format A named character vector containing the default/standard column names
#' for data catalogs, where the element names are the corresponding `create_catalog`
#' arguments. 
"std.catcols"