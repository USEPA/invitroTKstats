#This file is used by roxygen2 to generate man files (documentation) for data
#sets included in the package.

#' Mass Spectrometry for \insertCite{wambaugh2019assessing;textual}{invitroTKstats} Hepatocyte Incubations
#'
#' \insertCite{wambaugh2019assessing;textual}{invitroTKstats} includes measurements for intrinsic hepatic clearance
#' Clint measured using in vitro suspensions of pooled primary human hepatocytes
#' \insertCite{shibata2002prediction}{invitroTKstats}.
#'
#' @name wambaugh2019.clint
#' @aliases wambaugh2019.clint
#' @docType data
#' @format A data.frame 23,021 rows and 15 variables: \describe{
#' \item{\code{Preferred.Name}}{Preferred compound name from the CompTox Chemicals Dashboard (CCD).}
#' \item{\code{CAS}}{CAS Registry Number of the test compound.}
#' \item{\code{DTXSID}}{EPA's DSSTox Structure ID.}
#' \item{\code{Sample.Name}}{Sample name used by the laboratory.}
#' \item{\code{Name}}{Compound names used by the laboratory.}
#' \item{\code{Area}}{Mass spectrometry peak area for the test compound.}
#' \item{\code{ISTD.Area}}{Mass spectrometry peak area for the internal standard.}
#' \item{\code{ISTDResponseRatio}}{The ratio between the \code{Area} (area for the tested compound) and the \code{ISTD.Area} (area for the internal standard).}
#' \item{\code{Time..mins.}}{Time in minutes between the start of incubation and when the sample was taken.}
#' \item{\code{Conc}}{Concentration used in the assay.}
#' \item{\code{FileName}}{Name of the level 0 file containing the data.}
#' }
#'
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' \insertRef{wambaugh2019assessing}{invitroTKstats}
#' @source Wambaugh et al. (2019)
#' @keywords data
"wambaugh2019.clint"

#' Mass Spectrometry Methods for \insertCite{wambaugh2019assessing;textual}{invitroTKstats} Chemicals
#'
#' As reported in \insertCite{wambaugh2019assessing;textual}{invitroTKstats}, contractor Cyprotex performed in vitro
#' toxicokinetic experiments requiring the development of analytical chemistry
#' methods to quantitative concentration. 
#' Liquid chromatography/mass spectrometry (LC/MS) was faster and 
#' attempted first. Gas chromatography/mass spectrometry (GC/MS) was attempted
#' when LC/MS failed. For some chemicals the attempts to develop an
#' analysis method failed for both LC/MS and GC/MS and no
#' in vitro experiments were conducted.
#'
#' @name wambaugh2019.methods
#' @aliases Wambaugh2019.methods
#' @docType data
#' @format A data.frame 520 rows and 14 variables: \describe{
#' \item{\code{DTXSID}}{EPA's DSSTox Structure ID.}
#' \item{\code{PREFERRED_NAME}}{Preferred compound name from the CompTox Chemicals Dashboard (CCD).}
#' \item{\code{CASRN}}{CAS Registry Number of the test compound.}
#' \item{\code{MOLECULAR_FORMULA}}{Molecular formula of the test compound.}
#' \item{\code{AVERAGE_MASS}}{Molecular weight of the compound in daltons.}
#' \item{\code{QSAR_READY_SMILES}}{SMILES (Simplified molecular-input line-entry system) chemical structure description.}
#' \item{\code{LC}}{Logical variable indicating whether liquid chromatography-mass spectrometry was used.}
#' \item{\code{Agilent.QQQ}}{A column of yes/no/blank factor identifying whether the compound could/could not/wasn't measured with the Agilent.QQQ instrument.}
#' \item{\code{Water.s.Xevo}}{A column of yes/no/blank factor identifying whether the compound could/could not/wasn't measured with the Water.s.Xevo instrument.}
#' \item{\code{AB.Sciex.Qtrap}}{A column of yes/no/blank factor identifying whether the compound could/could not/wasn't measured with the AB.Sciex.Qtrap instrument.}
#' \item{\code{GC}}{Logical variable indicating whether gas chromatography-mass spectrometry was used.}
#' \item{\code{Agilent.GCMS}}{A column of yes/no/blank factor identifying whether the compound could/could not/wasn't measured with the Agilent.GCMS instrument.}
#' \item{\code{GCTOF}}{A column of yes/no/blank factor identifying whether the compound could/could not/wasn't measured with the GCTOF instrument.}
#' \item{\code{Comment}}{Any additional comments for the test compound.}
#'   }
#' @references
#' \insertRef{wambaugh2019assessing}{invitroTKstats}
#' @source Wambaugh et al. (2019)
#' @keywords data
"wambaugh2019.methods"

#' Mass Spectrometry for \insertCite{wambaugh2019assessing;textual}{invitroTKstats} Protein Binding
#'
#' \insertCite{wambaugh2019assessing;textual}{invitroTKstats} includes chemical-specific measurements of plasma
#' protein binding using the method of rapid equilibriuim dialysis
#' \insertCite{waters2008validation}{invitroTKstats}.
#'
#' @name wambaugh2019.red
#' @aliases Wambaugh2019.red
#' @docType data
#' @format A data.frame 17,689 rows and 14 variables: \describe{
#' \item{\code{Preferred.Name}}{Preferred compound name from the CompTox Chemicals Dashboard (CCD).}
#' \item{\code{CAS}}{CAS Registry Number of the test compound.}
#' \item{\code{DTXSID}}{EPA's DSSTox Structure ID.}
#' \item{\code{SampleName}}{Sample name used by the laboratory.}
#' \item{\code{CompoundName}}{Compound name used by the laboratory.}
#' \item{\code{Area}}{Mass spectrometry peak area for the test compound.}
#' \item{\code{ISTD.Area}}{Mass spectrometry peak area for the internal standard.}
#' \item{\code{ISTDResponseRatio}}{The ratio between the \code{Area} (area for the tested compound) and the \code{ISTD.Area} (area for the internal standard).}
#' \item{\code{Dilution.Factor}}{The number of times a sample was diluted.}
#' \item{\code{Protein}}{The percent of physiological protein concentration used for the assay.}
#' \item{\code{RawDataSet}}{Name of the level 0 file containing the data.}
#' }
#'
#' @references
#' \insertRef{waters2008validation}{invitroTKstats}
#'
#' \insertRef{wambaugh2019assessing}{invitroTKstats}
#' @source Wambaugh et al. (2019)
#' @keywords data
"wambaugh2019.red"

#' \insertCite{smeltz2023plasma;textual}{invitroTKstats} Ultracentrifugation Data Set
#'
#' Mass Spectrometry measurements of plasma protein binding measured by
#' ultracentrifugation for per- and poly-fluorinated alkyl substance
#' (PFAS) samples from experiments led by Dr.s Marci Smeltz and Barbara Wetmore.
#'
#' @name smeltz2023.uc
#' @aliases Smeltz2023.uc
#' @docType data
#' @format A level 2 data.frame with 10133 rows and 23 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample description used in the laboratory}
#' \item{\code{Date}}{Date sample was acquired}
#' \item{\code{Compound.Name}}{Compound name}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of UC sample}
#' \item{\code{Dilution.Factor}}{Number of times sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{Standard.Conc}}{Concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{TC.Assay.T1.Conc}}{Intended concentration of chemical intended in T1 sample (uM)}
#' \item{\code{ISTD.Name}}{Name of compound used as internal standard (ISTD)}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (pixels)}
#' \item{\code{Series}}{Identier for replicate series of UC measurements}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis)}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample)}
#' \item{\code{Level0.File}}{Name of data file from laboratory that was used to compile level0 data table)}
#' \item{\code{Level0.Sheet}}{Name of "sheet" (for Excel workbooks) from which the laboratory data were read)}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' \item{\code{Verified}}{If ="Y" then this sample is included in the analysis. Any other value leads to the data being ignored.)}
#' }
#'
#' @references
#' \insertRef{howard2010plasma}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#' @keywords data
"smeltz2023.uc"

#' \insertCite{smeltz2023plasma;textual}{invitroTKstats} rapid equlibriation dialysis data set
#'
#' Mass Spectrometry measurements of plasma protein binding measured by
#' rapid equlibriation dialysis for per- and poly-fluorinated alkyl substance
#' (PFAS) samples from experiments led by Dr.s Marci Smeltz and Barbara Wetmore.
#'
#' @name smeltz2023.red
#' @aliases Smeltz2023.red
#' @docType data
#' @format A level 2 data.frame with 3897 rows and 25 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample description used in the laboratory}
#' \item{\code{Date}}{Date sample was acquired}
#' \item{\code{Compound.Name}}{Compound name}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of RED sample}
#' \item{\code{Dilution.Factor}}{Number of times sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{Standard.Conc}}{Concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{Nominal.Test.Conc}}{Intended concentration of chemical introduced into RED plate (uM)}
#' \item{\code{Percent.Physiologic.Plasma}}{Percent of physiological plasma concentration in RED plate (in percent)}
#' \item{\code{Time}}{Time of sample measurment (h)}
#' \item{\code{ISTD.Name}}{Name of compound used as internal standard (ISTD)}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (pixels)}
#' \item{\code{Series}}{Identier for replicate series of UC measurements}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis)}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample)}
#' \item{\code{Level0.File}}{Name of data file from laboratory that was used to compile level0 data table)}
#' \item{\code{Level0.Sheet}}{Name of "sheet" (for Excel workbooks) from which the laboratory data were read)}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' \item{\code{Verified}}{If ="Y" then this sample is included in the analysis. Any other value leads to the data being ignored.)}
#' }
#'
#' @references
#' \insertRef{waters2008validation}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#' @keywords data
"smeltz2023.red"

#' \insertCite{smeltz2023plasma;textual}{invitroTKstats} Intinsic hepatic clearnace data set
#'
#' Mass Spectrometry measurements of intrinsic hepatic clearance for
#' cryopreserved pooled human hepatocytes. Chemicals were per- and
#' 'poly-fluorinated alkyl substance
#' (PFAS) samples. The experiments were
#' led by Dr.s Marci Smeltz and Barbara Wetmore.
#'
#' @name smeltz2023.clint
#' @aliases Smeltz2023.clint
#' @docType data
#' @format A level 2 data.frame with 625 rows and 24 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample description used in the laboratory}
#' \item{\code{Date}}{Date sample was acquired}
#' \item{\code{Compound.Name}}{Compound name}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of UC sample}
#' \item{\code{Dilution.Factor}}{Number of times sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{Standard.Conc}}{Concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{Clint.Assay.Conc}}{Intended initial concentration of chemical (uM)}
#' \item{\code{Time}}{Time point sample measured (h)}
#' \item{\code{ISTD.Name}}{Name of compound used as internal standard (ISTD)}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (pixels)}
#' \item{\code{Series}}{Identier for replicate series of UC measurements}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis)}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample)}
#' \item{\code{Level0.File}}{Name of data file from laboratory that was used to compile level0 data table)}
#' \item{\code{Level0.Sheet}}{Name of "sheet" (for Excel workbooks) from which the laboratory data were read)}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' \item{\code{Verified}}{If ="Y" then this sample is included in the analysis. Any other value leads to the data being ignored.)}
#' }
#'
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#' @keywords data
"smeltz2023.clint"

#' \insertCite{kreutz2023category;textual}{invitroTKstats} Ultracentrifugation Data Set
#'
#' Mass Spectrometry measurements of plasma protein binding measured by
#' ultracentrifugation for per- and poly-fluorinated alkyl substance
#' (PFAS) samples from experiments led by Dr.s Anna Kreutz and Barbara Wetmore.
#'
#' @name kreutz2023
#' @aliases Kreutz2023
#' @docType data
#' @format A data.frame with 2928 rows and 23 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample description used in the laboratory}
#' \item{\code{Date}}{Date sample was acquired}
#' \item{\code{Compound.Name}}{Compound name}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of UC sample}
#' \item{\code{Dilution.Factor}}{Number of times sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{Standard.Conc}}{Concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{TC.Assay.T1.Conc}}{Intended concentration of chemical intended in T1 sample (uM)}
#' \item{\code{ISTD.Name}}{Name of compound used as internal standard (ISTD)}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (pixels)}
#' \item{\code{Series}}{Identier for replicate series of UC measurements}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis)}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample)}
#' \item{\code{Level0.File}}{Name of data file from laboratory that was used to compile level0 data table)}
#' \item{\code{Level0.Sheet}}{Name of "sheet" (for Excel workbooks) from which the laboratory data were read)}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' \item{\code{Verified}}{If ="Y" then this sample is included in the analysis. Any other value leads to the data being ignored.)}
#' }
#' @references
#' \insertRef{howard2010plasma}{invitroTKstats}
#'
#' \insertRef{kreutz2023category}{invitroTKstats}
#' @keywords data
"kreutz2023.uc"

#' \insertCite{kreutz2023category;textual}{invitroTKstats} intrinsic hepatic clearance data set
#'
#' Mass Spectrometry measurements of intrinsic hepatic clearance for
#' cryopreserved pooled human hepatocytes. Chemicals were per- and
#' 'poly-fluorinated alkyl substance
#' (PFAS) samples. The experiments were
#' led by Dr.s Anna Kreutz and Barbara Wetmore.
#'
#' @name kreutz2023.clint
#' @aliases Kreutz2023.clint
#' @docType data
#' @format A level 2 data.frame with 5800 rows and 24 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample description used in the laboratory}
#' \item{\code{Date}}{Date sample was acquired}
#' \item{\code{Compound.Name}}{Compound name}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of UC sample}
#' \item{\code{Dilution.Factor}}{Number of times sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{Standard.Conc}}{Concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{Clint.Assay.Conc}}{Intended initial concentration of chemical (uM)}
#' \item{\code{Time}}{Time point sample measured (h)}
#' \item{\code{ISTD.Name}}{Name of compound used as internal standard (ISTD)}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (pixels)}
#' \item{\code{Series}}{Identier for replicate series of UC measurements}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis)}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample)}
#' \item{\code{Level0.File}}{Name of data file from laboratory that was used to compile level0 data table)}
#' \item{\code{Level0.Sheet}}{Name of "sheet" (for Excel workbooks) from which the laboratory data were read)}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' \item{\code{Verified}}{If ="Y" then this sample is included in the analysis. Any other value leads to the data being ignored.)}
#' }
#'
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' \insertRef{kreutz2023category}{invitroTKstats}
"kreutz2023.clint"

#' Caco2 Level 0 Example Data set
#'
#' A subset of mass spectrometry (MS) measurements of membrane permeability from Caco2 cells. 
#' This subset contains samples for 3 test analytes/compounds.
#' 
#' @name caco2_L0
#' @aliases caco2_L0
#' @docType data
#' @format A level 0 data.frame with 48 rows and 17 variables: \describe{
#' \item{\code{Compound}}{Compound name}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard)}
#' \item{\code{Lab.Compound.ID}}{Compound ID used in the laboratory}
#' \item{\code{Date}}{Date sample was acquired}
#' \item{\code{Sample}}{Sample Name}
#' \item{\code{Type}}{Type of Caco2 sample}
#' \item{\code{Compound.Conc}}{Concentration of test chemical (for calibration curve) (uM)}
#' \item{\code{Peak.Area}}{Peak area of analyte (target compound)}
#' \item{\code{ISTD.Peak.Area}}{Peak area of internal standard (pixels)}
#' \item{\code{ISTD.Name}}{Name of compound used as internal standard (ISTD)}
#' \item{\code{Analysis.Params}}{General description of chemical analysis method}
#' \item{\code{Level0.File}}{Name of data file from laboratory that was used to compile level0 data table)}
#' \item{\code{Level0.Sheet}}{Name of "sheet" (for Excel workbooks) from which the laboratory data were read)}
#' \item{\code{Direction}}{Direction of the Caco-2 permeability experiment}
#' \item{\code{Vol.Donor}}{The volume (in cm^3) of the donor portion of the Caco-2 experimental well}
#' \item{\code{Vol.Receiver}}{The volume (in cm^3) of the receiver portion of the Caco-2 experimental well}
#' \item{\code{Dilution.Factor}}{Number of times sample was diluted}
#' \item{\code{Note}}{Additional information}
#' }
#'
"caco2_L0"

#' Caco2 Level 1 Example Data set
#'
#' A subset of mass spectrometry (MS) measurements of membrane permeability from Caco2 cells. 
#' This subset contains samples for 3 test analytes/compounds.
#' 
#' @name caco2_L1
#' @aliases caco2_L1
#' @docType data
#' @format A level 1 data.frame with 48 rows and 24 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample name as described in the laboratory}
#' \item{\code{Date}}{Date sample was acquired}
#' \item{\code{Compound.Name}}{Compound name}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of Caco2 sample}
#' \item{\code{Direction}}{Direction of the Caco-2 permeability experiment}
#' \item{\code{Dilution.Factor}}{Number of times sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{Biological.Replicates}}{Replicate series of measurements with the same analyte}
#' \item{\code{Nominal.Conc}}{Concentration of test chemical (for calibration curve) (uM)}
#' \item{\code{Test.Target.Conc}}{Test concentration of chemical (uM) at time zero}
#' \item{\code{Time}}{Time point sample measured (h)}
#' \item{\code{ISTD.Name}}{Name of compound used as internal standard (ISTD)}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (pixels)}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Membrane.Area}}{The area of the Caco-2 monolayer.}
#' \item{\code{Vol.Donor}}{The volume (in cm^3) of the donor portion of the Caco-2 experimental well}
#' \item{\code{Vol.Receiver}}{The volume (in cm^3) of the receiver portion of the Caco-2 experimental well}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis)}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Additional information}
#' \item{\code{Level0.File}}{Name of data file from laboratory that was used to compile level0 data table)}
#' \item{\code{Level0.Sheet}}{Name of "sheet" (for Excel workbooks) from which the laboratory data were read)}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' }
#'
"caco2_L1"

#' Caco2 Level 2 Example Data set
#'
#' A subset of mass spectrometry (MS) measurements of membrane permeability from Caco2 cells. 
#' This subset contains samples for 3 test analytes/compounds.
#' 
#' @name caco2_L2
#' @aliases caco2_L2
#' @docType data
#' @format A level 2 data.frame with 48 rows and 25 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample name as described in the laboratory}
#' \item{\code{Date}}{Date sample was acquired}
#' \item{\code{Compound.Name}}{Compound name}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of Caco2 sample}
#' \item{\code{Direction}}{Direction of the Caco-2 permeability experiment}
#' \item{\code{Dilution.Factor}}{Number of times sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{Biological.Replicates}}{Replicate series of measurements with the same analyte}
#' \item{\code{Nominal.Conc}}{Concentration of test chemical (for calibration curve) (uM)}
#' \item{\code{Test.Target.Conc}}{Test concentration of chemical (uM) at time zero}
#' \item{\code{Time}}{Time point sample measured (h)}
#' \item{\code{ISTD.Name}}{Name of compound used as internal standard (ISTD)}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (pixels)}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Membrane.Area}}{The area of the Caco-2 monolayer.}
#' \item{\code{Vol.Donor}}{The volume (in cm^3) of the donor portion of the Caco-2 experimental well}
#' \item{\code{Vol.Receiver}}{The volume (in cm^3) of the receiver portion of the Caco-2 experimental well}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis)}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Additional information}
#' \item{\code{Level0.File}}{Name of data file from laboratory that was used to compile level0 data table)}
#' \item{\code{Level0.Sheet}}{Name of "sheet" (for Excel workbooks) from which the laboratory data were read)}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' \item{\code{Verified}}{If ="Y" then this sample is included in the analysis. Any other value leads to the data being ignored.)}
#' }
#'
"caco2_L2"

#' Fup-UC Level 0 Example Data set
#' 
#' Mass Spectrometry measurements of plasma protein binding (PPB) via 
#' ultracentrifugation (UC) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples from experiments led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. 
#' This data set is a subset of experimental data containing samples for 
#' 3 test analytes/compounds.
#'
#' @name fup_uc_L0
#' @aliases fup_uc_L0
#' @docType data
#' @format A level 0 data.frame with 240 rows and 17 variables: \describe{
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
#' \item{\code{Analysis.Params}}{The column contains the retention time}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level0 sample data was extracted}
#' \item{\code{Sample.Text}}{Additional notes on the sample}
#' \item{\code{Sample.Type}}{Type of UC sample in the package's annotations}
#' \item{\code{Dilution.Factor}}{Number of times the sample was diluted}
#' \item{\code{Replicate}}{Identifier for parallel measurements of multiple samples of a compound}
#' }
#'
#' @references
#' \insertRef{howard2010plasma}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#' @keywords data
"fup_uc_L0"

#' Fup-UC Level 1 Example Data set
#' 
#' Mass Spectrometry measurements of plasma protein binding (PPB) via 
#' ultracentrifugation (UC) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples from experiments led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. 
#' This data set is a subset of experimental data containing samples for 
#' 3 test analytes/compounds.
#'
#' @name fup_uc_L1
#' @aliases fup_uc_L1
#' @docType data
#' @format A level 1 data.frame with 240 rows and 23 variables: \describe{
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
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis)}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample)}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level0 sample data was extracted}
#' \item{\code{Test.Compound.Conc}}{Expected (or nominal) concentration of analyte (for calibration curve)}
#' \item{\code{UC.Assay.T1.Conc}}{Intended concentration of chemical intended in T1 sample (uM)}
#' \item{\code{Biological.Replicates}}{Identifier for parallel measurements of multiple samples of a compound}
#' \item{\code{Technical.Replicates}}{Identifier for repeated measurements of a sample of a compound}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' }
#'
#' @references
#' \insertRef{howard2010plasma}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#' @keywords data
"fup_uc_L1"

#' Fup-UC Level 2 Example Data set
#' 
#' Mass Spectrometry measurements of plasma protein binding (PPB) via 
#' ultracentrifugation (UC) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples from experiments led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}. 
#' This data set is a subset of experimental data containing samples for 
#' 3 test analytes/compounds.
#'
#' @name fup_uc_L2
#' @aliases fup_uc_L2
#' @docType data
#' @format A level 2 data.frame with 240 rows and 24 variables: \describe{
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
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis)}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample)}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level0 sample data was extracted}
#' \item{\code{Test.Compound.Conc}}{Expected (or nominal) concentration of analyte (for calibration curve)}
#' \item{\code{UC.Assay.T1.Conc}}{Intended concentration of chemical intended in T1 sample (uM)}
#' \item{\code{Biological.Replicates}}{Identifier for parallel measurements of multiple samples of a compound}
#' \item{\code{Technical.Replicates}}{Identifier for repeated measurements of a sample of a compound}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' \item{\code{Verified}}{If ="Y" then this sample is included in the analysis. Any other value leads to the data being ignored.)}
#' }
#'
#' @references
#' \insertRef{howard2010plasma}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#' @keywords data
"fup_uc_L2"

#' Fup RED Level 0 Example Data set
#'
#' Mass Spectrometry measurements of plasma protein binding (PPB) via rapid 
#' equilibrium dialysis (RED) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples from experiments led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}.
#' This data set is a subset of experimental data containing samples for 
#' 3 test analytes/compounds.
#' 
#' @name fup_red_L0
#' @aliases fup_red_L0
#' @docType data
#' @format A level 0 data.frame with 660 rows and 18 variables: \describe{
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
#' \item{\code{Analysis.Params}}{The column contains the retention time}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level0 sample data was extracted}
#' \item{\code{Sample Text}}{Additional notes on the sample}
#' \item{\code{Sample.Type}}{Type of RED sample in the package's annotations}
#' \item{\code{Replicate}}{Identifier for parallel measurements of multiple samples of a compound}
#' \item{\code{Time}}{Time point the sample was measured - in hours (h)}
#' \item{\code{Dilution.Factor}}{Number of times the sample was diluted}
#' }
#' 
#' @references
#' \insertRef{waters2008validation}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#'
"fup_red_L0"

#' Fup RED Level 1 Example Data set
#'
#' Mass Spectrometry measurements of plasma protein binding (PPB) via rapid 
#' equilibrium dialysis (RED) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples from experiments led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}.
#' This data set is a subset of experimental data containing samples for 
#' 3 test analytes/compounds.
#' 
#' @name fup_red_L1
#' @aliases fup_red_L1
#' @docType data
#' @format A level 1 data.frame with 492 rows and 24 variables: \describe{
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
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis)}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample)}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level0 sample data was extracted}
#' \item{\code{Time}}{Time point the sample was measured - in hours (h)}
#' \item{\code{Test.Compound.Conc}}{Expected (or nominal) concentration of analyte (for calibration curve)}
#' \item{\code{Test.Nominal.Conc}}{Intended concentration of chemical introduced into RED plate (uM)}
#' \item{\code{Percent.Physiologic.Plasma}}{Percent of physiological plasma concentration in RED plate (in percent)}
#' \item{\code{Technical.Replicates}}{Identifier for repeated measurements of a sample of a compound}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' }
#'
#' @references
#' \insertRef{howard2010plasma}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#' @keywords data
"fup_red_L1"

#' Fup RED Level 2 Example Data set
#'
#' Mass Spectrometry measurements of plasma protein binding (PPB) via rapid 
#' equilibrium dialysis (RED) for per- and poly-fluorinated alkyl substance
#' (PFAS) samples from experiments led by Dr.s Marci Smeltz and Barbara Wetmore 
#' \insertCite{smeltz2023plasma}{invitroTKstats}.
#' This data set is a subset of experimental data containing samples for 
#' 3 test analytes/compounds.
#' 
#' @name fup_red_L2
#' @aliases fup_red_L2
#' @docType data
#' @format A level 2 data.frame with 492 rows and 25 variables: \describe{
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
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis)}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample)}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level0 sample data was extracted}
#' \item{\code{Time}}{Time point the sample was measured - in hours (h)}
#' \item{\code{Test.Compound.Conc}}{Expected (or nominal) concentration of analyte (for calibration curve)}
#' \item{\code{Test.Nominal.Conc}}{Intended concentration of chemical introduced into RED plate (uM)}
#' \item{\code{Percent.Physiologic.Plasma}}{Percent of physiological plasma concentration in RED plate (in percent)}
#' \item{\code{Technical.Replicates}}{Identifier for repeated measurements of a sample of a compound}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' \item{\code{Verified}}{If ="Y" then this sample is included in the analysis. Any other value leads to the data being ignored.)}
#' }
#'
#' @references
#' 
#' \insertRef{waters2008validation}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
"fup_red_L2"

#' Clint Level 0 Example Data set
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
#' @format A level 0 data.frame with 247 rows and 16 variables: \describe{
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
#' \item{\code{Analysis.Params}}{The column contains the retention time}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level0 sample data was extracted}
#' \item{\code{Sample.Text}}{Additional notes on the sample}
#' \item{\code{Time}}{Time point the sample was measured - in hours (h)}
#' \item{\code{Dilution.Factor}}{Number of times the sample was diluted}
#' }
#' 
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#'
"clint_L0"

#' Clint Level 1 Example Data set
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
#' @format A level 1 data.frame with 229 rows and 24 variables: \describe{
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
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis)}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample)}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level0 sample data was extracted}
#' \item{\code{Time}}{Time point the sample was measured - in hours (h)}
#' \item{\code{Test.Compound.Conc}}{Expected concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{Clint.Assay.Conc}}{Intended initial concentration of chemical (uM)}
#' \item{\code{Hep.Density}}{The density (units of millions of hepatocytes per mL) hepatocytes in the in vitro incubation}
#' \item{\code{Biological.Replicates}}{Identifier for parallel measurements of multiple samples of a compound}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' }
#' 
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#'
"clint_L1"

#' Clint Level 2 Example Data set
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
#' @format A level 2 data.frame with 229 rows and 25 variables: \describe{
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
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis)}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample)}
#' \item{\code{Level0.File}}{Name of the laboratory data file from which the level0 sample data was extracted}
#' \item{\code{Level0.Sheet}}{Name of the Excel workbook 'sheet' from which the level0 sample data was extracted}
#' \item{\code{Time}}{Time point the sample was measured - in hours (h)}
#' \item{\code{Test.Compound.Conc}}{Expected concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{Clint.Assay.Conc}}{Intended initial concentration of chemical (uM)}
#' \item{\code{Hep.Density}}{The density (units of millions of hepatocytes per mL) hepatocytes in the in vitro incubation}
#' \item{\code{Biological.Replicates}}{Identifier for parallel measurements of multiple samples of a compound}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' \item{\code{Verified}}{If ="Y" then this sample is included in the analysis. Any other value leads to the data being ignored.)}
#' }
#' 
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#'
#' \insertRef{smeltz2023plasma}{invitroTKstats}
#'
"clint_L2"

#' Common Columns in Level-1
#' 
#' Common column names across the various *in vitro* assays used for collecting
#' HTTK relevant physiological parameters.
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