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
#' \item{\code{Sample.Type}}{Type of UC sample}
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
#' \item{\code{Compound.Conc}}{Test concentration (uM)}
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
#' \item{\code{Series}}{Identier for replicate series of measurements}
#' \item{\code{Standard.Conc}}{Concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{Test.Target.Conc}}{Test concentration of chemical (uM)}
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
#' \item{\code{Series}}{Identier for replicate series of measurements}
#' \item{\code{Standard.Conc}}{Concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{Test.Target.Conc}}{Test concentration of chemical (uM)}
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
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' \item{\code{Verified}}{If ="Y" then this sample is included in the analysis. Any other value leads to the data being ignored.)}
#' }
#'
"caco2_L2"

#' Fup-UC Level 0 Example Data set
#' 
#' Mass Spectrometry measurements of plasma protein binding measured by
#' ultracentrifugation for per- and poly-fluorinated alkyl substance
#' (PFAS) samples from experiments led by Dr.s Marci Smeltz and Barbara Wetmore. 
#' This data set is a subset of the experiment data which contains samples for 
#' 3 test analytes/compounds.
#'
#' @name fup_uc_L0
#' @aliases fup_uc_L0
#' @docType data
#' @format A level 0 data.frame with 224 rows and 32 variables: \describe{
#' \item{\code{Name}}{Sample description used in the laboratory}
#' \item{\code{Sample.Text}}{Additional notes on the sample}
#' \item{\code{Type}}{Sample types}
#' \item{\code{RT}}{Retention time}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Height}}{Height of the peak}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (pixels)}
#' \item{\code{Response}}{Response factor (calculated from analyte and ISTD peaks)}
#' \item{\code{Coeff..Of.Determination}}{Coefficient of determination/R sqaured}
#' \item{\code{Std..Conc}}{}
#' \item{\code{uM}}{}
#' \item{\code{X.Dev}}{}
#' \item{\code{Primary.Flags}}{Flags the libratory used to label samples}
#' \item{\code{Acq.Date}}{Date sample was acquired}
#' \item{\code{Date}}{Date sample was acquired}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard)}
#' \item{\code{ISTD.Name}}{Name of compound used as internal standard (ISTD)}
#' \item{\code{Compound.Name}}{Compound Name}
#' \item{\code{Std.Units}}{Unit for standard}
#' \item{\code{Level0.Sheet}}{Name of "sheet" (for Excel workbooks) from which the laboratory data were read)}
#' \item{\code{Sample.Type}}{Type of UC sample in invitroTKstats annotaions}
#' \item{\code{Std.Conc}}{Concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{Series}}{Identier for replicate series of UC measurements}
#' \item{\code{Cal}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{Dilution.Factor}}{Number of times sample was diluted}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis)}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{Test.Target.Conc}}{Test chemical concentration (uM)}
#' \item{\code{Replicate}}{indicator for technical replcaites of UC measurments}
#' \item{\code{Level0.File}}{Name of data file from laboratory that was used to compile level0 data table)}
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
#' Mass Spectrometry measurements of plasma protein binding measured by
#' ultracentrifugation for per- and poly-fluorinated alkyl substance
#' (PFAS) samples from experiments led by Dr.s Marci Smeltz and Barbara Wetmore. 
#' This data set is a subset of the experiment data which contains samples for 
#' 3 test analytes/compounds.
#'
#' @name fup_uc_L1
#' @aliases fup_uc_L1
#' @docType data
#' @format A level 1 data.frame with 224 rows and 22 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample description used in the laboratory}
#' \item{\code{Date}}{Date sample was acquired}
#' \item{\code{Compound.Name}}{Compound name}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of UC sample}
#' \item{\code{Dilution.Factor}}{Number of times sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{ISTD.Name}}{Name of compound used as internal standard (ISTD)}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (pixels)}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis)}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample)}
#' \item{\code{Level0.File}}{Name of data file from laboratory that was used to compile level0 data table)}
#' \item{\code{Level0.Sheet}}{Name of "sheet" (for Excel workbooks) from which the laboratory data were read)}
#' \item{\code{Test.Compound.Conc}}{Concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{UC.Assay.T1.Conc}}{Intended concentration of chemical intended in T1 sample (uM)}
#' \item{\code{Technical.Replicates}}{Identier for technical replicates of UC measurements}
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
#' Mass Spectrometry measurements of plasma protein binding measured by
#' ultracentrifugation for per- and poly-fluorinated alkyl substance
#' (PFAS) samples from experiments led by Dr.s Marci Smeltz and Barbara Wetmore. 
#' This data set is a subset of the experiment data which contains samples for 
#' 3 test analytes/compounds.
#'
#' @name fup_uc_L2
#' @aliases fup_uc_L2
#' @docType data
#' @format A level 1 data.frame with 224 rows and 23 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{Sample description used in the laboratory}
#' \item{\code{Date}}{Date sample was acquired}
#' \item{\code{Compound.Name}}{Compound name}
#' \item{\code{DTXSID}}{DSSTox Substance Identifier (CompTox Chemicals Dashboard)}
#' \item{\code{Lab.Compound.Name}}{Compound as described in the laboratory}
#' \item{\code{Sample.Type}}{Type of UC sample}
#' \item{\code{Dilution.Factor}}{Number of times sample was diluted}
#' \item{\code{Calibration}}{Identifier for mass spectrometry calibration -- usually the date}
#' \item{\code{ISTD.Name}}{Name of compound used as internal standard (ISTD)}
#' \item{\code{ISTD.Conc}}{Concentration of ISTD (uM)}
#' \item{\code{ISTD.Area}}{Peak area of internal standard (pixels)}
#' \item{\code{Area}}{Peak area of analyte (target compound)}
#' \item{\code{Analysis.Method}}{General description of chemical analysis method}
#' \item{\code{Analysis.Instrument}}{Instrument(s) used for chemical analysis)}
#' \item{\code{Analysis.Parameters}}{Parameters for identifing analyte peak (for example, retention time)}
#' \item{\code{Note}}{Any laboratory notes about sample)}
#' \item{\code{Level0.File}}{Name of data file from laboratory that was used to compile level0 data table)}
#' \item{\code{Level0.Sheet}}{Name of "sheet" (for Excel workbooks) from which the laboratory data were read)}
#' \item{\code{Test.Compound.Conc}}{Concentration of analytic standard (for calibration curve) (uM)}
#' \item{\code{UC.Assay.T1.Conc}}{Intended concentration of chemical intended in T1 sample (uM)}
#' \item{\code{Technical.Replicates}}{Identier for technical replicates of UC measurements}
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

