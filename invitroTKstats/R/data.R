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
#' \item{\code{ISTDResponseRatio}}{The ratio between \code{Area} (area for the tested compound) to \code{ISTD.Area} (area for the internal standard).}
#' \item{\code{Time..mins.}}{Time in minues between when the incubation started and when the sample was taken.}
#' \item{\code{Conc}}{Concentration used in the assay.}
#' \item{\code{FileName}}{Name of the level 0 file that contained the data.}
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
#' \item{\code{QSAR_READY_SMILES}}{Short description of the chemical structure.}
#' \item{\code{LC}}{Logical variable if liquid chromatography-mass spectrometry was used.}
#' \item{\code{Agilent.QQQ}}{A column of yes/no/blank factor identifying whether the compound could/could not/wasn't measured with this instrument.}
#' \item{\code{Water.s.Xevo}}{A column of yes/no/blank factor identifying whether the compound could/could not/wasn't measured with this instrument.}
#' \item{\code{AB.Sciex.Qtrap}}{A column of yes/no/blank factor identifying whether the compound could/could not/wasn't measured with this instrument.}
#' \item{\code{GC}}{Logical variable if gas chromatography-mass spectrometry was used.}
#' \item{\code{Agilent.GCMS}}{A column of yes/no/blank factor identifying whether the compound could/could not/wasn't measured with this instrument.}
#' \item{\code{GCTOF}}{A column of yes/no/blank factor identifying whether the compound could/could not/wasn't measured with this instrument.}
#' \item{\code{Comment}}{Any additional comments for a test compound.}
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
#' \item{\code{ISTDResponseRatio}}{The ratio between \code{Area} (area for the tested compound) to \code{ISTD.Area} (area for the internal standard).}
#' \item{\code{Dilution.Factor}}{How many times a sample was diluted.}
#' \item{\code{Protein}}{The percent of physiological protein concentration used for the assay.}
#' \item{\code{RawDataSet}}{Name of the level 0 file that contained the data.}
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
#' @format A level 2 data.frame with 10133 rows 23 variables: \describe{
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
#' @format A level 2 data.frame with 3897 rows 25 variables: \describe{
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
#' @format A level 2 data.frame with 625 rows 24 variables: \describe{
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
#' @format A level 2 data.frame with 5800 rows 24 variables: \describe{
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

