#This file is used by roxygen2 to generate man files (documentation) for data
#sets included in the package.

#' Mass Spectrometry for Wambaugh et al. (2019) Hepatocyte Incubations
#' 
#' Wambaugh et al. (2019) includes measurments for intrinsic hepatic clearance
#' Clint measured using in vitro suspoensions of pooled primary human hepatocytes
#' (Shibata, et al. 2002).)
#' 
#' @name wambaugh2019.clint
#' @aliases wambaugh2019.clint
#' @docType data
#' @format A data.frame 23,021 rows and 15 variables: \describe{
#' \item{\code{Preferred.Name}}{}
#' \item{\code{CAS}}{}
#' \item{\code{DTXSID}}{}
#' \item{\code{Sample.Name}}{}
#' \item{\code{Name}}{}
#' \item{\code{Transition}}{}
#' \item{\code{Area}}{}
#' \item{\code{ISTD.Area}}{}
#' \item{\code{ISTDResponseRatio}}{}
#' \item{\code{ln...Remaining}}{}
#' \item{\code{Time..mins.}}{}
#' \item{\code{Conc}}{}
#' \item{\code{TaskOrder}}{}
#' \item{\code{FileName}}{}
#' \item{\code{X}}{} 
#' } 
#'
#' @references
#' \insertRef{shibata2002prediction}{invitroTKstats}
#' \insertRef{wambaugh2019assessing}{invitroTKstats}
#' @source Wambaugh et al. (2019)
#' @keywords data
"wambaugh2019.clint"   

#' Mass Spectrometry Methods for Wambaugh et al. (2019) Chemicals
#' 
#' As reported in Wambaugh et al. (2019), contractor Cyprotex performed in vitro
#' toxicokinetic experiments requiring the development of analytical chemistry
#' methods to quantitate concentration. 
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
#' \item{\code{DTXSID}}{}
#' \item{\code{PREFERRED_NAME}}{}
#' \item{\code{CASRN}}{}
#' \item{\code{MOLECULAR_FORMULA}}{}
#' \item{\code{AVERAGE_MASS}}{}
#' \item{\code{QSAR_READY_SMILES}}{}
#' \item{\code{LC}}{}
#' \item{\code{Agilent.QQQ}}{}
#' \item{\code{Water.s.Xevo}}{}
#' \item{\code{AB.Sciex.Qtrap}}{}
#' \item{\code{GC}}{}
#' \item{\code{Agilent.GCMS}}{}
#' \item{\code{GCTOF}}{}
#' \item{\code{Comment}}{}
#'   }
#' @references
#' \insertRef{wambaugh2019assessing}{invitroTKstats}
#' @source Wambaugh et al. (2019)
#' @keywords data
"wambaugh2019.methods" 

#' Mass Spectrometry for Wambaugh et al. (2019) Protein Binding
#' 
#' Wambaugh et al. (2019) includes chemical-specific measurments of plasma
#' protein binding using the method of rapid equilibriuim dialysis (RED, 
#' Waters et al. 2008)
#' 
#' @name wambaugh2019.red
#' @aliases Wambaugh2019.red
#' @docType data
#' @format A data.frame 17,689 rows and 14 variables: \describe{
#' \item{\code{Preferred.Name}}{}
#' \item{\code{CAS}}{}
#' \item{\code{DTXSID}}{}             
#' \item{\code{SampleName}}{}
#' \item{\code{CompoundName}}{}
#' \item{\code{Transition}}{}          
#' \item{\code{Area}}{}
#' \item{\code{ISTD.Area}}{}
#' \item{\code{ISTDResponseRatio}}{}  
#' \item{\code{Dilution.Factor}}{}
#' \item{\code{Protein}}{}
#' \item{\code{Task.Order}}{}         
#' \item{\code{RawDataSet}}{}
#' \item{\code{CyprotexEPASetNumber}}{}
#' }
#'
#' @references
#' \insertRef{waters2008validation}{invitroTKstats}
#' \insertRef{wambaugh2019assessing}{invitroTKstats}
#' @source Wambaugh et al. (2019)
#' @keywords data
"wambaugh2019.red"  

#' Smeltz 2020 Ultracentrifucation Data Set 
#' 
#' 
#' Mass Spectrometry measurments of plasma protein binding measured by 
#' ultracentrifucation for per- and poly-fluorinated alkyl subtance
#' (PFAS) samples from experiments conducted by Dr. Marci Smeltz
#' 
#' @name smeltz2020
#' @aliases Smeltz2020
#' @docType data
#' @format A data.frame with 289 rows and 19 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{}
#' \item{\code{Date}}{}
#' \item{\code{Compound.Name}}{}
#' \item{\code{DTXSID}}{}
#' \item{\code{Lab.Compound.Name}}{}
#' \item{\code{Sample.Type}}{}
#' \item{\code{Dilution.Factor}}{}
#' \item{\code{Cal}}{}
#' \item{\code{Standard.Conc}}{}
#' \item{\code{Test.Target.Conc}}{}
#' \item{\code{ISTD.Name}}{}
#' \item{\code{ISTD.Conc}}{}
#' \item{\code{IS.Area}}{}
#' \item{\code{Series}}{}
#' \item{\code{Area}}{}
#' \item{\code{Response}}{}
#' \item{\code{Analysis.Method}}{}
#' \item{\code{Analysis.Instrument}}{}
#' \item{\code{Analysis.Parameters}}{}
#' }
#'
#' @references
#' \insertRef{howard2010plasma}{invitroTKstats}
#' @keywords data
"smeltz2020"   

#' Kreutz 2020 Ultracentrifucation Data Set 
#' 
#' Mass Spectrometry measurments of plasma protein binding measured by 
#' ultracentrifucation for per- and poly-fluorinated alkyl subtance
#' (PFAS) samples from experiments conducted by Dr. Anna Kreutz
#' 
#' @name kreutz2020
#' @aliases Kreutz2020
#' @docType data
#' @format A data.frame with 228 rows and 19 variables: \describe{
#' \item{\code{Lab.Sample.Name}}{}
#' \item{\code{Date}}{}
#' \item{\code{Compound.Name}}{}
#' \item{\code{DTXSID}}{}
#' \item{\code{Lab.Compound.Name}}{}
#' \item{\code{Sample.Type}}{}
#' \item{\code{Dilution.Factor}}{}
#' \item{\code{Cal}}{}
#' \item{\code{Standard.Conc}}{}
#' \item{\code{Test.Target.Conc}}{}
#' \item{\code{ISTD.Name}}{}
#' \item{\code{ISTD.Conc}}{}
#' \item{\code{IS.Area}}{}
#' \item{\code{Series}}{}
#' \item{\code{Area}}{}
#' \item{\code{Response}}{}
#' \item{\code{Analysis.Method}}{}
#' \item{\code{Analysis.Instrument}}{}
#' \item{\code{Analysis.Parameters}}{}
#'   }
#' @references
#' \insertRef{howard2010plasma}{invitroTKstats}
#' @keywords data
"kreutz2020"   

#' Cyprotex 2021 Caco2 Data Set 
#' 
#' Mass Spectrometry measurments of Caco2 membrane permeability as measured by
#' Cyprotex, a contractor to the U.S. EPA.
#' 
#' @name TO1caco2
#' @aliases to1caco2
#' @docType data
#' @format A data.table with 3110 rows 18 variables: \describe{
#' \item{\code{SampleName}}{}
#' \item{\code{CompoundName}}{}
#' \item{\code{Feature}}{}          
#' \item{\code{Area}}{}
#' \item{\code{ISTD.Area}}{}
#' \item{\code{ISTDResponseRatio}}{}
#' \item{\code{TO}}{}
#' \item{\code{FileName}}{}
#' \item{\code{SheetName}}{}        
#' \item{\code{Type}}{}
#' \item{\code{Dilution.Factor}}{}
#' \item{\code{Direction}}{}       
#' \item{\code{Vol.Receiver}}{}
#' \item{\code{Vol.Donor}}{}
#' \item{\code{Date}}{}            
#' \item{\code{ISTD.Name}}{}
#' \item{\code{ISTD.Conc}}{}
#' \item{\code{Test.Target.Conc}}{} 
#'   }
#' @references
#' \insertRef{hubatsch2007determination}{invitroTKstats}
#' @keywords data
"TO1caco2"   
