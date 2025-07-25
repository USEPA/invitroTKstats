% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_caco2.R
\name{format_caco2}
\alias{format_caco2}
\title{Creates a Standardized Data Frame with Caco-2 Data (Level-1)}
\usage{
format_caco2(
  FILENAME = "MYDATA",
  data.in,
  sample.col = "Lab.Sample.Name",
  lab.compound.col = "Lab.Compound.Name",
  dtxsid.col = "DTXSID",
  date = NULL,
  date.col = "Date",
  compound.col = "Compound.Name",
  area.col = "Area",
  istd.col = "ISTD.Area",
  type.col = "Type",
  direction.col = "Direction",
  membrane.area = NULL,
  membrane.area.col = "Membrane.Area",
  receiver.vol.col = "Vol.Receiver",
  donor.vol.col = "Vol.Donor",
  test.conc = NULL,
  test.conc.col = "Test.Compound.Conc",
  cal = NULL,
  cal.col = "Cal",
  dilution = NULL,
  dilution.col = "Dilution.Factor",
  time = NULL,
  time.col = "Time",
  istd.name = NULL,
  istd.name.col = "ISTD.Name",
  istd.conc = NULL,
  istd.conc.col = "ISTD.Conc",
  test.nominal.conc = NULL,
  test.nominal.conc.col = "Test.Target.Conc",
  biological.replicates = NULL,
  biological.replicates.col = "Biological.Replicates",
  technical.replicates = NULL,
  technical.replicates.col = "Technical.Replicates",
  analysis.method = NULL,
  analysis.method.col = "Analysis.Method",
  analysis.instrument = NULL,
  analysis.instrument.col = "Analysis.Instrument",
  analysis.parameters = NULL,
  analysis.parameters.col = "Analysis.Parameters",
  note.col = "Note",
  level0.file = NULL,
  level0.file.col = "Level0.File",
  level0.sheet = NULL,
  level0.sheet.col = "Level0.Sheet",
  output.res = FALSE,
  save.bad.types = FALSE,
  sig.figs = 5,
  INPUT.DIR = NULL,
  OUTPUT.DIR = NULL
)
}
\arguments{
\item{FILENAME}{(Character) A string used to identify the output level-1 file.
"<FILENAME>-Caco-2-Level1.tsv", and/or used to identify the input level-0 file,
"<FILENAME>-Caco-2-Level0.tsv" if importing from a .tsv file. (Defaults to "MYDATA".)}

\item{data.in}{(Data Frame) A level-0 data frame containing
mass-spectrometry peak areas, indication of chemical identity,
and measurement type. The data frame should
contain columns with names specified by the following arguments:}

\item{sample.col}{(Character) Column name of \code{data.in} containing the unique mass
spectrometry (MS) sample name used by the laboratory. (Defaults to
"Lab.Sample.Name".)}

\item{lab.compound.col}{(Character) Column name of \code{data.in} containing the test compound
name used by the laboratory. (Defaults to "Lab.Compound.Name".)}

\item{dtxsid.col}{(Character) Column name of \code{data.in} containing EPA's DSSTox Structure
ID (\url{http://comptox.epa.gov/dashboard}). (Defaults to "DTXSID".)}

\item{date}{(Character) The laboratory measurement date, format "MMDDYY" where 
"MM" = 2 digit month, "DD" = 2 digit day, and "YY" = 2 digit year. (Defaults to \code{NULL}.) 
(Note: Single entry only, use only if all data were collected on the same date.)}

\item{date.col}{(Character) Column name containing \code{date} information. (Defaults to "Date".) (Note: \code{data.in} does not
necessarily have this field. If this field is missing, it can be auto-filled with the value 
specified in \code{date}.)}

\item{compound.col}{(Character) Column name of \code{data.in} containing the test compound.
(Defaults to "Compound.Name".)}

\item{area.col}{(Character) Column name of \code{data.in} containing the target analyte (that
is, the test compound) MS peak area. (Defaults to "Area".)}

\item{istd.col}{(Character) Column name of \code{data.in} containing the
MS peak area for the internal standard. (Defaults to "ISTD.Area".)}

\item{type.col}{(Character) Column name of \code{data.in} containing the sample type (see table
under Details). (Defaults to "Type".)}

\item{direction.col}{(Character) Column name of \code{data.in} containing the direction of
the Caco-2 permeability experiment: either apical donor to basolateral receiver (AtoB), or 
basolateral donor to apical receiver (BtoA). (Defaults to "Direction".)}

\item{membrane.area}{(Numeric) The area of the Caco-2 monolayer (in cm^2). 
(Defaults to \code{NULL}.) (Note: Single entry only, use only if all tested compounds 
have the same area for the Caco-2 monolayer.)}

\item{membrane.area.col}{(Character) Column name containing \code{membrane.area} 
information. (Defaults to "Membrane.Area".) (Note: \code{data.in} does not
necessarily have this field. If this field is missing, it can be auto-filled with the value 
specified in \code{membrane.area}.)}

\item{receiver.vol.col}{(Character) Column name of \code{data.in} containing the media volume
(in cm^3) of the receiver portion of the Caco-2 experimental well. 
(Defaults to "Vol.Receiver".)}

\item{donor.vol.col}{(Character) Column name of \code{data.in} containing the media volume
(in cm^3) of the donor portion of the Caco-2 experimental well where the
test chemical is added. (Defaults to "Vol.Donor".)}

\item{test.conc}{(Numeric) The standard test chemical concentration for the
Caco-2 assay. (Defaults to \code{NULL}.) (Note: Single entry only, 
use only if the same standard concentration was used for all tested compounds.)}

\item{test.conc.col}{(Character) Column name containing \code{test.conc} 
information. (Defaults to "Test.Compound.Conc".) (Note: \code{data.in} does not
necessarily have this field. If this field is missing, it can be auto-filled with the value 
specified in \code{test.conc}.)}

\item{cal}{(Character) MS calibration the samples were based on. Typically, this uses 
indices or dates to represent if the analyses were done on different machines on 
the same day or on different days with the same MS analyzer. (Defaults to \code{NULL}.) 
(Note: Single entry only, 
use only if all data were collected based on the same calibration.)}

\item{cal.col}{(Character) Column name containing \code{cal} 
information. (Defaults to "Cal".) (Note: \code{data.in} does not
necessarily have this field. If this field is missing, it can be auto-filled with the value 
specified in \code{cal}.)}

\item{dilution}{(Numeric) Number of times the sample was diluted before MS 
analysis. (Defaults to \code{NULL}.) (Note: Single entry only, use only if all 
samples underwent the same number of dilutions.)}

\item{dilution.col}{(Character) Column name containing \code{dilution} 
information. (Defaults to "Dilution.Factor".) (Note: \code{data.in} does not
necessarily have this field. If this field is missing, it can be auto-filled with the value 
specified in \code{dilution}.)}

\item{time}{(Numeric) The amount of time (in hours) before the receiver and donor 
compartments are measured. (Defaults to \code{NULL}.)}

\item{time.col}{(Character) Column name containing \code{meas.time} 
information. (Defaults to "Time".) (Note: \code{data.in} does not
necessarily have this field. If this field is missing, it can be auto-filled with the value 
specified in \code{meas.time}.)}

\item{istd.name}{(Character) The identity of the internal standard. (Defaults to \code{NULL}.) 
(Note: Single entry only, use only if all tested compounds use the same internal standard.)}

\item{istd.name.col}{(Character) Column name containing \code{istd.name} information. (Defaults to "ISTD.Name".) 
(Note: \code{data.in} does not necessarily have this field. If this field is missing, 
it can be auto-filled with the value specified in \code{istd.name}.)}

\item{istd.conc}{(Numeric) The concentration for the internal standard. (Defaults to \code{NULL}.) 
(Note: Single entry only, use only if all tested compounds have the same 
internal standard concentration.)}

\item{istd.conc.col}{(Character) Column name containing \code{istd.conc} information. (Defaults to "ISTD.Conc".) 
(Note: \code{data.in} does not necessarily have this field. If this field is missing, 
it can be auto-filled with the value specified in \code{istd.conc}.)}

\item{test.nominal.conc}{(Numeric) The nominal concentration added to the donor
compartment at time 0. (Defaults to \code{NULL}.) (Note: Single entry only,
use only if all tested compounds used the same concentration at time 0.}

\item{test.nominal.conc.col}{(Character) Column name containing \code{test.nominal.conc} 
information. (Defaults to "Test.Target.Conc".) (Note: \code{data.in} does not
necessarily have this field. If this field is missing, it can be auto-filled with the value 
specified in \code{test.nominal.conc}.)}

\item{biological.replicates}{(Character) Replicates with the same analyte. Typically, this uses 
numbers or letters to index. (Defaults to \code{NULL}.) (Note: Single entry only, 
use only if none of the test compounds have replicates.)}

\item{biological.replicates.col}{(Character) Column name of \code{data.in} containing the number or 
the indices of replicates with the same analyte. (Defaults to "Biological.Replicates".)
(Note: \code{data.in} does not necessarily have this field. If this field is missing, it can be auto-filled with the value 
specified in \code{biological.replicates}.)}

\item{technical.replicates}{(Character) Repeated measurements from one sample. Typically, this uses 
numbers or letters to index. (Defaults to \code{NULL}.) (Note: Single entry only, 
use only if none of the test compounds have replicates.)}

\item{technical.replicates.col}{(Character) Column name of \code{data.in} containing the number or 
the indices of replicates taken from the one sample. (Defaults to "Technical.Replicates".) 
(Note: \code{data.in} does not necessarily have this field. If this field is missing, it can be auto-filled with the value 
specified in \code{technical.replicates}.)}

\item{analysis.method}{(Character) The analytical chemistry analysis method, 
typically "LCMS" or "GCMS", liquid chromatography or gas chromatography–mass
spectrometry, respectively. (Defaults to \code{NULL}.)
(Note: Single entry only, use only if the same method was used for all tested compounds.)}

\item{analysis.method.col}{(Character) Column name containing \code{analysis.method} 
information. (Defaults to "Analysis.Method".) (Note: \code{data.in} does not
necessarily have this field. If this field is missing, it can be auto-filled with the value 
specified in \code{analysis.method}.)}

\item{analysis.instrument}{(Character) The instrument used for chemical analysis, 
for example "Agilent 6890 GC with model 5973 MS". (Defaults to \code{NULL}.) 
(Note: Single entry only, use only if the same instrument was used for all tested compounds.)}

\item{analysis.instrument.col}{(Character) Column name containing \code{analysis.instrument} 
information. (Defaults to "Analysis.Instrument".) (Note: \code{data.in} does not
necessarily have this field. If this field is missing, it can be auto-filled with the value 
specified in \code{analysis.instrument}.)}

\item{analysis.parameters}{(Character) The parameters used to identify the 
compound on the chemical analysis instrument, for example
"Negative Mode, 221.6/161.6, -DPb=26, FPc=-200, EPd=-10, CEe=-20, CXPf=-25.0". (Defaults to \code{NULL}.) 
(Note: Single entry only, use only if the same parameters were used for all tested compounds.)}

\item{analysis.parameters.col}{(Character) Column name containing \code{analysis.parameters} 
information. (Defaults to "Analysis.Parameters".) (Note: \code{data.in} does not
necessarily have this field. If this field is missing, it can be auto-filled with the value 
specified in \code{analysis.parameters}.)}

\item{note.col}{(Character) Column name of \code{data.in} containing additional notes on 
test compounds. (Defaults to "Note").}

\item{level0.file}{(Character) The level-0 file from which the \code{data.in} were obtained.
(Defaults to \code{NULL}.) (Note: Single entry only, use only if all rows in data.in
were obtained from the same level-0 file.)}

\item{level0.file.col}{(Character) Column name containing \code{level0.file} information. 
(Defaults to "Level0.File".) (Note: \code{data.in} does not
necessarily have this field. If this field is missing, it can be auto-filled with the value 
specified in \code{level0.file}.)}

\item{level0.sheet}{(Character) The specific sheet name of level-0 file from which the 
\code{data.in} is obtained from, if the level-0 file is an Excel workbook. 
(Defaults to \code{NULL}.) (Note: Single entry only, use only if all rows in \code{data.in}
were obtained from the same sheet in the same level-0 file.)}

\item{level0.sheet.col}{(Character) Column name containing \code{level0.sheet} information.
(Defaults to "Level0.Sheet".) (Note: \code{data.in} does not
necessarily have this field. If this field is missing, it can be auto-filled with the value 
specified in \code{level0.sheet}.)}

\item{output.res}{(Logical) When set to \code{TRUE}, the result 
table (level-1) will be exported to the user's per-session temporary 
directory or \code{OUTPUT.DIR} (if specified) as a .tsv file. 
(Defaults to \code{FALSE}.)}

\item{save.bad.types}{(Logical) When set to \code{TRUE}, export data removed 
due to inappropriate sample types. See the Detail section for the required sample types. 
(Defaults to \code{FALSE}.)}

\item{sig.figs}{(Numeric) The number of significant figures to round the exported result table (level-1). 
(Defaults to \code{5}.)}

\item{INPUT.DIR}{(Character) Path to the directory where the input level-0 file exists. 
If \code{NULL}, looking for the input level-0 file in the current working
directory. (Defaults to \code{NULL}.)}

\item{OUTPUT.DIR}{(Character) Path to the directory to save the output file. 
If \code{NULL}, the output file will be saved to the user's per-session temporary
directory or \code{INPUT.DIR} if specified. (Defaults to \code{NULL}.)}
}
\value{
A level-1 data frame with a standardized format containing a  
standardized set of columns and column names with membrane permeability data
from a Caco-2 assay.
}
\description{
This function formats data describing mass spectrometry (MS) peak areas
from samples collected as part of \emph{in vitro} measurements of membrane
permeability using Caco-2 cells \insertCite{hubatsch2007determination}{invitroTKstats}.
The input data frame is organized into a standard set of columns and is
written to a tab-separated text file.
}
\details{
In this experiment an
\emph{in vitro} well is separated into two by a membrane composed of a monolayer of
Caco-2 cells. A test chemical is added to either the apical or basolateral side of
of the monolayer at time 0, and after a set time samples are taken from both
the "donor" (side where the test chemical was added) and the "receiver" side.
Depending on the direction of the test the donor side can be either apical or
basolateral.

The data frame of observations should be annotated according to direction
(either apical to basolateral -- "AtoB" -- or basolateral to apical -- "BtoA") and type
of concentration measured:
\tabular{rr}{
  Blank with no chemical added \tab Blank \cr
  Target concentration added to donor compartment at time 0 (C0) \tab D0\cr
  Donor compartment at end of experiment \tab D2\cr
  Receiver compartment at end of experiment\tab R2\cr
}

Chemical concentration is calculated qualitatively as a response and 
returned as a column in the output data frame:

Response <- AREA / ISTD.AREA * ISTD.CONC

If the output level-1 result table is chosen to be exported and an output 
directory is not specified, it will be exported to the user's R session
temporary directory. This temporary directory is a per-session directory 
whose path can be found with the following code: \code{tempdir()}. For more 
details, see \url{https://www.collinberke.com/til/posts/2023-10-24-temp-directories/}.

As a best practice, \code{INPUT.DIR} and/or \code{OUTPUT.DIR} should be 
specified to simplify the process of importing and exporting files. This 
practice ensures that the exported files can easily be found and will not be 
exported to a temporary directory.
}
\examples{
## Load example level-0 data and do not export the result table
level0 <- invitroTKstats::caco2_L0
level1 <- format_caco2(data.in = level0,
                       sample.col = "Sample",
                       lab.compound.col = "Lab.Compound.ID",
                       compound.col = "Compound",
                       area.col = "Peak.Area",
                       istd.col = "ISTD.Peak.Area",
                       membrane.area = 0.11,
                       test.conc.col = "Compound.Conc",
                       cal = 1, 
                       time = 2, 
                       istd.conc = 1, 
                       test.nominal.conc = 10, 
                       biological.replicates = 1, 
                       technical.replicates = 1,
                       analysis.method.col = "Analysis.Params",
                       analysis.instrument = "Agilent.GCMS",
                       analysis.parameters = "Unknown",
                       note.col = NULL,
                       output.res = FALSE
)

}
\references{
\insertRef{hubatsch2007determination}{invitroTKstats}
}
\author{
John Wambaugh
}
