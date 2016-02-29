\name{msqc1-package}
\alias{msqc1-package}
\alias{msqc1}
\docType{package}
\title{
\packageTitle{msqc1}
}
\description{
\packageDescription{msqc1}


}
\details{
To  obtain  the  MSQC1  stock  solution,  MS  Qual/Quant  QC  Mix  (Catalog  #:  MSQC1,  Lot  #  081M6281)  was  purchased from SigmaLAldrich (Buchs, Switzerland). 

The first sample \code{msqc1_8rep} is based on eight techniqal replicates.  
To  derive  a  complex  sample  matrix,  10  mul  of  MSQC1  working  solution  was  mixed  with  10mul  (approx 10  mug)  tryptic  yeast  digest  (Saccharomyces  cerevisiae),  1mul  iRT  peptide  stock  solution  (Biognosys,  Schlieren,  Switzerland),  and  19  mul  0.1\%  FA. This  sample  is  denoted  as  standard  sample.

The \code{msqc1_dil} sample was derived from a six  point  dilution  series  containing  relative  MSQC1  amounts  of  0.025,  0.05,  0.2,  1,  2,  and  5  with  respect  to  the  standard  sample.  The  amount  of  yeast  digest  and  iRT peptides was kept constant across all dilution steps. The reference L:H ratio vs. the on column amount of  SIL peptide can be found in the data set \code{peptides}.

    Both samples were meassured on five mass spectrometers.
    The mass spectrometers and operation modes are:

    \itemize{

    \item \code{QTRAP} (ABSciex,  Concord,  Canada)  - {SRM mode}

    \item \code{TSQvantage} (ThermoScientific,  Bremen,  Germany) - {SRM mode}

    \item \code{Qexactive} (ThermoScientific,  Bremen,  Germany) - {PRM mode}
    
    \item \code{QExactiveHF} (ThermoScientific,  Bremen,  Germany) - {DIA mode}

    \item \code{TRIPLETOF} 5600 (ABSciex,  Concord,  Canada) - {SWATH and MS1 mode}
    }
    
All Raw  LC-MS  data  from  all  platforms  were  imported  into  Skyline  3.1 (\doi{10.1093/bioinformatics/btq054}).
A  single  expert  user  adjusted  the  automatic  peak  group  selection  and  integration  border  positioning  (skyline  legacy  peak  scoring).

The data were exported from skyline as csv files and included in this package as \code{data.frame} in RData files. 
}

\author{
\packageAuthor{msqc1}

Maintainer: \packageMaintainer{msqc1}

}

\references{
 \itemize{
 
 \item \url{www.sigmaLaldrich.com}
 
 \item skyline software - \doi{10.1093/bioinformatics/btq054}
 
  \item peptide mix \url{http://www.sigmaaldrich.com/life-science/proteomics/mass-spectrometry/ms-qual-quant-qc-mix.html}, 2014
  
  \item Functional Genomic Center Zurich data repository (internal project p1959) \url{http://fgcz-bfabric.uzh.ch/bfabric/userlab/show-project.html?projectId=1959} - contains all mass spectrometer generated raw files. access for registered users.
  
  \item \url{https://panoramaweb.org/labkey/project/Panorama Public/} search for msqc1 (not published yet) - contains all sklyline files and the export template
  }
  
}

\keyword{msqc1, proteomics, mass spectrometry, quantification, PRM, SWATH-MS, DIA}

%\seealso{

%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
 # Have Fun!
browseVignettes('msqc1')
}
