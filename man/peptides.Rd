\name{peptides}
\alias{peptides}
\docType{data}
\title{
sigma mix peptides
}
\value{
	          A \code{data.frame}.
}

\description{
This table contains msqc1 peptide proberties.
}
\usage{data("peptides")}
\format{
  A data frame with 14 observations on the following 5 variables.
  \describe{
    \item{\code{Peptide.Sequence}}{a factor with levels \code{ALIVLAHSER} \code{AVQQPDGLAVLGIFLK} \code{EGHLSPDIVAEQK} \code{ESDTSYVSLK} \code{FEDENFILK} \code{FSTVAGESGSADTVR} \code{GAGAFGYFEVTHDITK} \code{GGPFSDSYR} \code{GYSIFSYATK} \code{NLSVEDAAR} \code{SADFTNFDPR} \code{TAENFR} \code{VLDALQAIK} \code{VSFELFADK}}
    \item{\code{SIL.peptide.per.vial}}{a numeric vector}
    \item{\code{LH.ratio}}{a numeric vector}
    \item{\code{actual.LH.ratio}}{a numeric vector}
    \item{\code{Protein.Name}}{a factor with levels \code{CAH1/P00915} \code{CAH2/P00918} \code{Catalase/P04040} \code{CRP/P02741} \code{DT-Diaphorase(NQO1)/P15559} \code{PPIA/Q13427}}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
\url{http://www.sigmaaldrich.com/life-science/proteomics/mass-spectrometry/ms-qual-quant-qc-mix.html}, 2014
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(peptides)
## maybe str(peptides) ; plot(peptides) ...
}
\keyword{datasets}
