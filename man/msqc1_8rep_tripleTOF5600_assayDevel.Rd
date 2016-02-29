\name{msqc1_8rep_tripleTOF5600_assayDevel}
\alias{msqc1_8rep_tripleTOF5600_assayDevel}
\docType{data}
\title{
tripleTOF5600 peptide signal response
}
\value{
	          A \code{data.frame} and \code{msqc1} S3 object.
}

\description{
This data set contains peptide signal response for the tripleTOF5600 mass spec device.
The peptide level signals were measured  in SWATH-MS mode.

%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("msqc1_8rep_tripleTOF5600_assayDevel")}
\format{
  A data frame with 1616 observations on the following 17 variables.
  \describe{
    \item{\code{Replicate.Name}}{a factor with levels \code{20131204_001_HS_yeast} \code{20131204_002_HS_yeast} \code{20131204_003_HS_yeast} \code{20131204_004_HS_yeast} \code{20131204_005_HS_yeast} \code{20131204_006_HS_yeast} \code{20131204_007_HS_yeast} \code{20131204_008_HS_yeast}}
    \item{\code{File.Name}}{a factor with levels \code{20131204_001_HS_yeast.wiff} \code{20131204_002_HS_yeast.wiff} \code{20131204_003_HS_yeast.wiff} \code{20131204_004_HS_yeast.wiff} \code{20131204_005_HS_yeast.wiff} \code{20131204_006_HS_yeast.wiff} \code{20131204_007_HS_yeast.wiff} \code{20131204_008_HS_yeast.wiff}}
    \item{\code{Protein.Name}}{a factor with levels \code{CAH1/P00915} \code{CAH2/P00918} \code{Catalase/P04040} \code{CRP/P02741} \code{DT-Diaphorase(NQO1)/P15559} \code{iRT-C18 Standard Peptides} \code{P00915MSQC1|CAH1_HUMAN_MSQC1} \code{P00918MSQC1|CAH2_HUMAN_MSQC1} \code{P02741MSQC1|CRP_HUMAN_MSQC1} \code{P04040MSQC1|CATA_HUMAN_MSQC1} \code{P15559MSQC1|NQO1_HUMAN_MSQC1} \code{P62937MSQC1|PPIA_HUMAN_MSQC1} \code{PPIA/Q13427}}
    \item{\code{Peptide.Sequence}}{a factor with levels \code{ADGLAVIGVLMK} \code{ADVLTTGAGNPVGDK} \code{ADVTPADFSEWSK} \code{AFYVNVLNEEQR} \code{ALIVLAHSER} \code{APLTKPLK} \code{AVQQPDGLAVLGIFLK} \code{DAQIFIQK} \code{DGLDAASYYAPVR} \code{DLFNAIATGK} \code{DPILFPSFIHSQK} \code{DYPLIPVGK} \code{EGHLSPDIVAEQK} \code{EGMNIVEAMER} \code{EPISVSSEQVLK} \code{ESDTSYVSLK} \code{ESISVSSEQLAQFR} \code{FEDENFILK} \code{FGLSVGHHLGK} \code{FNTANDDNVTQVR} \code{FSTVAGESGSADTVR} \code{GAGAFGYFEVTHDITK} \code{GAGSSEPVTGLDAK} \code{GGPFSDSYR} \code{GGPLDGTYR} \code{GTFIIDPAAVIR} \code{GTFIIDPGGVIR} \code{GYSIFSYATK} \code{HDTSLKPISVSYNPATAK} \code{ILNNGHAFNVEFDDSQDK} \code{IQILEGWK} \code{LFAYPDTHR} \code{LFLQFGAQGSPFLK} \code{LGGNEQVTR} \code{LNVITVGPR} \code{LSQEDPDYGIR} \code{LYPIANGNNQSPVDIK} \code{NAIHTFVQSGSHLAAR} \code{NGPEQWSK} \code{NLSVEDAAR} \code{QSPVDIDTHTAK} \code{SADFTNFDPR} \code{TEWLDGK} \code{TPVISGGPYEYR} \code{TPVITGAPYEYR} \code{TSFNYAMK} \code{VEATFGVDESNAK} \code{VGSAKPGLQK} \code{VLDALQAIK} \code{VSFELFADK} \code{VVDVLDSIK} \code{YAAELHLVHWNTK} \code{YDPSLKPLSVSYDQATSLR} \code{YILAGVENSK} \code{YSAELHVAHWNSAK} \code{YSSLAEAASK}}
    \item{\code{Isotope.Label.Type}}{a factor with levels \code{light}}
    \item{\code{Precursor.Mz}}{a numeric vector}
    \item{\code{Precursor.Charge}}{a numeric vector}
    \item{\code{Product.Mz}}{a numeric vector}
    \item{\code{Product.Charge}}{a numeric vector}
    \item{\code{Fragment.Ion}}{a factor with levels \code{b2} \code{b3} \code{b7} \code{b8} \code{precursor} \code{precursor [M+1]} \code{precursor [M+2]} \code{y10} \code{y11} \code{y12} \code{y13} \code{y14} \code{y2} \code{y3} \code{y4} \code{y5} \code{y6} \code{y7} \code{y8} \code{y9}}
    \item{\code{Retention.Time}}{a numeric vector}
    \item{\code{Area}}{a numeric vector}
    \item{\code{Background}}{a numeric vector}
    \item{\code{relative.amount}}{a numeric vector}
    \item{\code{flag}}{a factor with levels \code{False} \code{True}}
    \item{\code{instrument}}{a factor with levels \code{TRIPLETOF}}
    \item{\code{File.Name.Id}}{a numeric vector}
  }
}
\details{
The signal response of the peptide *TAENF[R]* was not measured due to the lower m/z limit used by the DIA/SWATH-MS methods.
}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(msqc1_8rep_tripleTOF5600_assayDevel)

#filter <- !(msqc1_8rep_tripleTOF5600_assayDevel$Protein.Name %in% peptides$Protein.Name)

#msqc1:::.signal_response_ident_plot(msqc1_8rep_tripleTOF5600_assayDevel[filter, ], 
#  instrument="tripleTOF560",  
#  prot="P02741")

}
\keyword{datasets}
