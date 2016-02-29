\name{msqc1_8rep_QEXACTIVEHF_assayDevel}
\alias{msqc1_8rep_QEXACTIVEHF_assayDevel}
\docType{data}
\title{
QEXACTIVEHF peptide signal response. 
}
\value{
	          A \code{data.frame} and \code{msqc1} S3 object.
}

\description{
This data set contains peptide signal response for the QEXACTIVEHF mass spec device.
The peptide level signals as measured by DIA. 
}
\usage{data("msqc1_8rep_QEXACTIVEHF_assayDevel")}
\format{
  A data frame with 1864 observations on the following 16 variables.
  \describe{
    \item{\code{Replicate.Name}}{a factor with levels \code{04_MSQC1_in_yd_incl_iRT_1in40_rep1} \code{05_MSQC1_in_yd_incl_iRT_1in40_rep2} \code{06_MSQC1_in_yd_incl_iRT_1in40_rep3} \code{07_MSQC1_in_yd_incl_iRT_1in40_rep4} \code{08_MSQC1_in_yd_incl_iRT_1in40_rep5} \code{09_MSQC1_in_yd_incl_iRT_1in40_rep6} \code{10_MSQC1_in_yd_incl_iRT_1in40_rep7} \code{11_MSQC1_in_yd_incl_iRT_1in40_rep8}}
    \item{\code{File.Name}}{a factor with levels \code{20150526_04_MSQC1_in_yd_incl_iRT_1in40_rep1.raw} \code{20150526_05_MSQC1_in_yd_incl_iRT_1in40_rep2.raw} \code{20150526_06_MSQC1_in_yd_incl_iRT_1in40_rep3.raw} \code{20150526_07_MSQC1_in_yd_incl_iRT_1in40_rep4.raw} \code{20150526_08_MSQC1_in_yd_incl_iRT_1in40_rep5.raw} \code{20150526_09_MSQC1_in_yd_incl_iRT_1in40_rep6.raw} \code{20150526_10_MSQC1_in_yd_incl_iRT_1in40_rep7.raw} \code{20150526_11_MSQC1_in_yd_incl_iRT_1in40_rep8.raw}}
    \item{\code{Protein.Name}}{a factor with levels \code{CAH1/P00915} \code{CAH2/P00918} \code{Catalase/P04040} \code{CRP/P02741} \code{DT-Diaphorase(NQO1)/P15559} \code{iRT-C18 Standard Peptides} \code{P00915MSQC1|CAH1_HUMAN_MSQC1} \code{P00918MSQC1|CAH2_HUMAN_MSQC1} \code{P02741MSQC1|CRP_HUMAN_MSQC1} \code{P04040MSQC1|CATA_HUMAN_MSQC1} \code{P15559MSQC1|NQO1_HUMAN_MSQC1} \code{P62937MSQC1|PPIA_HUMAN_MSQC1} \code{PPIA/Q13427}}
    \item{\code{Peptide.Sequence}}{a factor with levels \code{ADVLTTGAGNPVGDK} \code{ADVTPADFSEWSK} \code{AFYVNVLNEEQR} \code{ALIVLAHSER} \code{APLTKPLK} \code{AVQQPDGLAVLGIFLK} \code{DAQIFIQK} \code{DGLDAASYYAPVR} \code{DLFNAIATGK} \code{DPANFQYPAESVLAYK} \code{DPILFPSFIHSQK} \code{DYPLIPVGK} \code{EGHLSPDIVAEQK} \code{EIINVGHSFHVNFEDNDNR} \code{EPISVSSEQVLK} \code{ESDTSYVSLK} \code{ESISVSSEQLAQFR} \code{FEDENFILK} \code{FGLSVGHHLGK} \code{FNTANDDNVTQVR} \code{FSTVAGESGSADTVR} \code{GAGAFGYFEVTHDITK} \code{GAGSSEPVTGLDAK} \code{GGPFSDSYR} \code{GGPLDGTYR} \code{GTFIIDPAAVIR} \code{GTFIIDPGGVIR} \code{GYSIFSYATK} \code{HDTSLKPISVSYNPATAK} \code{HNGPEHWHK} \code{ILNNGHAFNVEFDDSQDK} \code{IQILEGWK} \code{LCENIAGHLK} \code{LFAYPDTHR} \code{LFLQFGAQGSPFLK} \code{LFQFHFHWGSTNEHGSEHTVDGVK} \code{LGGNEQVTR} \code{LGPNYLHIPVNCPYR} \code{LIQFHFHWGSLDGQGSEHTVDK} \code{LNVITVGPR} \code{LSQEDPDYGIR} \code{LVNANGEAVYCK} \code{LYPIANGNNQSPVDIK} \code{NAIHTFVQSGSHLAAR} \code{NGPEQWSK} \code{NLSVEDAAR} \code{QDNEILIFWSK} \code{QSPVDIDTHTAK} \code{SADFTNFDPR} \code{SIPTDNQIK} \code{TEWLDGK} \code{TPVISGGPYEYR} \code{TPVITGAPYEYR} \code{VEATFGVDESNAK} \code{VFEHIGK} \code{VGSAKPGLQK} \code{VLDALQAIK} \code{VSFELFADK} \code{VVDVLDSIK} \code{YAAELHLVHWNTK} \code{YDPSLKPLSVSYDQATSLR} \code{YEVQGEVFTKPQLWP} \code{YILAGVENSK} \code{YNAEKPK} \code{YSAELHVAHWNSAK} \code{YSSLAEAASK}}
    \item{\code{Isotope.Label.Type}}{a factor with levels \code{light}}
    \item{\code{Precursor.Mz}}{a numeric vector}
    \item{\code{Precursor.Charge}}{a numeric vector}
    \item{\code{Product.Mz}}{a numeric vector}
    \item{\code{Product.Charge}}{a numeric vector}
    \item{\code{Fragment.Ion}}{a factor with levels \code{b3} \code{b4} \code{b5} \code{b6} \code{b7} \code{b8} \code{y10} \code{y11} \code{y12} \code{y3} \code{y4} \code{y5} \code{y6} \code{y7} \code{y8} \code{y9}}
    \item{\code{Retention.Time}}{a numeric vector}
    \item{\code{Area}}{a numeric vector}
    \item{\code{Background}}{a numeric vector}
    \item{\code{relative.amount}}{a numeric vector}
    \item{\code{QEXACTIVEHF}}{a factor with levels \code{QEXACTIVEHF}}
    \item{\code{File.Name.Id}}{a numeric vector}
  }
}
%\details{
%%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(msqc1_8rep_QEXACTIVEHF_assayDevel)
## maybe str(msqc1_8rep_QEXACTIVEHF_assayDevel) ; plot(msqc1_8rep_QEXACTIVEHF_assayDevel) ...
}
\keyword{datasets}
