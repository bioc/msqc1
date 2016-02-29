\name{msqc1_8rep}
\alias{msqc1_8rep}
\docType{data}
\title{
MSQC1 8 technical replicates
}
\value{
	          A \code{data.frame} and \code{msqc1} S3 object.
}

\description{
 The data set contains a 8 technical replicates meassured on five mass spec devices.
}
\usage{data("msqc1_8rep")}
\format{
  A data frame with 6272 observations on the following 16 variables.
  \describe{
    \item{\code{Replicate.Name}}{a factor with levels \code{04_MSQC1_in_yd_incl_iRT_1in40_rep1} \code{05_MSQC1_in_yd_incl_iRT_1in40_rep2} \code{06_MSQC1_in_yd_incl_iRT_1in40_rep3} \code{07_MSQC1_in_yd_incl_iRT_1in40_rep4} \code{08_MSQC1_in_yd_incl_iRT_1in40_rep5} \code{09_MSQC1_in_yd_incl_iRT_1in40_rep6} \code{10_MSQC1_in_yd_incl_iRT_1in40_rep7} \code{11_MSQC1_in_yd_incl_iRT_1in40_rep8} \code{20131204_001_HS_yeast} \code{20131204_002_HS_yeast} \code{20131204_003_HS_yeast} \code{20131204_003_HS_yeast_01} \code{20131204_004_HS_yeast} \code{20131204_004_HS_yeast_02} \code{20131204_005_HS_yeast} \code{20131204_005_HS_yeast_03} \code{20131204_006_HS_yeast} \code{20131204_006_HS_yeast_04} \code{20131204_007_HS_yeast} \code{20131204_008_HS_yeast} \code{20131204_008_HS_yeast_05} \code{20131204_009_HS_yeast_06} \code{20131204_010_HS_yeast_07} \code{20131204_011_HS_yeast_08} \code{20140212_006_MSQC1_QTrap_01} \code{20140212_007_MSQC1_QTrap_02} \code{20140212_008_MSQC1_QTrap_03} \code{20140212_009_MSQC1_QTrap_04} \code{20140212_011_MSQC1_QTrap_05} \code{20140212_012_MSQC1_QTrap_06} \code{20140212_013_MSQC1_QTrap_07} \code{20140212_014_MSQC1_QTrap_08} \code{20140326_01_MSQC1_Yeast_iRT_12MS2_1e5_70k} \code{20140326_02_MSQC1_Yeast_iRT_12MS2_1e5_70k} \code{20140326_03_MSQC1_Yeast_iRT_12MS2_1e5_70k} \code{20140326_04_MSQC1_Yeast_iRT_12MS2_1e5_70k} \code{20140326_05_MSQC1_Yeast_iRT_12MS2_1e5_70k} \code{20140326_06_MSQC1_Yeast_iRT_12MS2_1e5_70k} \code{20140326_07_MSQC1_Yeast_iRT_12MS2_1e5_70k} \code{20140326_08_MSQC1_Yeast_iRT_12MS2_1e5_70k}}
    \item{\code{File.Name}}{a factor with levels \code{20131204_001_HS_yeast.wiff} \code{20131204_002_HS_yeast.wiff} \code{20131204_003_HS_yeast_01.raw} \code{20131204_003_HS_yeast.wiff} \code{20131204_004_HS_yeast_02.raw} \code{20131204_004_HS_yeast.wiff} \code{20131204_005_HS_yeast_03.raw} \code{20131204_005_HS_yeast.wiff} \code{20131204_006_HS_yeast_04.raw} \code{20131204_006_HS_yeast.wiff} \code{20131204_007_HS_yeast.wiff} \code{20131204_008_HS_yeast_05.raw} \code{20131204_008_HS_yeast.wiff} \code{20131204_009_HS_yeast_06.raw} \code{20131204_010_HS_yeast_07.raw} \code{20131204_011_HS_yeast_08.raw} \code{20140212_006_MSQC1_QTrap_01.wiff} \code{20140212_007_MSQC1_QTrap_02.wiff} \code{20140212_008_MSQC1_QTrap_03.wiff} \code{20140212_009_MSQC1_QTrap_04.wiff} \code{20140212_011_MSQC1_QTrap_05.wiff} \code{20140212_012_MSQC1_QTrap_06.wiff} \code{20140212_013_MSQC1_QTrap_07.wiff} \code{20140212_014_MSQC1_QTrap_08.wiff} \code{20140326_01_MSQC1_Yeast_iRT_12MS2_1e5_70k.raw} \code{20140326_02_MSQC1_Yeast_iRT_12MS2_1e5_70k.raw} \code{20140326_03_MSQC1_Yeast_iRT_12MS2_1e5_70k.raw} \code{20140326_04_MSQC1_Yeast_iRT_12MS2_1e5_70k.raw} \code{20140326_05_MSQC1_Yeast_iRT_12MS2_1e5_70k.raw} \code{20140326_06_MSQC1_Yeast_iRT_12MS2_1e5_70k.raw} \code{20140326_07_MSQC1_Yeast_iRT_12MS2_1e5_70k.raw} \code{20140326_08_MSQC1_Yeast_iRT_12MS2_1e5_70k.raw} \code{20150526_04_MSQC1_in_yd_incl_iRT_1in40_rep1.raw} \code{20150526_05_MSQC1_in_yd_incl_iRT_1in40_rep2.raw} \code{20150526_06_MSQC1_in_yd_incl_iRT_1in40_rep3.raw} \code{20150526_07_MSQC1_in_yd_incl_iRT_1in40_rep4.raw} \code{20150526_08_MSQC1_in_yd_incl_iRT_1in40_rep5.raw} \code{20150526_09_MSQC1_in_yd_incl_iRT_1in40_rep6.raw} \code{20150526_10_MSQC1_in_yd_incl_iRT_1in40_rep7.raw} \code{20150526_11_MSQC1_in_yd_incl_iRT_1in40_rep8.raw}}
    \item{\code{Protein.Name}}{a factor with levels \code{CAH1/P00915} \code{CAH2/P00918} \code{Catalase/P04040} \code{CRP/P02741} \code{DT-Diaphorase(NQO1)/P15559} \code{iRT-C18 Standard Peptides} \code{PPIA/Q13427}}
    \item{\code{Peptide.Sequence}}{a factor with levels \code{ADVTPADFSEWSK} \code{ALIVLAHSER} \code{AVQQPDGLAVLGIFLK} \code{DGLDAASYYAPVR} \code{EGHLSPDIVAEQK} \code{ESDTSYVSLK} \code{FEDENFILK} \code{FSTVAGESGSADTVR} \code{GAGAFGYFEVTHDITK} \code{GAGSSEPVTGLDAK} \code{GGPFSDSYR} \code{GTFIIDPAAVIR} \code{GTFIIDPGGVIR} \code{GYSIFSYATK} \code{LFLQFGAQGSPFLK} \code{LGGNEQVTR} \code{NLSVEDAAR} \code{SADFTNFDPR} \code{TAENFR} \code{TPVISGGPYEYR} \code{TPVITGAPYEYR} \code{VEATFGVDESNAK} \code{VLDALQAIK} \code{VSFELFADK} \code{YILAGVENSK}}
    \item{\code{Isotope.Label.Type}}{a factor with levels \code{heavy} \code{light}}
    \item{\code{Precursor.Mz}}{a numeric vector}
    \item{\code{Precursor.Charge}}{a numeric vector}
    \item{\code{Product.Mz}}{a numeric vector}
    \item{\code{Product.Charge}}{a numeric vector}
    \item{\code{Fragment.Ion}}{a factor with levels \code{b8} \code{precursor} \code{precursor [M+1]} \code{precursor [M+2]} \code{y10} \code{y11} \code{y12} \code{y4} \code{y5} \code{y6} \code{y7} \code{y8} \code{y9}}
    \item{\code{Retention.Time}}{a numeric vector}
    \item{\code{Area}}{a numeric vector}
    \item{\code{Background}}{a numeric vector}
    \item{\code{relative.amount}}{a numeric vector}
    \item{\code{instrument}}{a factor with levels \code{QExactive} \code{QExactiveHF} \code{QTRAP} \code{TRIPLETOF} \code{TSQVantage}}
    \item{\code{File.Name.Id}}{a numeric vector}
  }
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(msqc1_8rep)
## maybe str(msqc1_8rep) ; plot(msqc1_8rep) ...
}
\keyword{datasets}
