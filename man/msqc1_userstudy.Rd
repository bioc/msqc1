\name{msqc1_userstudy}
\alias{msqc1_userstudy}
\docType{data}
\title{
msqc1 user study
}
\value{
	          A \code{data.frame} and \code{msqc1} S3 object.
}

\description{
contains a user study
}
\usage{data("msqc1_userstudy")}
\format{
  A data frame with 21204 observations on the following 17 variables.
  \describe{
    \item{\code{Replicate.Name}}{a factor with levels \code{01_QCdil_40dil_1} \code{02_QCdil_40dil_2} \code{03_MSQC1_1_40dil_1} \code{03_QCdil_40dil_3} \code{04_MSQC1_1_40dil_1} \code{04_QCdil_20dil_1} \code{05_MSQC1_1_40dil_3} \code{05_QCdil_20dil_2} \code{06_QCdil_20dil_3} \code{07_MSQC1_1_20dil_1} \code{07_QCdil_5dil_1} \code{08_MSQC1_1_20dil_2} \code{08_QCdil_5dil_2} \code{09_MSQC1_1_20dil_3} \code{09_QCdil_5dil_3} \code{10_nodil_1} \code{11_MSQC1_NoDil_1} \code{11_nodil_2} \code{12_nodil_3} \code{13_2conc_1_re} \code{13_MSQC1_NoDil_2} \code{14_2conc_2_re} \code{15_2conc_3} \code{15_MSQC1_NoDil_3} \code{16_5conc_1} \code{17_5conc_2} \code{17_MSQC1_x2_1} \code{18_5conc_3} \code{19_MSQC1_x2_2} \code{20140818_003_MSQC1_1_40dil_1} \code{20140818_004_MSQC1_1_40dil_1} \code{20140818_005_MSQC1_1_40dil_3} \code{20140818_007_MSQC1_1_20dil_1} \code{20140818_008_MSQC1_1_20dil_2} \code{20140818_009_MSQC1_1_20dil_3} \code{20140818_011_MSQC1_NoDil_1} \code{20140818_013_MSQC1_NoDil_2} \code{20140818_015_MSQC1_NoDil_3} \code{20140818_017_MSQC1_x2_1} \code{20140818_019_MSQC1_x2_2} \code{20140818_021_MSQC1_x2_3} \code{20140818_023_MSQC1_x5_1} \code{20140818_025_MSQC1_x5_2} \code{20140818_027_MSQC1_x5_3} \code{20140818_030_MSQC1_1_5dil_1} \code{20140818_032_MSQC1_1_5dil_2} \code{20140818_034_MSQC1_1_5dil_3} \code{20150601_001_QCdil_40dil_1} \code{20150601_002_QCdil_40dil_2} \code{20150601_003_QCdil_40dil_3} \code{20150601_004_QCdil_20dil_1} \code{20150601_005_QCdil_20dil_2} \code{20150601_006_QCdil_20dil_3} \code{20150601_007_QCdil_5dil_1} \code{20150601_008_QCdil_5dil_2} \code{20150601_009_QCdil_5dil_3} \code{20150601_010_nodil_1} \code{20150601_011_nodil_2} \code{20150601_012_nodil_3} \code{20150601_013_2conc_1_re} \code{20150601_014_2conc_2_re} \code{20150601_015_2conc_3} \code{20150601_016_5conc_1} \code{20150601_017_5conc_2} \code{20150601_018_5conc_3} \code{21_MSQC1_x2_3} \code{23_MSQC1_x5_1} \code{25_MSQC1_x5_2} \code{27_MSQC1_x5_3} \code{30_MSQC1_1_5dil_1} \code{32_MSQC1_1_5dil_2} \code{34_MSQC1_1_5dil_3}}
    \item{\code{File.Name}}{a factor with levels \code{20140818_003_MSQC1_1_40dil_1.wiff} \code{20140818_004_MSQC1_1_40dil_1.wiff} \code{20140818_005_MSQC1_1_40dil_3.wiff} \code{20140818_007_MSQC1_1_20dil_1.wiff} \code{20140818_008_MSQC1_1_20dil_2.wiff} \code{20140818_009_MSQC1_1_20dil_3.wiff} \code{20140818_011_MSQC1_NoDil_1.wiff} \code{20140818_013_MSQC1_NoDil_2.wiff} \code{20140818_015_MSQC1_NoDil_3.wiff} \code{20140818_017_MSQC1_x2_1.wiff} \code{20140818_019_MSQC1_x2_2.wiff} \code{20140818_021_MSQC1_x2_3.wiff} \code{20140818_023_MSQC1_x5_1.wiff} \code{20140818_025_MSQC1_x5_2.wiff} \code{20140818_027_MSQC1_x5_3.wiff} \code{20140818_030_MSQC1_1_5dil_1.wiff} \code{20140818_032_MSQC1_1_5dil_2.wiff} \code{20140818_034_MSQC1_1_5dil_3.wiff} \code{20150601_001_QCdil_40dil_1.wiff} \code{20150601_002_QCdil_40dil_2.wiff} \code{20150601_003_QCdil_40dil_3.wiff} \code{20150601_004_QCdil_20dil_1.wiff} \code{20150601_005_QCdil_20dil_2.wiff} \code{20150601_006_QCdil_20dil_3.wiff} \code{20150601_007_QCdil_5dil_1.wiff} \code{20150601_008_QCdil_5dil_2.wiff} \code{20150601_009_QCdil_5dil_3.wiff} \code{20150601_010_nodil_1.wiff} \code{20150601_011_nodil_2.wiff} \code{20150601_012_nodil_3.wiff} \code{20150601_013_2conc_1_re.wiff} \code{20150601_014_2conc_2_re.wiff} \code{20150601_015_2conc_3.wiff} \code{20150601_016_5conc_1.wiff} \code{20150601_017_5conc_2.wiff} \code{20150601_018_5conc_3.wiff}}
    \item{\code{Protein.Name}}{a factor with levels \code{CAH1/P00915} \code{CAH2/P00918} \code{Catalase/P04040} \code{CRP/P02741} \code{DT-Diaphorase(NQO1)/P15559} \code{iRT-C18 Standard Peptides} \code{PPIA/Q13427}}
    \item{\code{Peptide.Sequence}}{a factor with levels \code{ADVTPADFSEWSK} \code{ALIVLAHSER} \code{AVQQPDGLAVLGIFLK} \code{DGLDAASYYAPVR} \code{EGHLSPDIVAEQK} \code{ESDTSYVSLK} \code{FEDENFILK} \code{FSTVAGESGSADTVR} \code{GAGAFGYFEVTHDITK} \code{GAGSSEPVTGLDAK} \code{GGPFSDSYR} \code{GTFIIDPAAVIR} \code{GTFIIDPGGVIR} \code{GYSIFSYATK} \code{LFLQFGAQGSPFLK} \code{LGGNEQVTR} \code{NLSVEDAAR} \code{SADFTNFDPR} \code{TAENFR} \code{TPVISGGPYEYR} \code{TPVITGAPYEYR} \code{VEATFGVDES#N/AK} \code{VLDALQAIK} \code{VSFELFADK} \code{YILAGVENSK}}
    \item{\code{Isotope.Label.Type}}{a factor with levels \code{heavy} \code{light}}
    \item{\code{Precursor.Mz}}{a numeric vector}
    \item{\code{Precursor.Charge}}{a numeric vector}
    \item{\code{Product.Mz}}{a numeric vector}
    \item{\code{Product.Charge}}{a numeric vector}
    \item{\code{Fragment.Ion}}{a factor with levels \code{b8} \code{y10} \code{y11} \code{y12} \code{y4} \code{y5} \code{y6} \code{y7} \code{y8} \code{y9}}
    \item{\code{Retention.Time}}{a numeric vector}
    \item{\code{Area}}{a numeric vector}
    \item{\code{Background}}{a numeric vector}
    \item{\code{relative.amount}}{a numeric vector}
    \item{\code{instrument}}{a factor with levels \code{QTRAP} \code{TRIPLETOF}}
    \item{\code{user}}{a factor with levels \code{algorithm_1} \code{algorithm_2} \code{algorithm_3} \code{user_1} \code{user_2} \code{user_3} \code{user_4} \code{user_5} \code{user_6}}
    \item{\code{attempt}}{a factor with levels \code{beginner} \code{expert} \code{legacy} \code{mProphet} \code{second_best}}
  }
}
\details{
The data were derived by a user study.
Therfore each individual was given a skyline document containing the dilution series data recorded on the  QTRAP and TripleTOF 5600. Starting from skyline legacy results, users curated the data.

In  addition, two mProphet \doi{10.1038/nmeth.1584} models (as implemented in skyline) were trained on the 8 replicate data from the same instrument.
}
\examples{
data(msqc1_userstudy)
## maybe str(msqc1_userstudy) ; plot(msqc1_userstudy) ...
}
\keyword{datasets}
