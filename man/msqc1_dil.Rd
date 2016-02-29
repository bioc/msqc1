\name{msqc1_dil}
\alias{msqc1_dil}
\docType{data}
\title{
msqc1 dilution series
}
\value{
	  A \code{data.frame} and \code{msqc1} S3 object.
}

\description{
MSQC1 3 technical replicates dilution series meassured on five mass spec devices.
}
\usage{data("msqc1_dil")}
\format{
  A data frame with 15300 observations on the following 15 variables.
  \describe{
    \item{\code{Replicate.Name}}{a factor with levels \code{20140807_01_MSQC1_1_40dil_1_140808162146} \code{20140807_02_MSQC1_1_40dil_2_140808175141} \code{20140807_02_MSQC1_1_40dil_2_140808175142} \code{20140807_03_MSQC1_1_40dil_3_140808192115} \code{20140807_04_MSQC1_1_20dil_1_140808214050} \code{20140807_05_MSQC1_1_20dil_2_140808231024} \code{20140807_06_MSQC1_1_20dil_3_140809004002} \code{20140807_07_MSQC1_1_5dil_1_140809025941} \code{20140807_08_MSQC1_1_5dil_2_140809051908} \code{20140807_09_MSQC1_1_5dil_3} \code{20140807_10_MSQC1_nodil_1} \code{20140807_11_MSQC1_nodil_2} \code{20140807_12_MSQC1_nodil_3} \code{20140807_13_MSQC1_x2_1} \code{20140807_14_MSQC1_x2_2} \code{20140807_15_MSQC1_x2_3} \code{20140807_16_MSQC1_x5_1} \code{20140807_17_MSQC1_x5_2} \code{20140807_18_MSQC1_x5_3} \code{20140808_02_MSQC1_1_40dil_2_140808175141} \code{20140812_004_QC_dil1_40_1} \code{20140812_005_QC_dil1_40_2} \code{20140812_006_QC_dil1_40_3} \code{20140812_008_QC_dil1_20_1} \code{20140812_009_QC_dil1_20_2} \code{20140812_010_QC_dil1_20_3} \code{20140812_012_QC_dil1_5_1_re} \code{20140812_013_QC_dil1_5_2_re} \code{20140812_014_QC_dil1_5_3} \code{20140812_016_QC_nodil_1} \code{20140812_017_QC_nodil_2} \code{20140812_018_QC_nodil_3} \code{20140812_020_QC_conc2_1} \code{20140812_021_QC_conc2_2} \code{20140812_022_QC_conc_2_3} \code{20140812_024_QC_conc5_1} \code{20140812_025_QC_conc5_2} \code{20140812_026_QC_conc5_3} \code{20140818_003_MSQC1_1_40dil_1} \code{20140818_004_MSQC1_1_40dil_1} \code{20140818_005_MSQC1_1_40dil_3} \code{20140818_007_MSQC1_1_20dil_1} \code{20140818_008_MSQC1_1_20dil_2} \code{20140818_009_MSQC1_1_20dil_3} \code{20140818_011_MSQC1_NoDil_1} \code{20140818_013_MSQC1_NoDil_2} \code{20140818_015_MSQC1_NoDil_3} \code{20140818_017_MSQC1_x2_1} \code{20140818_019_MSQC1_x2_2} \code{20140818_021_MSQC1_x2_3} \code{20140818_023_MSQC1_x5_1} \code{20140818_025_MSQC1_x5_2} \code{20140818_027_MSQC1_x5_3} \code{20140818_030_MSQC1_1_5dil_1} \code{20140818_032_MSQC1_1_5dil_2} \code{20140818_034_MSQC1_1_5dil_3} \code{20150526_16_MSCQ1_40xdil_rep1} \code{20150526_17_MSCQ1_40xdil_rep2} \code{20150526_18_MSCQ1_40xdil_rep3} \code{20150526_20_MSCQ1_20xdil_rep1} \code{20150526_21_MSCQ1_20xdil_rep2} \code{20150526_22_MSCQ1_20xdil_rep3} \code{20150526_24_MSCQ1_5xdil_rep1} \code{20150526_25_MSCQ1_5xdil_rep2} \code{20150526_26_MSCQ1_5xdil_rep3} \code{20150526_28_MSCQ1_nodil_rep1} \code{20150526_29_MSCQ1_nodil_rep2} \code{20150526_30_MSCQ1_nodil_rep3} \code{20150526_32_MSCQ1_2x_rep1} \code{20150526_33_MSCQ1_2x_rep2} \code{20150526_34_MSCQ1_2x_rep3} \code{20150526_36_MSCQ1_5x_rep1} \code{20150526_37_MSCQ1_5x_rep2} \code{20150526_38_MSCQ1_5x_rep3} \code{20150601_001_QCdil_40dil_1} \code{20150601_002_QCdil_40dil_2} \code{20150601_003_QCdil_40dil_3} \code{20150601_004_QCdil_20dil_1} \code{20150601_005_QCdil_20dil_2} \code{20150601_006_QCdil_20dil_3} \code{20150601_007_QCdil_5dil_1} \code{20150601_008_QCdil_5dil_2} \code{20150601_009_QCdil_5dil_3} \code{20150601_010_nodil_1} \code{20150601_011_nodil_2} \code{20150601_012_nodil_3} \code{20150601_013_2conc_1_re} \code{20150601_014_2conc_2_re} \code{20150601_015_2conc_3} \code{20150601_016_5conc_1} \code{20150601_017_5conc_2} \code{20150601_018_5conc_3}}
    \item{\code{File.Name}}{a factor with levels \code{20140807_01_MSQC1_1_40dil_1_140808162146.raw} \code{20140807_02_MSQC1_1_40dil_2_140808175141.raw} \code{20140807_03_MSQC1_1_40dil_3_140808192115.raw} \code{20140807_04_MSQC1_1_20dil_1_140808214050.raw} \code{20140807_05_MSQC1_1_20dil_2_140808231024.raw} \code{20140807_06_MSQC1_1_20dil_3_140809004002.raw} \code{20140807_07_MSQC1_1_5dil_1_140809025941.raw} \code{20140807_08_MSQC1_1_5dil_2_140809051908.raw} \code{20140807_09_MSQC1_1_5dil_3.raw} \code{20140807_10_MSQC1_nodil_1.raw} \code{20140807_11_MSQC1_nodil_2.raw} \code{20140807_12_MSQC1_nodil_3.raw} \code{20140807_13_MSQC1_x2_1.raw} \code{20140807_14_MSQC1_x2_2.raw} \code{20140807_15_MSQC1_x2_3.raw} \code{20140807_16_MSQC1_x5_1.raw} \code{20140807_17_MSQC1_x5_2.raw} \code{20140807_18_MSQC1_x5_3.raw} \code{20140812_004_QC_dil1_40_1.raw} \code{20140812_005_QC_dil1_40_2.raw} \code{20140812_006_QC_dil1_40_3.raw} \code{20140812_008_QC_dil1_20_1.raw} \code{20140812_009_QC_dil1_20_2.raw} \code{20140812_010_QC_dil1_20_3.raw} \code{20140812_012_QC_dil1_5_1_re.raw} \code{20140812_013_QC_dil1_5_2_re.raw} \code{20140812_014_QC_dil1_5_3.raw} \code{20140812_016_QC_nodil_1.raw} \code{20140812_017_QC_nodil_2.raw} \code{20140812_018_QC_nodil_3.raw} \code{20140812_020_QC_conc2_1.raw} \code{20140812_021_QC_conc2_2.raw} \code{20140812_022_QC_conc_2_3.raw} \code{20140812_024_QC_conc5_1.raw} \code{20140812_025_QC_conc5_2.raw} \code{20140812_026_QC_conc5_3.raw} \code{20140818_003_MSQC1_1_40dil_1.wiff} \code{20140818_004_MSQC1_1_40dil_1.wiff} \code{20140818_005_MSQC1_1_40dil_3.wiff} \code{20140818_007_MSQC1_1_20dil_1.wiff} \code{20140818_008_MSQC1_1_20dil_2.wiff} \code{20140818_009_MSQC1_1_20dil_3.wiff} \code{20140818_011_MSQC1_NoDil_1.wiff} \code{20140818_013_MSQC1_NoDil_2.wiff} \code{20140818_015_MSQC1_NoDil_3.wiff} \code{20140818_017_MSQC1_x2_1.wiff} \code{20140818_019_MSQC1_x2_2.wiff} \code{20140818_021_MSQC1_x2_3.wiff} \code{20140818_023_MSQC1_x5_1.wiff} \code{20140818_025_MSQC1_x5_2.wiff} \code{20140818_027_MSQC1_x5_3.wiff} \code{20140818_030_MSQC1_1_5dil_1.wiff} \code{20140818_032_MSQC1_1_5dil_2.wiff} \code{20140818_034_MSQC1_1_5dil_3.wiff} \code{20150526_16_MSCQ1_40xdil_rep1.raw} \code{20150526_17_MSCQ1_40xdil_rep2.raw} \code{20150526_18_MSCQ1_40xdil_rep3.raw} \code{20150526_20_MSCQ1_20xdil_rep1.raw} \code{20150526_21_MSCQ1_20xdil_rep2.raw} \code{20150526_22_MSCQ1_20xdil_rep3.raw} \code{20150526_24_MSCQ1_5xdil_rep1.raw} \code{20150526_25_MSCQ1_5xdil_rep2.raw} \code{20150526_26_MSCQ1_5xdil_rep3.raw} \code{20150526_28_MSCQ1_nodil_rep1.raw} \code{20150526_29_MSCQ1_nodil_rep2.raw} \code{20150526_30_MSCQ1_nodil_rep3.raw} \code{20150526_32_MSCQ1_2x_rep1.raw} \code{20150526_33_MSCQ1_2x_rep2.raw} \code{20150526_34_MSCQ1_2x_rep3.raw} \code{20150526_36_MSCQ1_5x_rep1.raw} \code{20150526_37_MSCQ1_5x_rep2.raw} \code{20150526_38_MSCQ1_5x_rep3.raw} \code{20150601_001_QCdil_40dil_1.wiff} \code{20150601_002_QCdil_40dil_2.wiff} \code{20150601_003_QCdil_40dil_3.wiff} \code{20150601_004_QCdil_20dil_1.wiff} \code{20150601_005_QCdil_20dil_2.wiff} \code{20150601_006_QCdil_20dil_3.wiff} \code{20150601_007_QCdil_5dil_1.wiff} \code{20150601_008_QCdil_5dil_2.wiff} \code{20150601_009_QCdil_5dil_3.wiff} \code{20150601_010_nodil_1.wiff} \code{20150601_011_nodil_2.wiff} \code{20150601_012_nodil_3.wiff} \code{20150601_013_2conc_1_re.wiff} \code{20150601_014_2conc_2_re.wiff} \code{20150601_015_2conc_3.wiff} \code{20150601_016_5conc_1.wiff} \code{20150601_017_5conc_2.wiff} \code{20150601_018_5conc_3.wiff}}
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
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(msqc1_dil)
## maybe str(msqc1_dil) ; plot(msqc1_dil) ...
}
\keyword{datasets}
