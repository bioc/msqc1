#R

# $HeadURL$
# $Id$
# $Date$

.signal_response_ident_plot <- function(data, 
                           prot=c("P00915", "P00918", "P04040", "P02741", "P15559", "PPIA"), 
                           instrument="NA", peptides){
  
  x.peptide <- aggregate(Area ~ Peptide.Sequence + Isotope.Label.Type + File.Name.Id + Protein.Name, 
                         data=data, 
                         FUN=function(x){
                           sum(as.numeric(x), na.rm = TRUE)
                           }
                         )
  
  r <- lapply(prot, function(p){
    x.peptide.filtered <- droplevels(x.peptide[grepl(p, x.peptide$Protein.Name), ])
    x.a <- aggregate(Area ~ Peptide.Sequence, FUN=max, data=x.peptide.filtered)
    
    idx <- order(x.a[,2], decreasing = TRUE)
    
    x.peptide.filtered$Peptide.Sequence <- factor(x.peptide.filtered$Peptide.Sequence, 
                                                  levels = x.a[order(x.a[,2], 
                                                                     decreasing = TRUE), 1])
    
    xyplot(log(Area,10) ~ Peptide.Sequence | p * instrument, 
                 groups=x.peptide.filtered$Protein.Name, 
                 data=x.peptide.filtered,
                 scales = list(x = list(rot = 45)), 
                 panel = function(x,y,...){
                   panel.abline(v = x[x %in% peptides$Peptide.Sequence], col='grey', lwd=3)
                   panel.xyplot(x,y,...)
                 },
                 auto.key=list(space = "top", points = TRUE, lines = FALSE, cex=1),
                 levels=order(x.peptide.filtered$Peptide.Sequence, 
                              decreasing = TRUE)
    )
  })
}

.get_sensitivity_relative.amount <- function(x, stepsize=0.005, max=3 * 13, log2ratio.min=0, log2ratio.max=8){
  log2ratio.measured <- log(x$Area.light,2 ) - log(x$Area.heavy,2)
  log2ratio.theoretic <- unlist(lapply(x$Peptide.Sequence, function(p){
    idx <- which(as.character(p) == as.character(msqc1::peptides$Peptide.Sequence)); 
    log2(msqc1::peptides$actual.LH.ratio[idx])
  }
  ))
  
  rr <- lapply(unique(x$relative.amount), function(relative.amount){
    r <-lapply(unique(x$instrument), function(instrument){
      eps <- seq(log2ratio.min, log2ratio.max, by=stepsize)
      
      tt <-unlist(lapply(eps, function(eps){
        filter <- (x$instrument == instrument & x$relative.amount == relative.amount)
        sum(abs(log2ratio.measured[filter] - log2ratio.theoretic[filter]) < eps) / max 
      }))
      
      data.frame(log2ratio.cutoff=eps, sensitivity=tt, instrument=instrument, relative.amount=relative.amount)
    })
    
    do.call('rbind', r)
  })
  do.call('rbind', rr)
 }

.get_sensitivity <-function(x, stepsize=0.005, 
                            max = nrow(x) / length(unique(x$instrument)), 
                            log2ratio.min=0, log2ratio.max=8, lh.aggregate=FALSE){
 
  if (lh.aggregate){
    
    xx.heavy <- aggregate(Area.heavy ~ instrument + Peptide.Sequence, data=x, FUN=mean)
    xx.light <- aggregate(Area.light ~ instrument + Peptide.Sequence, data=x, FUN=mean)
    x <- cbind(xx.light, xx.heavy$Area.heavy)
    names(x)[4] <- "Area.heavy"
  }
  
  
  message(paste("number of elements is", nrow(x)))
  log2ratio.measured <- log(x$Area.light, 2) - log(x$Area.heavy,2)
  
  log2ratio.theoretic <- unlist(lapply(x$Peptide.Sequence, function(p){
    idx <- which(as.character(p) == as.character(msqc1::peptides$Peptide.Sequence)); 
    log2(msqc1::peptides$actual.LH.ratio[idx])
    }
    ))
  
  r <-lapply(unique(x$instrument), function(instrument){
    eps <- seq(log2ratio.min, log2ratio.max, by=stepsize)
    
    tt <-unlist(lapply(eps, function(eps){
      sum(abs(log2ratio.measured[x$instrument == instrument] - log2ratio.theoretic[x$instrument == instrument]) < eps) / max 
      }))
    
    #plot(eps,tt,type='b'); abline(v=1)
    data.frame(log2ratio.cutoff=eps, sensitivity=tt, instrument=instrument)
  })
  do.call('rbind', r)
  #xyplot(tt~eps, group=instrument, data=do.call('rbind', r),type='b', auto.key=list(space = "top", points = TRUE, lines = FALSE, cex=1))
}

.panel.msqc1_dil <- function(...){
  idx <- order(msqc1::peptides$SIL.peptide.per.vial[order(msqc1::peptides$Peptide.Sequence)])
  peptide.idx <- sort(msqc1::peptides$Peptide.Sequence)[idx]
  
  pep <- peptide.idx[panel.number()]
  pep.idx <- (which(as.character(pep) == as.character(msqc1::peptides$Peptide.Sequence)))
  Protein.Name <- (msqc1::peptides$Protein.Name[pep.idx])
  actual.LH.ratio <- (msqc1::peptides$actual.LH.ratio[pep.idx])
  
  SIL.peptide.per.vial <- (msqc1::peptides$SIL.peptide.per.vial[(which(as.character(pep) == as.character(msqc1::peptides$Peptide.Sequence)))])
  
  # log2 changes
  panel.rect(-100,log2(actual.LH.ratio)-1, 5, log2(actual.LH.ratio)+1 ,border = '#EEEEEEEE', col='#EEEEEEEE')
  panel.rect(-100,log2(actual.LH.ratio)-0.5, 5, log2(actual.LH.ratio)+0.5 ,border = '#CCCCCCCC', col='#CCCCCCCC')
  
  # zero line
  panel.abline(h=c(0, log2(actual.LH.ratio)), col=c("grey", "black"), lwd=c(1,1))
  
  # plot the data
  panel.xyplot(...)
  panel.xyplot(..., type='smooth')
  
  # legend
  SIL.onColumnAmount <- paste("SIL on column=", round(SIL.peptide.per.vial * 10, 1), "fmol",sep='')
  message(Protein.Name)
  x.x<-0.85
  x.cex<-0.7
  ltext(x=x.x, y=8.6, as.character(Protein.Name), cex=x.cex, pos=2)
  ltext(x=x.x, y=8, SIL.onColumnAmount, cex=x.cex, pos=2)
  #ltext(x=x.x, y=8, paste("SIL=",round(SIL.peptide.per.vial, 2),sep=''), cex=x.cex, pos=2)
  ltext(x=x.x, y=7.4, paste("actual.LH.ratio=",round(actual.LH.ratio, 2),sep=''), cex=x.cex, pos=2)
}

.panel.8rep <- function(x, y, ..., errorflag=TRUE) {
  
  idx <- order(msqc1::peptides$SIL.peptide.per.vial[order(msqc1::peptides$Peptide.Sequence)])
  peptide.idx <- sort(msqc1::peptides$Peptide.Sequence)[idx]
  
  pepStr <- peptide.idx[panel.number()] 
  
  TheoRatioLH <- log2(msqc1::peptides$actual.LH.ratio[(which(as.character(pepStr) == as.character(msqc1::peptides$Peptide.Sequence)))])
  SIL.peptide.per.vial <- msqc1::peptides$SIL.peptide.per.vial[(which(as.character(pepStr) == as.character(msqc1::peptides$Peptide.Sequence)))]
  pep.idx <- (which(as.character(pepStr) == as.character(msqc1::peptides$Peptide.Sequence)))
  Protein.Name <- (msqc1::peptides$Protein.Name[pep.idx])
  
  if (errorflag){
   panel.abline(h=0, lwd=2, col='grey')
###################
   a.mean <- aggregate((y - TheoRatioLH) ~ x, FUN = mean)
   
   a.sd <- aggregate((y - TheoRatioLH) ~ x, FUN = sd)
###################
   
   #panel.segments(a.mean[,1], a.mean[,2] - a.sd[,2], a.mean[,1], a.mean[,2] + a.sd[,2],  lwd=3, col='black')
   
   #panel.points(a.mean[,1], a.mean[,2] - a.sd[,2], pch="-", col='black', cex=3)
   #panel.points(a.mean[,1], a.mean[,2] + a.sd[,2], pch="-", col='black', cex=3)

   #panel.points(a.mean[,1], a.mean[,2], col='cornflowerblue', pch='.', cex=3)
   panel.bwplot(x,y- TheoRatioLH,...)
  }else{
  
   panel.rect(0, TheoRatioLH - 1, 6, TheoRatioLH + 1, border = '#EEEEEEEE', col='#EEEEEEEE')
   panel.rect(0, TheoRatioLH - 0.5, 6, TheoRatioLH + 0.5, border = '#CCCCCCCC', col='#CCCCCCCC')
   panel.abline(h=TheoRatioLH, col="black", lwd=1)
  
   trust.ratios <- TheoRatioLH + (c(-0.1, 0.1) * TheoRatioLH)
   panel.abline(h=0, lwd=2, col='grey')
  
   panel.xyplot(x, y, ...)
  
   SILonColumnAmount <- paste("SIL on column=", round(SIL.peptide.per.vial * 10, 1), "fmol",sep='')
   x.x <- 5.5
   x.y <-5
  
   if (TheoRatioLH > 3.5){x.y<-1}
    x.cex<-0.7
    ltext(x=x.x, y=x.y+.6, as.character(Protein.Name), cex=x.cex, pos=2)
    ltext(x=x.x, y=x.y, SILonColumnAmount, cex=x.cex, pos=2)
  }
}


.load_msqc1 <- function(name=''){
  f <- system.file("extdata", name, package = "msqc1")
  if (name == ''){
    print(list.files(path = f))
  }else{
  
  s <- read.csv(f, sep=",", 
                header=TRUE, 
                na.strings='#N/A')
  class(s)=c('msqc1',class(s))
  return (s)}
}

.compute_cvs <- function(s, ...){
  s.light <- s[s$Isotope.Label.Type == "light",]
  s.heavy <- s[s$Isotope.Label.Type == "heavy",]
  m <- merge(s.heavy, s.light, by=c('instrument', 'File.Name', 'Peptide.Sequence', 'Fragment.Ion', 'File.Name.Id'))
  
  m.a <- aggregate(cbind(heavy=log2(Area.x), light=log2(Area.y)) ~ instrument + Peptide.Sequence + Fragment.Ion, 
                   data=m, FUN=function(x){print (paste(x[1,], x[2,]))})
  
  print(xyplot(Fragment.Ion ~ m.a[,4] | instrument * Peptide.Sequence , data=m.a, xlab="coefficient of variation", 
               scales = list(y = list(rot = 45))))
  print(bwplot(~m.a[,4] | instrument * Peptide.Sequence, data=m.a, xlab="coefficient of variation"))
  print(bwplot(~m.a[,4] | instrument, data=m.a, xlab="coefficient of variation", layout=c(1,5)))
  
}


## plots used in the publication
.figure_setup <- function(){
  
  t <- trellis.par.get()
  t$superpose.symbol$pch <- 1:7
  t$superpose.line$col <- c("#0080ffAA", "#ff00ffAA",
  	"#00640080", "#ff0000AA", "#FFA50080", 
	"#00ff00AA", "#A52A2A80")
  
  yg16 <- c("#EDF8E9","#E0F3DA", "#D2EDCC", "#C4E8BE",
  	"#B5E2AF", "#A3D99E", "#91D18E", "#7EC87E",
	"#6CC071", "#5CB768", "#4AAE5F", "#37A556",
	"#29984C", "#1D8941", "#107B36", "#006D2C")
  
  t$regions$col <- yg16
  
  ot <- trellis.par.set(t)
  
}

.figure1 <- function(data, peptides){
  .figure_setup()
  
  s.8rep = data
  # consider only b and y (e.g., no MS1!!!)
  s.8rep <- s.8rep[grep("[by]", s.8rep$Fragment.Ion), ]

  s.8rep.peptide <- aggregate(Area ~ Peptide.Sequence + Isotope.Label.Type + instrument + File.Name.Id, 
                              data=s.8rep, 
                              FUN=sum)
  
  s.8rep.peptide.light <- s.8rep.peptide[s.8rep.peptide$Isotope.Label.Type == "light",]
  s.8rep.peptide.heavy <- s.8rep.peptide[s.8rep.peptide$Isotope.Label.Type == "heavy",]
  m.8rep.peptide <- merge(s.8rep.peptide.heavy, s.8rep.peptide.light, 
                          by=c('instrument', 'Peptide.Sequence', 'File.Name.Id'))
  
  names(m.8rep.peptide) <- c("instrument", "Peptide.Sequence", "File.Name.Id", 
                             "Isotope.Label.Type.x", "Area.heavy", 
                             "Isotope.Label.Type.y", "Area.light")
  
  idx <- order(peptides$SIL.peptide.per.vial[order(peptides$Peptide.Sequence)])
  peptide.idx <- sort(peptides$Peptide.Sequence)[idx]
  
  xyplot((log2(Area.light) - log2(Area.heavy)) ~ instrument | Peptide.Sequence, 
         index.cond = list(idx),
         group = m.8rep.peptide$instrument,
         panel = function(...){.panel.8rep(..., errorflag=FALSE)},
         auto.key = list(space = "top", points = TRUE, lines = FALSE, cex=1),
         scales = list(x = list(rot = 45)),
         layout = c(7,2),
         sub = "log2 light heavy ratios of 8 replicates on 5 MS plattforms",
         data=m.8rep.peptide)
}

.figure2 <- function(data){
  .figure_setup()
  #s.8rep <- load_msqc1("8Rep.csv")
  s.8rep <- data
  s.8rep <- s.8rep[grep("[by]", s.8rep$Fragment.Ion), ]
  s.8rep.peptide <- aggregate(Area ~ Peptide.Sequence + Isotope.Label.Type + instrument + File.Name.Id, 
                              data=s.8rep, FUN=sum)
  
  
  s.8rep.peptide.light <- s.8rep.peptide[s.8rep.peptide$Isotope.Label.Type == "light",]
  s.8rep.peptide.heavy <- s.8rep.peptide[s.8rep.peptide$Isotope.Label.Type == "heavy",]
  m.8rep.peptide <- merge(s.8rep.peptide.heavy, s.8rep.peptide.light, 
                          by=c('instrument', 'Peptide.Sequence', 'File.Name.Id'))
  
  names(m.8rep.peptide) <- c("instrument", "Peptide.Sequence", "File.Name.Id", 
                             "Isotope.Label.Type.x", "Area.heavy", 
                             "Isotope.Label.Type.y", "Area.light")
  
  xx.table <- table(s.8rep$Peptide.Sequence)
  
  sensitivity.data <- m.8rep.peptide[ !m.8rep.peptide$Peptide.Sequence  %in%  c('TAENFR','GYSIFSYATK'), ]
  sensitivites <- .get_sensitivity(sensitivity.data, 
                                           log2ratio.max = 1,
                                           lh.aggregate = FALSE,
                                           max = max(table(sensitivity.data$instrument)))
  AUC <- aggregate(sensitivity ~ instrument, 
                   data=sensitivites, 
                   FUN=function(x){round(sum(x)/length(x),2)})
  
  xyplot(sensitivity ~ log2ratio.cutoff, 
         group=sensitivites$instrument, 
         #xlab=paste("log2-scaled L:H cutoff value", expression(epsilon)),
         xlab=expression(epsilon),
         ylab='relative amount correctly quantified replicates',
         #sub=list(paste("not considered peptides: ",peptides.not.considered ,sep=''), cex=0.5),
         panel=function(x,y,...){
           
           panel.rect(0,0,0.5,1,col='#CCCCCCCC', border='#CCCCCCCC')
           panel.rect(0.5,0,1,1,col='#EEEEEEEE', border='#EEEEEEEE')
           panel.abline(h=0.9)
           panel.xyplot(x,y,...)
           panel.text(0.50,0.55,"AUC:", pos=4)
           panel.text(0.51, seq(0.3,0.5,length=5), paste(AUC$instrument, AUC$sensitivity, sep=" = "), pos=4, cex=0.75)
         },
         #ylab = list('relative number of data items having a distance to the theoretical log2ratio cutoff',  cex=.75),
         data=sensitivites,
         type = 'l', 
         auto.key=list(space = "top", 
                       points = TRUE, 
                       lines = FALSE, cex=1), 
         #sub=as.character(paste(paste(AUC$instrument, AUC$sensitivity, sep="="), col=',',sep='')),
         lwd = 3)
}

.figure3 <- function(data, peptides){
  .figure_setup()

  x = data
  x <- x[grepl("[by]", x$Fragment.Ion) & x$Peptide.Sequence %in% peptides$Peptide.Sequence, ]
  
  x.sum <- aggregate(Area ~ instrument + Isotope.Label.Type + relative.amount + Peptide.Sequence + File.Name.Id, 
                     data=x, FUN=sum)
  
  x.sum.cv <-  aggregate(Area ~ instrument + Isotope.Label.Type + Peptide.Sequence, 
                         data = x.sum, 
                         FUN = function(x){100 * sd(x) / mean(x)})

  bwplot(Area ~ instrument,
         data = x.sum.cv, 
         panel = function(...){
           panel.abline(h = log(c(5,10,20), base = 10), col.line ='grey')
           panel.violin(..., fill = NULL, col = NULL, adjust = 1.0, varwidth = FALSE)
           panel.bwplot(..., fill = "#AAAAAA88")
         },
         group=x.sum.cv$Isotope.Label.Type, 
         scales = list(x = list(rot = 45), 
                       y = list(log = TRUE, at = c(1, 2, 5, 10, 20, 50 ,100))),
         ylab = 'coefficient of variation [%]'
  )
  
}

.figure4 <- function(data, peptides){
  .figure_setup()
  # load the package data
  # s <- msqc1::load_msqc1('dilSeries.csv')
  
  s <- data
  s <- s[grep("[by]", s$Fragment.Ion), ]
  # aggregate the haevy and light transition areas
  s.peptide_areas <- aggregate(Area ~ instrument + Isotope.Label.Type + relative.amount + Peptide.Sequence + File.Name, 
                               data=s, FUN=function(x){sum(x, na.rm=TRUE)})
  
  s.light <- s.peptide_areas[s.peptide_areas$Isotope.Label.Type=='light',] 
  s.heavy <- s.peptide_areas[s.peptide_areas$Isotope.Label.Type=='heavy',]
  
  s.peptie_areas_hl <- merge(s.heavy, s.light, by=c('instrument', 'Peptide.Sequence', 'relative.amount', 'File.Name'))
  
  names(s.peptie_areas_hl) <- c("instrument", "Peptide.Sequence", 
                                "relative.amount", "File.Name", 
                                "Isotope.Label.Type.x", "Area.heavy", 
                                "Isotope.Label.Type.y", "Area.light")
  
  idx <- order(peptides$SIL.peptide.per.vial[order(peptides$Peptide.Sequence)])
  peptide.idx <- sort(peptides$Peptide.Sequence)[idx]
  
  xyplot(log2(Area.light) - log2(Area.heavy) ~ relative.amount | Peptide.Sequence, 
         groups=s.peptie_areas_hl$instrument,
         data=s.peptie_areas_hl, 
         panel=.panel.msqc1_dil,
         layout=c(7,2),
         auto.key=list(space = "top", points = TRUE, lines = FALSE, cex=1),
         index.cond=list(idx),
         scales=list(x = list(rot = 45, log=TRUE, at=sort(unique(s.peptie_areas_hl$relative.amount)) )),
         sub="log2 light heavy ratios of 6 dilutions on 5 MS plattforms",
         main='sigma mix peptide level signal')
}

.figure5 <- function(data, data_reference){
  .figure_setup()
  
  s <- data
  s <- s[grep("[by]", s$Fragment.Ion), ]
  s <- s[!s$Peptide.Sequence %in%  c('TAENFR','GYSIFSYATK'), ]
  
  s.peptide_areas <- aggregate(Area ~ instrument + Isotope.Label.Type + relative.amount + Peptide.Sequence + File.Name,  data=s, FUN=function(x){sum(x, na.rm=TRUE)})
  
  s.light <- s.peptide_areas[s.peptide_areas$Isotope.Label.Type=='light',] 
  s.heavy <- s.peptide_areas[s.peptide_areas$Isotope.Label.Type=='heavy',]
  
  s.peptie_areas_hl <- merge(s.heavy, s.light, by=c('instrument', 'Peptide.Sequence', 'relative.amount', 'File.Name'))
  
  names(s.peptie_areas_hl) <- c("instrument", "Peptide.Sequence", 
                                "relative.amount", "File.Name", 
                                "Isotope.Label.Type.x", "Area.heavy", 
                                "Isotope.Label.Type.y", "Area.light")
  
 
  s.8rep <- data_reference
  xx.table <- table(s.8rep$Peptide.Sequence)
  
  
  xxx <- do.call('rbind', 
                 lapply(c(0), function(cutoff){
                   xx.240 <- names(xx.table[xx.table > cutoff * 240 | xx.table == 160])
                   
                   sensitivity.dil <- .get_sensitivity_relative.amount(s.peptie_areas_hl[s.peptie_areas_hl$Peptide.Sequence %in% xx.240,], 
                                                                               max=42/14*12)
                   
                   sensitivity.dil$cutoff <- cutoff
                   sensitivity.dil
                 }))
                 
  
  AUC <- aggregate(sensitivity ~ instrument + relative.amount, 
                   data=xxx[xxx$log2ratio.cutoff<=1,], 
                   FUN=function(x){round(sum(x)/length(x),2)})
  
  AUC.relative.amount <- unique(AUC$relative.amount)
  
  xyplot(sensitivity ~ log2ratio.cutoff | relative.amount , 
         group=xxx$instrument, 
         #xlab="log2-scaled L:H cutoff value ",
         xlab=expression(epsilon),
         ylab='relative amount correctly quantified replicates',
         panel=function(...){
           pn = panel.number()
           panel.rect(0,0,0.5,1,col='#CCCCCCCC', border='#CCCCCCCC')
           panel.rect(0.5,0,1,1,col='#EEEEEEEE', border='#EEEEEEEE')
           panel.abline(h=c(0.5,0.9,1), col='darkgrey')
           panel.xyplot(...)
           message(paste(pn, AUC.relative.amount[pn]))
           AUC.panel <- AUC[AUC$relative.amount == AUC.relative.amount[pn], ]
           panel.text(1.50,0.31,"AUC [0,1]:", pos=4)
           panel.text(1.51, seq(0.05,0.25,length=5), 
                      paste(AUC.panel$instrument, AUC.panel$sensitivity, sep=" = "), 
                      pos=4, cex=0.75)
         },
         sub = list('sensitivity (relative number of data items having a distance to the theoretical log2ratio cutoff)', cex=1),
         data=xxx,
         xlim=c(0,5),
         type='l', 
         strip = strip.custom(strip.names = TRUE, strip.levels = TRUE),
         auto.key=list(space = "top", points = TRUE, lines = FALSE, cex=1),
         lwd=3
  )
}

.figure6 <- function(data, peptides){
  .figure_setup()
  x_user <- data
  x_user <- droplevels(x_user[x_user$Peptide.Sequence %in% peptides$Peptide.Sequence, ])
  
  x_user.light <- x_user[x_user$Isotope.Label.Type == "light",]
  x_user.heavy <- x_user[x_user$Isotope.Label.Type == "heavy",]
  
  m_user <- merge(x_user.heavy, x_user.light, 
                  by=c('instrument', 'File.Name', 'Peptide.Sequence', 'Replicate.Name',
                       'Fragment.Ion', 'relative.amount', 'user', 'attempt', 'Protein.Name', 'Precursor.Charge'))
  
  pp <- data.frame(ratio=m_user$Area.y / m_user$Area.x, 
                   user=m_user$user, 
                   attempt=m_user$attempt,
                   instrument=m_user$instrument, 
                   Peptide.Sequence=m_user$Peptide.Sequence)
  
  m_user <- merge(peptides, pp, by='Peptide.Sequence')
  
  idx <- which(is.infinite(m_user$ratio))
  m_user$ratio[idx] <- NA
  
  user_n <- aggregate((actual.LH.ratio - ratio) ~ user + instrument + attempt, 
                      data=m_user, 
                      FUN=function(x){length(x)})
  names(user_n)[4] <- 'n'
  user_sd <- aggregate((actual.LH.ratio - ratio) ~ user + instrument + attempt, 
                       data=m_user, 
                       FUN=function(x){sd(x,na.rm=TRUE)})
  names(user_sd)[4] <- 'sd'
  
  user_sd_n <- cbind(user_sd, user_n)
  xyplot(sd ~ n | instrument, 
         group=user_sd_n$attempt, 
         data=user_sd_n,
         xlab='number of valid ratios',
         ylab="Standard deviation of Error",
         panel = function(...){
           
           auto.sd <- user_sd$sd[user_sd$attempt=='legacy']
           auto.n <- user_n$n[user_n$attempt=='legacy']
           panel.abline(h=auto.sd[panel.number()],col='grey')
           panel.abline(v=auto.n[panel.number()],col='grey')
           panel.xyplot(..., cex=1.4, lwd=1.4, type='p')
         },
         auto.key=list(space = "right", points = TRUE, lines = FALSE, cex=1.5)
  )
}


.get_log2LHratio <- function(s){
  s.peptide_areas <- aggregate(Area ~ instrument + Isotope.Label.Type + relative.amount + Peptide.Sequence + File.Name, 
                               data=s, FUN=function(x){sum(x, na.rm=TRUE)})
  
  s.light <- s.peptide_areas[s.peptide_areas$Isotope.Label.Type=='light',] 
  s.heavy <- s.peptide_areas[s.peptide_areas$Isotope.Label.Type=='heavy',]
  
  s.peptie_areas_hl <- merge(s.heavy, s.light, by=c('instrument', 'Peptide.Sequence', 'relative.amount', 'File.Name'))
  
  names(s.peptie_areas_hl) <- c("instrument", "Peptide.Sequence", 
                                "relative.amount", "File.Name", 
                                "Isotope.Label.Type.x", "Area.heavy", 
                                "Isotope.Label.Type.y", "Area.light")
  
  return(droplevels(s.peptie_areas_hl))
}


.supplement_figure6 <- function(data, peptides){
  .figure_setup()
  s <- data
  s <- s[s$instrument %in% c('QExactive', 'QExactiveHF', 'TRIPLETOF'), ]
  s.product.ions <- s[grep("[by]", s$Fragment.Ion), ]
  s.ms1 <- s[!grepl("^[by]", s$Fragment.Ion), ]
  
  s.product.ions.log2LHratio <- .get_log2LHratio(s.product.ions)
  s.product.ions.log2LHratio$ms.level <- as.factor("product.ions")
  
  s.ms1.log2LHratio <- .get_log2LHratio(s.ms1)
  s.ms1.log2LHratio$ms.level <- as.factor("MS1")
  
  ##########################################################################
  pp <- msqc1::peptides[msqc1::peptides$Peptide.Sequence!="TAENFR", ]
  
  idx <- order(pp$SIL.peptide.per.vial[order(pp$Peptide.Sequence)])
  peptide.idx <- sort(pp$Peptide.Sequence)[idx]
  
  xx<-droplevels(do.call('rbind', list(s.product.ions.log2LHratio, s.ms1.log2LHratio)))
  
  xyplot(log2(Area.light)- log2(Area.heavy) ~ relative.amount |  Peptide.Sequence * ms.level,
         index.cond=list(idx, c(2,1)),
         groups=xx$instrument, 
         data = xx, 
         scales=list(x = list(rot = 45, 
                              log=TRUE, 
                              at=sort(unique(xx$relative.amount)) )), 
         type='p', 
         auto.key=list(space = "top", points = TRUE, lines = FALSE, cex=1),
         panel=function(...){  
           pep <- peptide.idx[((panel.number() - 1) %% length(peptide.idx)) + 1]
           pep.idx <- (which(as.character(pep) == as.character(peptides$Peptide.Sequence)))
           
           Protein.Name <- peptides$Protein.Name[pep.idx]
           
           theoretic.log2.lh.ratio <- log2(msqc1::peptides$actual.LH.ratio[pep.idx])
           
           SIL.peptide.per.vial <- (peptides$SIL.peptide.per.vial[(which(as.character(pep) == as.character(peptides$Peptide.Sequence)))])
           
           # log2 changes
           panel.rect(-100,theoretic.log2.lh.ratio - 1, 5, theoretic.log2.lh.ratio + 1,
                      border = '#EEEEEEEE', col='#EEEEEEEE')
           panel.rect(-100,theoretic.log2.lh.ratio - 0.5, 5, theoretic.log2.lh.ratio+0.5,
                      border = '#CCCCCCCC', col='#CCCCCCCC')
           
           panel.abline(h=c(0, theoretic.log2.lh.ratio), col=c("grey", "black"), lwd=c(1,1))
           
           panel.xyplot(...)
           #message(panel.number())
         }
  )
  
  
}

.supplement_figure7 <- function(data){
  s <- data
  s <- s[s$instrument %in% c('QExactive', 'QExactiveHF', 'TRIPLETOF'), ]
  
  s.product.ions <- s[grep("[by]", s$Fragment.Ion), ]
  s.ms1 <- s[!grepl("^[by]", s$Fragment.Ion), ]
  
  s.product.ions.log2LHratio <- .get_log2LHratio(s.product.ions)
  s.product.ions.log2LHratio$ms.level <- as.factor("product.ions")
  
  s.ms1.log2LHratio <- .get_log2LHratio(s.ms1)
  s.ms1.log2LHratio$ms.level <- as.factor("MS1")
  
  sensitivity.product.ions <- .get_sensitivity_relative.amount(s.product.ions.log2LHratio[s.product.ions.log2LHratio$Peptide.Sequence !='GYSIFSYATK', ], max=1)
  
  sensitivity.product.ions$ms.level <- as.factor("product.ions")
  
  sensitivity.ms1 <- .get_sensitivity_relative.amount(s.ms1.log2LHratio[s.ms1.log2LHratio$Peptide.Sequence != 'GYSIFSYATK',], max=1)
  sensitivity.ms1$ms.level <- "MS1"
  
  sensitivity <- do.call('rbind', list(sensitivity.ms1, sensitivity.product.ions))
  
  
  AUC.s <- aggregate(sensitivity / max(sensitivity) ~ instrument + relative.amount + ms.level, 
                     data = sensitivity[sensitivity$log2ratio.cutoff <= 1, ], 
                     FUN = function(x){
                       round(sum(x) / length(x), 2)
                     })
  AUC.relative.amount <- rep(unique(AUC.s$relative.amount), 2)
  AUC.ms.level <- c(rep("MS1", 6), rep("product.ions", 6))
  group.instrument <- sensitivity$instrument
  xyplot(sensitivity / max(sensitivity) ~ log2ratio.cutoff | relative.amount * ms.level, 
         data=sensitivity, 
         group=group.instrument,
         type='l',
         auto.key=list(space = "top", points = TRUE, lines = FALSE, cex=1),
         xlim=c(0,5),
         panel = function(...){
           pn = panel.number()
           AUC.panel.s <- AUC.s[(AUC.s$relative.amount == AUC.relative.amount[pn] & AUC.s$ms.level == AUC.ms.level[pn]), ]
           panel.rect(0,0,0.5,1,col='#CCCCCCCC', border='#CCCCCCCC')
           panel.rect(0.5,0,1,1,col='#EEEEEEEE', border='#EEEEEEEE')
           panel.abline(h=c(0.5,0.9,1), col='darkgrey')
           panel.xyplot(...)
           panel.text(1.50,0.31,"AUC [0,1]:", pos=4)
           AUC.values <-  paste(AUC.panel.s[,1], AUC.panel.s[,4], sep=" = ")
           message(AUC.values)
           panel.text(1.51, seq(0.10,0.23, length=length(AUC.values)), 
                      paste(AUC.panel.s[,1], AUC.panel.s[,4], sep=" = "), 
                      pos=4, cex=0.75)
         },
         strip = strip.custom(strip.names = TRUE, strip.levels = TRUE),
         led=3,
         xlab=expression(epsilon),
         ylab='relative amount correctly quantified replicates'
  )
}


.supplement_figure8 <- function(data){
  .figure_setup()
  
  peptide.filter <- c('GGPFSDSYR', 'VLDALQAIK', 'EGHLSPDIVAEQK', 'ALIVLAHSER', 'ESDTSYVSLK', 'GYSIFSYATK', 'FEDENFILK',
                      'FSTVAGESGSADTVR',  'NLSVEDAAR', 'SADFTNFDPR')
  
  s <- data[data$Peptide.Sequence %in% peptide.filter, ]
  s <- s[grep("[by]", s$Fragment.Ion), ]
  
  s.agg <- aggregate(Retention.Time ~ Peptide.Sequence * File.Name.Id * instrument, FUN=mean, data=s)
  s.model <- lm(Retention.Time ~ Peptide.Sequence * instrument, data=s.agg)
  summary(s.model)
  
  s.mean <- aggregate(Retention.Time ~ Peptide.Sequence *  instrument, FUN=mean, data=s.agg)
  s.sd <- aggregate(Retention.Time ~ Peptide.Sequence *  instrument, FUN=sd, data=s.agg)
  s.m<-merge(s.mean, s.sd, by=c('Peptide.Sequence', 'instrument'))
  
  xyplot(Retention.Time.y/Retention.Time.x ~ Peptide.Sequence, 
         data=s.m, 
         type='b',
         group=s.m$instrument, 
         auto.key=list(space = "top", points = TRUE, lines = FALSE, cex=1),
         scales = list(x = list(rot = 45)), 
         ylab='coefficient of variation (CV) -  between 8 replicate retention time')
  
  xyplot(Retention.Time.y/Retention.Time.x ~ Peptide.Sequence, 
         data=s.m, 
         type='b',
         group=s.m$instrument, 
         auto.key=list(space = "top", points = TRUE, lines = FALSE, cex=1),
         scales = list(x = list(rot = 45), y = list(log = TRUE, at=c(0.001, 0.01, 0.1, 1))), 
         ylab='coefficient of variation (CV) -  between 8 replicate retention time')
  
  
  # ok; ordering is not working
  xyplot(Retention.Time ~ File.Name.Id | Peptide.Sequence, 
         group=s.agg$instrument, 
         data=s.agg, 
         auto.key=list(space = "top", points = TRUE, lines = FALSE, cex=1),
         type='b',
         scales = list(x = list(rot = 45)))
}
  
.supplement_figure9 <- function(data){
    .figure_setup()
    
    s <- data
    s <- s[grep("[by]", s$Fragment.Ion), ]
    s <- aggregate(Retention.Time ~ Peptide.Sequence * File.Name.Id * instrument, FUN=mean, data=s)
   
    tt<-table(s$Peptide.Sequence) 
    # use only the peptides which were meassured on all plattforms on all eight repititions
    peptide.filter <- names(tt[tt==40])
    s<- droplevels(s[s$Peptide.Sequence %in% peptide.filter,])
    
    s.min <- aggregate(Retention.Time ~ instrument * File.Name.Id , FUN=min, data=s)
    s.max <- aggregate(Retention.Time ~ instrument * File.Name.Id , FUN=function(x){ max(x)}, data=s)
  
    s <- merge(merge(s, s.min, by=c('File.Name.Id', 'instrument')), s.max, by=c('File.Name.Id', 'instrument'))
    s$Retention.Time.Delta <- (s$Retention.Time - s$Retention.Time.y)
    s$Retention.Time.Normalized <- (s$Retention.Time.x - s$Retention.Time.y) / s$Retention.Time.Delta
    #s$Retention.Time.Normalized <- (s$Retention.Time.x ) / s$Retention.Time
    
    xyplot(Retention.Time.Normalized ~ File.Name.Id | Peptide.Sequence, 
         group=s$instrument, 
         data=s, 
         auto.key=list(space = "top", points = TRUE, lines = FALSE, cex=1),
         type='b',
         scales = list(x = list(rot = 45)), layout=c(length(unique(s$Peptide.Sequence)),1))
}

.normalize_rt <- function(S, S.training, reference_instrument = 'Retention.Time.QTRQP'){
  S.normalized <- S
  
  # build linear model to reference_instrument 
  # TODO(cp): use by
  for (instrument in unique(S.normalized$instrument)){
    
    i <- paste("Retention.Time", instrument, sep='.')
    
    rt.out <- S.training[, reference_instrument]
    rt.in <- S.training[, i]
    S.fit <- lm(rt.out ~ rt.in)
    
    S.normalized[S.normalized$instrument == instrument, 'Retention.Time'] <- predict(S.fit,                                                                                 data.frame(rt.in = S.normalized[S.normalized$instrument == instrument, 
                                                                                                                                                                                                            'Retention.Time']))
    
  }
  
  # scaling
  S.normalized.min <- min(S.normalized$Retention.Time, na.rm = TRUE)
  S.normalized.delta <- max(S.normalized$Retention.Time, na.rm = TRUE) - S.normalized.min 
  
  S.normalized$Retention.Time <- (S.normalized$Retention.Time - S.normalized.min) / S.normalized.delta
  
  return(S.normalized)
}

.reshape_rt <- function(S, peptides=peptides, plot=TRUE){
  S <- S[grep("[by]", S$Fragment.Ion), ]
  S <- S[S$Peptide.Sequence %in% peptides$Peptide.Sequence, ]
  S <- aggregate(Retention.Time ~ Peptide.Sequence * instrument,
                 FUN=mean, 
                 data=S)
  
  #tt <- (table(S$Peptide.Sequence) == frequencyTableCutOff)
  #S <- S[S$Peptide.Sequence %in% names(tt)[tt], ]
  S <- droplevels(S)
  S.training <- reshape(S, direction = 'wide', 
                        v.names = 'Retention.Time', 
                        timevar = c('instrument'), 
                        idvar = 'Peptide.Sequence')
  
  if (plot == TRUE){
    pairs(S.training, 
          pch=as.integer(S.training$Peptide.Sequence), 
          col=as.integer(S.training$Peptide.Sequence),
          lower.panel = NULL)
  }
  return(S.training)
}


.plot_rt_8rep <- function(S, peptides, ...){
  .figure_setup()
  S <- S[grep("[by]", S$Fragment.Ion), ]
  S <- S[S$Peptide.Sequence %in% peptides$Peptide.Sequence, ]
  S <- aggregate(Retention.Time ~ Peptide.Sequence * File.Name.Id * instrument 
                 * relative.amount * Isotope.Label.Type,
                 FUN=mean, 
                 data=S)
  S <- droplevels(S)
  
  xyplot(Retention.Time ~ File.Name.Id  | Isotope.Label.Type * instrument, 
         data = S, 
         layout = c(10, 1),
         group = S$Peptide.Sequence, 
         auto.key = list(space = "right", points = TRUE, lines = FALSE, cex=1),
         ...)
}

.plot_rt_dil <- function(S, peptides,...){
  .figure_setup()
  S <- S[grep("[by]", S$Fragment.Ion), ]
  S <- S[S$Peptide.Sequence %in% peptides$Peptide.Sequence, ]
  S <- aggregate(Retention.Time ~ Peptide.Sequence * File.Name * instrument 
                 * relative.amount * Isotope.Label.Type,
                 FUN=mean, 
                 data=S)
  S <- droplevels(S)
  
  xyplot(Retention.Time ~ relative.amount |  Isotope.Label.Type * instrument , 
         data = S, 
         layout = c(10, 1),
         group = S$Peptide.Sequence, 
         scales = list(x = list(rot = 45, log=TRUE, at=sort(unique(S$relative.amount)) )),
         auto.key = list(space = "right", points = TRUE, lines = FALSE, cex=1),...)
}

