rm(list=ls())
pheno <- read.csv('./data/metadata.csv',stringsAsFactors = FALSE)

pheno$time <- gsub('[0-9]','',pheno$Sample.ID)
pheno$time <- gsub(' ','',pheno$time)
pheno$time[pheno$time == 'A'] <- 'day0'
pheno$time[pheno$time == 'B'] <- 'day3'
pheno$time[pheno$time == 'C'] <- 'day7'
pheno$Sample.ID <- gsub(' ','',paste0('S',pheno$Sample.ID))

proteo.data <- read.csv(
  './data/norm-data.csv',
  stringsAsFactors = FALSE
)
colnames(proteo.data) <- gsub('^X','S',colnames(proteo.data))

proteinIDs <- proteo.data$proteinID
proteinID <- sapply(proteinIDs,function(x){strsplit(x,split = '\\|')[[1]][2]})
rownames(proteo.data) <- proteinID

proteo.data['P20851',-1][,59]
# 0
proteo.data[,-1] <- proteo.data[,-1] + 1
pro.log2 <- log2(proteo.data[,-1])
pro.log2.long <- reshape2::melt(pro.log2)

meta.time <- pheno$time
names(meta.time) <- pheno$Sample.ID
meta.pheno <- pheno$Label
names(meta.pheno) <- pheno$Sample.ID
pro.log2.long$time <- meta.time[as.character(pro.log2.long$variable)]
pro.log2.long$phenotype <- meta.pheno[as.character(pro.log2.long$variable)]
meta.batch <- paste0('0',as.numeric(as.factor(pheno$Date.of.run)))
names(meta.batch) <- pheno$Sample.ID
pheno$bacth <- meta.batch
pro.log2.long$batch <- meta.batch[as.character(pro.log2.long$variable)]
library('ggplot2')
library('ggpubr')
p.overview01 <- ggplot(pro.log2.long,aes(value)) +
  geom_density(aes(col=variable)) +
  guides(col=FALSE) +
  xlab('log2(Area Normalized Intensity)') +
  theme_pubr()
p.overview02 <- ggplot(pro.log2.long,aes(value)) +
  geom_density(aes(col=variable)) +
  guides(col=FALSE) +
  xlab('log2(Area Normalized Intensity)') +
  facet_grid(phenotype~time) +
  theme_pubr()
p.overview03 <- ggplot(pro.log2.long,aes(value)) +
  geom_density(aes(linetype=batch)) +
  guides(linetype=FALSE) +
  xlab('log2(Area Normalized Intensity)') +
  facet_grid(phenotype~time) +
  theme_pubr()

library(gridExtra)
dir.create('./figures',showWarnings = FALSE)
pdf('./figures/00-Raw-log2-without-batch-correction.pdf',
    width = 6,height = 8)
grid.arrange(p.overview01,p.overview02,p.overview03,nrow=3)
dev.off()

#------
# BATCH EFFECT 
library(sva)
library(snm)
raw.dat <- pro.log2
bio.var <- data.frame(
  label = as.factor(pheno$Label),
  time = as.factor(pheno$time)
)
adj.var <- data.frame(batch=as.factor(pheno$bacth))

int.var <-  data.frame(array=as.factor(1:nrow(bio.var)))
set.seed(7859)
snm.obj <- snm(raw.dat = as.matrix(raw.dat),
               bio.var = model.matrix(~.,bio.var),
               adj.var = model.matrix(~.,adj.var),
               int.var = int.var
)
norm.dat <- snm.obj$norm.dat
norm.long <- reshape2::melt(norm.dat)
norm.long$time <- meta.time[as.character(norm.long$Var2)]
norm.long$phenotype <- meta.pheno[as.character(norm.long$Var2)]
norm.long$batch <- meta.batch[as.character(norm.long$Var2)]

p.norm01 <- ggplot(norm.long,aes(value)) +
  geom_density(aes(col=Var2)) +
  guides(col=FALSE) +
  xlab('log2(Area Normalized Intensity)') +
  theme_pubr()
p.norm02 <- ggplot(norm.long,aes(value)) +
  geom_density(aes(col=Var2)) +
  guides(col=FALSE) +
  xlab('log2(Area Normalized Intensity)') +
  facet_grid(phenotype~time) +
  theme_pubr()
p.norm03 <- ggplot(norm.long,aes(value)) +
  geom_density(aes(linetype=batch)) +
  guides(linetype=FALSE) +
  xlab('log2(Area Normalized Intensity)') +
  facet_grid(phenotype~time) +
  theme_pubr()
pdf('./figures/00-Raw-log2-with-batch-correction.pdf',
    width = 6,height = 8)
grid.arrange(p.norm01,p.norm02,p.norm03,nrow=3)
dev.off()
norm.proteo <- norm.dat
save(pheno,norm.proteo,file='./data/normalized-data.RData')
rm(list=ls())
load('./data/normalized-data.RData')
pheno.split <- plyr::dlply(pheno,'Label',function(x){
  plyr::dlply(x,'time')
})
#---
# perform differential expression Asym day0, day3, day7
estDE <- function(label,grp1,grp2){
  # label = 'Asymptomatic'
  # grp1 = 'day0'
  # grp2 = 'day3'
  samp.grp1 <- pheno.split[label][[1]][grp1][[1]]$Sample.ID
  samp.grp2 <- pheno.split[label][[1]][grp2][[1]]$Sample.ID
  input <- norm.proteo[,c(samp.grp1,
                          samp.grp2)]
  deout <- apply(input,1,
                 function(x){
                   x1 <- x[1:length(samp.grp1)];
                   y1 <- x[(length(samp.grp2)+1):length(x)]
                   tout <- t.test(x = x1,
                                  y = y1,
                                  alternative = 'two.sided',
                                  paired=TRUE
                   )
                   xx <- data.frame(
                     diff = tout$estimate,
                     pval=tout$p.value
                   )
                 }
  )
  deout.df <- plyr::ldply(deout)
  cols <- ifelse(deout.df$pval<0.05,yes = 'red',no = 'black')
  deout.df$sig <- ifelse(deout.df$pval<0.05,yes = 'T',no = 'F')
  return(list(
    deout = deout.df,
    group=paste0(grp1,'-',grp2)
  )
  )
}
#' Asym, day0, day3, day07
asym.de0.3 <- estDE(label = 'Asymptomatic',grp1 = 'day0',grp2 = 'day3')
asym.de0.7 <- estDE(label = 'Asymptomatic',grp1 = 'day0',grp2 = 'day7')
asym.de3.7 <- estDE(label = 'Asymptomatic',grp1 = 'day3',grp2 = 'day7')

#' Sym, day0, day3, day07
sym.de0.3 <- estDE(label = 'Symptomatic',grp1 = 'day0',grp2 = 'day3')
sym.de0.7 <- estDE(label = 'Symptomatic',grp1 = 'day0',grp2 = 'day7')
sym.de3.7 <- estDE(label = 'Symptomatic',grp1 = 'day3',grp2 = 'day7')

de.all <- rbind(
  # asym
  cbind(asym.de0.3$deout,group=asym.de0.3$group,label='Asymp'),
  cbind(asym.de0.7$deout,group=asym.de0.7$group,label='Asymp'),
  cbind(asym.de3.7$deout,group=asym.de3.7$group,label='Asymp'),
  # sym
  cbind(sym.de0.3$deout,group=sym.de0.3$group,label='Symp'),
  cbind(sym.de0.7$deout,group=sym.de0.7$group,label='Symp'),
  cbind(sym.de3.7$deout,group=sym.de3.7$group,label='Symp')
)

idx.sig <- de.all$sig == 'T'
de.label2days <- de.all[idx.sig,]
save(de.label2days, file='./data/de.label2days.RData')


# perform differential expression Asym day0, day3, day7
pheno.day <- plyr::dlply(pheno,'time',function(x){
  plyr::dlply(x,'Label')
})
estDEday <- function(day){
  # label = 'Asymptomatic
  label = day
  grp1 = 'Asymptomatic'
  grp2 = 'Symptomatic'
  samp.grp1 <- pheno.day[label][[1]][grp1][[1]]$Sample.ID
  samp.grp2 <- pheno.day[label][[1]][grp2][[1]]$Sample.ID
  input <- norm.proteo[,c(samp.grp1,
                          samp.grp2)]
  deout <- apply(input,1,
                 function(x){
                   x1 <- x[1:length(samp.grp1)];
                   y1 <- x[(length(samp.grp2)+1):length(x)]
                   tout <- t.test(x = x1,
                                  y = y1,
                                  alternative = 'two.sided',
                                  paired=TRUE
                   )
                   xx <- data.frame(
                     diff = tout$estimate,
                     pval=tout$p.value
                   )
                 }
  )
  deout.df <- plyr::ldply(deout)
  cols <- ifelse(deout.df$pval<0.05,yes = 'red',no = 'black')
  deout.df$sig <- ifelse(deout.df$pval<0.05,yes = 'T',no = 'F')
  return(list(
    deout = deout.df,
    group=paste0(grp1,'-',grp2)
  )
  )
}
#' Asym, day0, day3, day07
day0.de <- estDEday(day = 'day0')
day3.de <- estDEday(day = 'day3')
day7.de <- estDEday(day = 'day7')
de.day <- rbind(
  cbind(day0.de$deout,group=day0.de$group,label='day0'),
  cbind(day3.de$deout,group=day3.de$group,label='day3'),
  cbind(day7.de$deout,group=day7.de$group,label='day7')
)
idx.sig.day <- de.day$sig == 'T'
table(de.day[idx.sig.day,c('group','label')])
table(de.all[idx.sig,c('group','label')])

de.out <- rbind(
  de.all[idx.sig,],
  de.day[idx.sig.day,]
)
colnames(de.out)[1] <- 'ProteinID'
colnames(de.out)[5:6] <- c('comparison','group')
write.csv(de.out,'./data/DE-proteins-GIMS-COVID19-060321.csv',
          row.names = FALSE)

