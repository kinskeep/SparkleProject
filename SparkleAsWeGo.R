
setwd('C:/Users/Katie/OneDrive/Documents/Feder Spring 2018/SparkleAnalysis')
reads <- read.csv('Sparklef_indiv_cov.csv')
hist (reads[,3],breaks=20)
plot (reads[,2], reads[,3])

DB <- reads[grep ("DB",reads[,2]),]
SB <- reads[grep ("SB",reads[,2]),]
hist (DB[,3],breaks=20)
hist (SB[,3],breaks=20)
plot (DB[,2], DB[,3])
plot (SB[,2], SB[,3])
nrow(reads[reads[,3]<500000,])

db_sc <- reads[grep ("DB_SC",reads[,2]),]
sb_sc <- reads[grep ("SB_SC",reads[,2]),]
db_al <- reads[grep ("DB_AL",reads[,2]),]
sb_al <- reads[grep ("SB_AL",reads[,2]),]
db_ar <- reads[grep ("DB_AR",reads[,2]),]
sb_ar <- reads[grep ("SB_AR",reads[,2]),]
hist (db_sc[,3],breaks=20)
hist (sb_sc[,3],breaks=20)
hist (db_al[,3],breaks=20)
hist (sb_al[,3],breaks=20)
hist (db_ar[,3],breaks=20)
hist (sb_ar[,3],breaks=20)

depth <- read.table('Sparklef.keep.50.3.30.bi.idepth',header=T,stringsAsFactors=F)
missyou <- read.table('Sparklef.keep.50.3.30.bi.imiss',header=T,stringsAsFactors=F)
hist (depth[,3],breaks=20)
hist (missyou[,5],breaks=20)
print (missyou[missyou[,5]>.6,])
Lowcov <- read.table('LowCovIndiv.75.txt')

HWidepth <- read.table('Sparklef.keep.50.3.30.bi.75.filtered.HWE.idepth',header=T,stringsAsFactors=F)
HWimiss <- read.table('Sparklef.keep.50.3.30.bi.75.filtered.HWE.imiss',header=T,stringsAsFactors=F)
HWldepthmean <- read.table('Sparklef.keep.50.3.30.bi.75.filtered.HWE.ldepth.mean',header=T,stringsAsFactors=F)
hist (HWidepth[,3],breaks=20)
hist (HWimiss[,5],breaks=20)
hist (HWldepthmean[,3],breaks=40,xlim=c(12,25),ylim=c(0,400))
mean(HWldepthmean[,3])+4*sqrt(mean(HWldepthmean[,3]))

DPimiss <- read.table('Sparklef.keep.50.3.30.bi.75.filtered.HWE.DP11.imiss',header=T,stringsAsFactors=F)
print(DPimiss[DPimiss[,5]>.7,])
print(DPimiss[DPimiss[,5]<.1,])
DPidepth <- read.table('Sparklef.keep.50.3.30.bi.75.filtered.HWE.DP11.idepth',header=T,stringsAsFactors=F)
print(DPidepth[DPidepth[,3]<2.5,])
print(DPidepth[DPidepth[,3]>8,])
DPldepthmean <- read.table('Sparklef.keep.50.3.30.bi.75.filtered.HWE.DP11.ldepth.mean',header=T,stringsAsFactors=F)
print(DPldepthmean[DPldepthmean[,3]<1.5,])
print(DPldepthmean[DPldepthmean[,3]>15,])

Finalimiss <- read.table('Sparklef.final.imiss',header=T,stringsAsFactors=F)
Finalidepth <- read.table('Sparklef.final.idepth',header=T,stringsAsFactors=F)
print(Finalimiss[Finalimiss[,5]>.6,])
print(Finalimiss[Finalimiss[,5]<.1,])
print(Finalidepth[Finalidepth[,3]<3,])
print(Finalidepth[Finalidepth[,3]>7,])
hist(Finalimiss[,5])
hist(Finalidepth[,3])

Het <- read.table('Sparklef.keep.50.3.30.bi.75.filtered.HWE.DP14.recode.het',header=T,stringsAsFactors=F)
hist(Het[,5],breaks=20)
Ordered <- Het[order(Het[,5]),]
