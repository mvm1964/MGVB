# PEP calculations for gMGVB

# Import KDE package
library(ks)

# KDE computing function, uses ks, copy/paste to R prompt to load
computePEP<- function(data, x, nmod) {
	data1<- data[data$nmod1 == nmod,]
	pep_f<- ks::kde(x = cbind(data1$score_1[data1$decoy1==0], data1$plen1[data1$decoy1==0]), eval.points = x, binned=F)$estimate
	pep_r<- ks::kde(x = cbind(data1$score_1[data1$decoy1==1], data1$plen1[data1$decoy1==1]), eval.points = x, binned=F)$estimate
	return(pep_r/(pep_r + pep_f))
}
	

# Import data
data<- read.csv("results.txt", sep="\t", stringsAsFactors=F)	

# Will compute separately for 0, 1, 2, and 3 nmod

# First for 0 nmod
df_nmod0<- data[data$nmod1==0, ]
x<- df_nmod0[, c(3,9)]
df_nmod0$PEP<- computePEP(data, x, 0)

# For 1 nmod
df_nmod1<- data[data$nmod1==1, ]
x<- df_nmod1[, c(3,9)]
df_nmod1$PEP<- computePEP(data, x, 1)

# For 2 nmod
df_nmod2<- data[data$nmod1==2, ]
x<- df_nmod2[, c(3,9)]
df_nmod2$PEP<- computePEP(data, x, 2)

# For 3 nmod
df_nmod3<- data[data$nmod1==3, ]
x<- df_nmod3[, c(3,9)]
df_nmod3$PEP<- computePEP(data, x, 3)

# Merge them
df_final<- rbind(df_nmod0, df_nmod1, df_nmod2, df_nmod3)

# Now use PEP to filter at 1% FDR
df_final_ordered<- df_final[order(df_final$PEP),]
df_final_ordered$cum_decoy<- cumsum(df_final_ordered$decoy1)
df_final_ordered$rownum<- 1:nrow(df_final_ordered)
df_final_ordered$FDR<- df_final_ordered$cum_decoy/df_final_ordered$rownum
df_sig<- df_final_ordered[df_final_ordered$FDR<=0.01,]
write.table(df_sig[,1:9], "sig_results.txt", sep = "\t", row.names=F)