#'''
#Copyright (c) 2015, David Edwards, Kat Holt
#All rights reserved. (see README.txt for more details)
#'''

getRecombHeatmap<-function(snps,ref,w=NULL,g=NULL,offset=0) {

# remove offset
snps[,1]<-as.numeric(snps[,1])-offset

# get window size
if (is.null(w)){
if (!is.null(g)){
w <- round(g/length(snps)*10,0) # default window size to ensure average of 10 SNPs per window
}
else(return("Can't run. Please specify window size or genome size"))
}

# get genome size

if (is.null(g)) { g <- max(as.numeric(snps[,1])) + w }

# generate density matrix
m<-vector()
i<-0 # strain count

for (strain in colnames(snps)[-1]) {
if (strain!=ref) {
i<-i+1

# snps for strain vs ref
s<-snps[,1][as.character(snps[,which(colnames(snps)==strain)]) != as.character(snps[,which(colnames(snps)==ref)])]

# get counts in windows of size w
h<-hist(s,breaks=seq(1,g,w),plot=F)

if (i==1) { m<-h$counts }
else { m<-rbind(m,h$counts) }

}
}

# return matrix
rownames(m) <- colnames(snps[,-which(colnames(snps)==ref)])[-1]
colnames(m) <- seq(1,g,w)[-1] - 1
m

}


getRecombBetweenStrains<-function(snps,strain1,strain2,w=NULL,g=NULL,offset=0,plotResult=T,multiplier=1) {
s<-snps[,1][as.character(snps[,which(colnames(snps)==strain1)]) != as.character(snps[,which(colnames(snps)==strain2)])]
result<-getRecomb(as.numeric(s)-offset,w=w,g=g,plotResult=plotResult,multiplier=multiplier,strain1,strain2)
return(list(p=result$p,blocks=result$blocks,meanCount=result$meanCount,s=result$s,r=result$r,w=result$w))
}

getRecomb<-function(snps,w=NULL,g=NULL,plotResult=T,multiplier=1,strain1,strain2){
if (is.null(w)){
if (!is.null(g)){
w <- round(g/length(snps)*10,0) # default window size to ensure average of 10 SNPs per window
}
else(return("Can't run. Please specify window size or genome size"))
}
# get counts in windows of size w
if (!is.null(g)){ h<-hist(snps,breaks=seq(1,g+w,w),plot=F) }
else { h<-hist(snps,breaks=seq(1,max(snps)+w,w),plot=F) }
# get p-values for each window
p<-cbind(h$mids,pbinom(h$counts,w,mean(h$counts)*multiplier/w,lower.tail=F),h$counts)   
# Bonferroni-corrected significant p-values
p.sig<-data.frame(p[ p[,2]< 0.05/length(snps), ])
plot(h,ylim=c(0,max(h$counts*1.5)))
if (dim(p.sig)[2]==3) {
blocks<-getBlocks(p.sig,w) 
try(for(i in 1:nrow(blocks)){lines(c(blocks[i,1],blocks[i,2]),rep(max(h$counts)*1.3,2),col=4,lwd=8,pch=15)})
colnames(p.sig)<-c("window_mid","pvalue","numSNPs")
}
else{
blocks<-NULL
if (dim(p.sig)[1]==3) {
lines(c(p.sig[1,1]-w,p.sig[1,1]+w-1),rep(max(h$counts)*1.2,2),col=4,lwd=8,pch=15)
p.sig<-t(matrix(t(p.sig)))
colnames(p.sig)<-c("window_mid","pvalue","numSNPs")
}
}
cleaning<-removeSNPsInBlocks(snps,blocks)
snps_to_remove<-cleaning$remove
clean_snps<-cleaning$keep
if (plotResult){ hist(clean_snps,breaks=seq(1,max(clean_snps)+w,w),add=T,fill=2,col=2,border=2,main=paste(strain1,strain2)) }
return(list(p=p.sig,blocks=blocks,meanCount=mean(h$counts),s=clean_snps,r=snps_to_remove,w=w))
}


getBlocks<-function(p.sig,w) {
window_start<-vector()
window_stop<-vector()
mean_pvalue<-vector()
total_counts<-vector()
# initialise with first window
prev_start<-p.sig[1,1] - w/2
prev_stop<-p.sig[1,1] + w/2 - 1
pvals<-p.sig[1,2]
counts<-p.sig[1,3]
block<-1
for (i in 2:nrow(p.sig)){
start<-p.sig[i,1] - w/2
stop<-p.sig[i,1] + w/2 - 1
p<-p.sig[i,2]
if (start == prev_stop + 1) {
# overlapping blocks
prev_stop <- stop # update to current stop
pvals <- c(pvals,p) # record p-value
counts <- counts + p.sig[i,3]
}
else {
# new block
# record old block info
window_start[block] <- prev_start
window_stop[block] <- prev_stop
mean_pvalue[block] <- mean(pvals)
total_counts[block] <- counts
# initialise new block
block <- block + 1
prev_start <- start
prev_stop <- stop
pvals <- p
counts <- p.sig[i,3]
}
}
# record last block info (!)
window_start[block] <- prev_start
window_stop[block] <- prev_stop
mean_pvalue[block] <- mean(pvals)
total_counts[block] <- counts

result<-cbind(window_start,window_stop,window_stop-window_start+1,total_counts,mean_pvalue)
colnames(result)<-c("start","stop","windowSize","counts","meanPvalue")
return(result)
}


removeSNPsInBlocks<-function(snps,blocks) {
snps_to_remove<-vector()
for (i in 1:nrow(blocks)){
start<-blocks[i,1]
stop<-blocks[i,2]
snps_to_remove <- c(snps_to_remove, snps[snps>start & snps<stop])
}
snps_to_keep <- snps[! (snps %in% snps_to_remove)]
return(list(keep=snps_to_keep,remove=snps_to_remove))
}

getDistAtSnpSet<-function(snps,strain,snpset) {
counts<-vector()
straincol <- which(colnames(snps)==strain)
for (i in 2:ncol(snps)) {
s<-snps[snps[,1] %in% snpset,c(i,straincol)]
snpcount<-c(1:length(s))[s[,1] == s[,2] && s[,2] != "-"]
counts<-c(counts,snpcount)
}
names(counts)<-colnames(snps)[2:ncol(snps)]
return(counts)
}