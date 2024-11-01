#library(missMDA)
library(maps)

loci <- c( "05", "08", "11", "23", "24", "26", "30")

data<-read.table("raw_genotypes.tsv", as.is=T, header=T)
data[data==0] <- NA
    

## include <- c(sample(sum(data$Country=="CAMEROON"), 25), which(data$Country!="CAMEROON"))
## data <- data[include,]

recode <- data[,1:2]

for(l in loci){
    locus1 <- data[,paste0("Locus_", l, "_1")]
    locus2 <- data[,paste0("Locus_", l, "_2")]

    alleles<-unique(c(locus1, locus2))
    alleles <- alleles[!is.na(alleles)]

    if(length(alleles)<2){ #if monomorphic, ignore
        cat(paste0("Ignoring monomorphic locus " ,l))
        next
    }

    counts <- sapply(alleles, function(x,y){sum(x==y, na.rm=TRUE)}, y=c(locus1, locus2))
    counts <- data.frame(allele=alleles, count=counts)
    counts <- counts[order(-counts$count),]

    major <- counts[1,1]
    for(i in 2:(NROW(counts))){
        this.allele.name <- paste0("Locus_", l, "_", counts[i,1])
        genotype <- (locus1==counts[i,1])+(locus2==counts[i,1])
        recode[,this.allele.name] <- genotype
    }
}


gt <-recode[,3:NCOL(recode)]

pops <- sort(unique(data$Country))
#cols <- brewer.pal(8, "Set1")
cols <- scan("cols.txt", "")
names(cols) <- pops

pop.complete <- data$Country[!apply(is.na(gt), 1, any)]
subpop.complete <- data$Population[!apply(is.na(gt), 1, any)]

gt.complete <- gt[!apply(is.na(gt), 1, any),]
pca.complete <- prcomp(t(gt.complete))$rotation

pdf("PCA_complete_obs.pdf", width=4, height=4)
                                        #plot(jitter(-pca.complete[,1], 100), jitter(-pca.complete[,2], 100), col=cols[pop.complete], pch=16, xlab="PC1", ylab="PC2", main="Complete observations")
par(mar=c(2,2,0.25,0.25))
plot(jitter(pca.complete[,1], 100), jitter(pca.complete[,2], 100), col=cols[pop.complete], pch=16, xlab="PC1", ylab="PC2", xaxt="n", yaxt="n")
                                        #legend("topright", names(cols), pch=16, col=cols, bty="n", cex=0.8)
mtext(side=1, "PC1", line=0.5)
mtext(side=2, "PC2", line=0.5)

dev.off()

## gt.imputed <- imputePCA(gt, mdtaethod="EM", ncp=2)$completeObs
## pca.imputed <- prcomp(t(gt.imputed))$rotation

## pdf("PCA_imputed_obs.pdf")
## plot(pca.imputed[,1], pca.imputed[,2], col=cols[data$Country], pch=16, xlab="PC1", ylab="PC2", main="Imputed observations")
## legend("bottomright", names(cols), pch=16, col=cols, bty="n", cex=0.8)
## dev.off()

pdf("Figure1.pdf", width=8)
par(mar=c(5.1, 6.1, 4.1, 2.1))
image(as.matrix(cbind(NA, gt)), col=c("lightblue", "#2D6BB4", "black"), xaxt="n", yaxt="n", bty="n")
mtext(c("", colnames(gt)), 2, at=(0:(NCOL(gt)))/(NCOL(gt)), las=2, line=0.5)
for(pop in pops){
    start=min(which(data$Country==pop)-1)/NROW(gt)
    end=max(which(data$Country==pop))/NROW(gt)
    rect(start, -0.02, end, 0.02, col=cols[pop], border=cols[pop])
    mtext(pop, 1, at=(start+end)/2, line=0, cex=0.75)
}
dev.off()

N <- NROW(recode)
pedfile <- data.frame(FID=data$Country, IID=data$Country, FID=rep(0,N), MID=rep(0,N), Sex=rep(0,N), Phe=rep(0,N))
gt <- recode[,3:NCOL(recode)]
pedfile <- cbind(pedfile, ifelse(is.na(gt), "0 0", ifelse(gt==0, "A A", ifelse(gt==1, "A T", ifelse(gt==2, "T T", NA)))))

M <- NCOL(gt)
map <- data.frame(chr=rep(1,M), colnames(gt), rep(0, M), 1:M)

write.table(pedfile, "faidherbia.ped", row.names=F, col.names=F, quote=F, sep= " " )
write.table(map, "faidherbia.map", row.names=F, col.names=F, quote=F, sep= " " )

system("plink --file faidherbia --make-bed --out faidherbia")
system("admixture --cv faidherbia.bed 2")

pdf("admixture_k2.pdf", width=8, height=2)
par(mar=c(3,1,1,1))

qmat <- read.table("faidherbia.2.Q", as.is=TRUE)
barplot(t(qmat), col=c("red", "blue"), yaxt="n", space=0, border="lightgrey")
par(xpd=TRUE)
for(pop in pops){
    start=min(which(data$Country==pop)-1)
    end=max(which(data$Country==pop))
    rect(start, -0.2, end, -0.1, col=cols[pop], border=cols[pop])
    mtext(substr(pop, 1, 3), 1, at=(start+end)/2, line=1, cex=0.8)
}
dev.off()

system("admixture --cv faidherbia.bed 3")

pdf("admixture_k3.pdf", width=8, height=2)
par(mar=c(3,1,1,1))
qmat <- read.table("faidherbia.3.Q", as.is=TRUE)
barplot(t(qmat), col=c("red", "blue", "green"), yaxt="n", space=0, border="lightgrey")
par(xpd=TRUE)
for(pop in pops){
q    start=min(which(data$Country==pop)-1)
    end=max(which(data$Country==pop))
    rect(start, -0.2, end, -0.1, col=cols[pop], border=cols[pop])
    mtext(substr(pop, 1, 3), 1, at=(start+end)/2, line=1, cex=0.8)
}
dev.off()

system("admixture --cv faidherbia.bed 4")

pdf("admixture_k4.pdf", width=8, height=2)
par(mar=c(3,1,1,1))
qmat <- read.table("faidherbia.4.Q", as.is=TRUE)
barplot(t(qmat), col=c("red", "blue", "green", "grey"), yaxt="n", space=0, border="lightgrey")
par(xpd=TRUE)
for(pop in pops){
    start=min(which(data$Country==pop)-1)
    end=max(which(data$Country==pop))
    rect(start, -0.2, end, -0.1, col=cols[pop], border=cols[pop])
    mtext(substr(pop, 1, 3), 1, at=(start+end)/2, line=1, cex=0.8)
}
dev.off()

locations <- read.table("locations.txt", as.is=T, header=TRUE)
yrng=c(-35,40)
xrng=c(-20, 50)
pdf("Map.pdf", width=4, height=4)
map('world',col="grey", fill=TRUE, bg="white", lwd=0.05, mar=rep(0.25,4),border=0, ylim=yrng, xlim=xrng )
points(locations$Lon, locations$Lat, pch=16, col=cols[locations$Country], cex=1)
legend("bottomleft", names(cols), pch=16, col=cols, bty="n", cex=1)

dev.off()
