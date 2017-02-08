
library(data.table)
AFR <- fread("zcat AFR.frq.gz", header=TRUE)
EUR <- fread("zcat EUR.frq.gz", header=TRUE)
EAS <- fread("zcat EAS.frq.gz", header=TRUE)
SAS <- fread("zcat SAS.frq.gz", header=TRUE)

dups <- unique(AFR$SNP[duplicated(AFR$SNP)])
AFR <- subset(AFR, ! SNP %in% dups)
EUR <- subset(EUR, ! SNP %in% dups)
EAS <- subset(EAS, ! SNP %in% dups)
SAS <- subset(SAS, ! SNP %in% dups)

aa <- fread("zcat ../aa_clean.txt.gz", header=FALSE)
dups <- unique(aa$V1[duplicated(aa$V1)])
aa <- subset(aa, ! V1 %in% dups)
aa <- subset(aa, V1 %in% AFR$SNP)
AFR <- subset(AFR, SNP %in% aa$V1)
EUR <- subset(EUR, SNP %in% aa$V1)
EAS <- subset(EAS, SNP %in% aa$V1)
SAS <- subset(SAS, SNP %in% aa$V1)

all(AFR$SNP == aa$V1)
all(EUR$SNP == aa$V1)
all(EAS$SNP == aa$V1)
all(SAS$SNP == aa$V1)

table(AFR$A1 == EUR$A1 | AFR$A1 == EUR$A2)
table(SAS$A1 == EUR$A1 | SAS$A1 == EUR$A2)
table(SAS$A1 == EAS$A1 | SAS$A1 == EAS$A2)
table(AFR$A2 == aa$V2 | AFR$A1 == aa$V2)

index <- (AFR$A2 == aa$V2 | AFR$A1 == aa$V2)

aa <- aa[index, ]
AFR <- AFR[index, ]
EUR <- EUR[index, ]
EAS <- EAS[index, ]
SAS <- SAS[index, ]

index <- AFR$A2 == aa$V2
AFR$AF <- 1 - AFR$MAF
AFR$AF[!index] <- 1 - AFR$AF[!index]

index <- EUR$A2 == aa$V2
EUR$AF <- 1 - EUR$MAF
EUR$AF[!index] <- 1 - EUR$AF[!index]

index <- EAS$A2 == aa$V2
EAS$AF <- 1 - EAS$MAF
EAS$AF[!index] <- 1 - EAS$AF[!index]

index <- SAS$A2 == aa$V2
SAS$AF <- 1 - SAS$MAF
SAS$AF[!index] <- 1 - SAS$AF[!index]

save(AFR, EUR, EAS, SAS, aa, file="freq_aa.rdata")


##


load("freq_aa.rdata")


calc_fst <- function(dat)
{
	m1 <- rowMeans(dat)
	m2 <- 1 - m1

	het <- apply(dat, 2, function(x) {
		2 * x * (1-x)
	})
	hs <- rowMeans(het)
	ht <- 2 * m1 * m2
	fst <- (ht - hs) / ht
	return(fst)
}

dat <- data.frame(AFR = AFR$AF, EUR = EUR$AF, EAS = EAS$AF, SAS = SAS$AF)

fst_all <- calc_fst(dat)
hist(fst_all, breaks=100)

dat <- data.frame(SNP=AFR$SNP, fst=fst_all, AFR=dat$AFR, EUR=dat$EUR, EAS=dat$EAS, SAS=dat$SAS)
save(dat, file="all_fst.rdata")
