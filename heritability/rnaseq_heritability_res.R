library("session")
options(warn=1)
sessionInfo()
options("scipen"=100, "digits"=4)

args <- commandArgs(trailingOnly = TRUE)
phenFp <- args[1]
statusFp <- args[2]
outFp <- args[3]
covStr <- args[4]

# assume name=phen

if (1 == 0){
    phenFp=
    statusFp=
    outFp=
    covStr="PEER_1+PEER_2+PEER_3"
}

t <- read.table(phenFp,header=1,sep="\t",stringsAsFactor=F)
head(t)

id=t[,2]

res <- c()
f <- as.formula(paste0("phen ~ ", covStr))
lm.obj <- lm(f, data=t)
res <- residuals(lm.obj)

t.out <- data.frame(FID=id,IID=id,res=res)
write.table(t.out, file=outFp, sep="\t", row.names=F, col.names=F, quote=F)

cat("END of script\n")
file.create(file=statusFp,overwrite=T)