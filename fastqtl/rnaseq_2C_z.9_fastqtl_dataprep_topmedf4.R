if ( 1 == 0){
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
    BiocManager::install("edgeR")
    
    # "R CMD INSTALL /wynton/home/burchard/maka/data_shared/workflow.sh/R/R_peer_source_1.3.tgz"
}


library("session")
library(peer)
library(fastDummies)

options(warn=1)
sessionInfo()
options("scipen"=100, "digits"=4)

args <- commandArgs(trailingOnly = TRUE)
outFpI <- args[1]
inFpI <- args[2]
k <- as.numeric(args[3])
# max_iter <- 1000
max_iter <- as.numeric(args[4])
nPeer <- as.numeric(args[5])
invnormFp <- args[6]
statusFp <- args[7]

for (i in 1:length(args)){
    cat(paste0(i," <- '",args[i],"'\n"))
}

if (1 == 0){
    basedir <- "/wynton/group/burchard/maka"
    homedir <- "/wynton/group/burchard/maka"
    
    source(paste0(homedir,"/data_shared/workflow.sh/bin.sub/general.qb3.wf.R"))
    inFpI <- "/wynton/group/burchard/maka/project/topmedp4.rnaseq/gtex.preprocessing/run3B.ERCC/run3B.ERCC"
    outFpI <- "/wynton/group/burchard/maka/project/topmedp4.rnaseq/gtex.preprocessing/run3B.ERCC/z.6/peer.100.5/out"
    k <- 4
    nPeer <- 5
    max_iter <- 100
    geneFp <- paste0("/wynton/group/burchard/maka/data_shared/GRCh38/gencode.v30.GRCh38.annotation.ERCC.gOnly.bed")
    invnormFp="/wynton/group/burchard/maka/project/topmedp4.rnaseq/gtex.preprocessing/run3B.ERCC/z.2C/out.All.count_matrix.invnorm.rds"
    
}

geneFp <- paste0("/wynton/group/burchard/maka/data_shared/GRCh38/gencode.v30.GRCh38.annotation.ERCC.gOnly.bed")
v.elong <- c("African American","Mexican","Puerto Rican","All","3pop")
v.e <- c("AA","MX","PR","All","3pop")

e <- v.e[k]
cat(e,"\n")

t.phen.e.r1 <- readRDS(file=paste0(inFpI,".",e,".phen.v1.rds"))
t.phen.e.r1$e <- factor(t.phen.e.r1$ethnicity.am,levels=c("African American","Puerto Rican","Mexican","Other Latino"),labels=c("AA","PR","MX","LA"))
t.phen.e.r1$e <- as.character(t.phen.e.r1$e)


# convert categorical to dummy (e)
cov.pop <- c()
if (length(unique(t.phen.e.r1$e)) > 1){
    t <- dummy_cols(t.phen.e.r1$e,remove_first_dummy = TRUE)
    colnames(t) <- paste0("e",colnames(t))
    cov.pop <- colnames(t)[-1]
    t.phen.e <- cbind(t.phen.e.r1,t[,-1,drop=F])
}else{
    t.phen.e <- t.phen.e.r1
}


# Major difference from script before: no batch variable

cov.base <- c("age","Sex","asthma","PC1","PC2","PC3","PC4","PC5")

cov <- c(cov.base,cov.pop)

print(cov)


t.cov <- t.phen.e[,cov]

for (var in c("Sex","asthma",cov.pop)){
    t.cov[,var] <- as.factor(t.cov[,var])
}

for (var in c("age","PC1","PC2","PC3","PC4","PC5")){
    t.cov[,var] <- as.vector(scale(t.cov[,var],center=T))
}

sample_id <- t.phen.e$NWD_ID
torid <- t.phen.e$TORID
rownames(t.cov) <- sample_id


invnorm.r1 <- readRDS(file=invnormFp)
invnorm <- invnorm.r1[,as.character(torid)]
colnames(invnorm) <- sample_id


# plot
if (1 == 0){
    identical(head(colnames(invnorm)),head(rownames(t.cov)))
    res_tmp <- t(invnorm)
    res_tmp2 <- cbind(res_tmp,t.cov)
    res <- res_tmp
    nGene <- ncol(invnorm)
    
    f <- as.formula(paste0("phen ~ ", paste0(colnames(t.cov),collapse=" + ")))
    
    for (i in 1:nGene){
        cat(i,"\n")
        res_tmp2$phen <- res_tmp2[,i]
        
        lm.obj <- lm(f, data=res_tmp2)
        res[,i] <- residuals(lm.obj)
    }
    
    for (i in 1:nGene){
        hist(invnorm[1,])
        # plot(res[1,])
    }
}


#  -------------------------------
## Part 2 : calculate PEER factors
#  -------------------------------

nSample <- ncol(invnorm)

# need invnorm from Part 1
alphaprior_a <- 0.001
alphaprior_b <- 0.01
epsprior_a <- 0.1
epsprior_b <- 10

if ( 1 == 0 ){
    if (nSample >= 350){
        nPeer=60
    }else if (nSample >= 250 & nSample < 350){
        nPeer=45
    }else if (nSample >= 150 & nSample < 250){
        nPeer=30
    }else if (nSample < 150){
        nPeer=15
    }
    cat("nPeer=",nPeer,"\n")
    nGene <- nrow(invnorm)
    write(c(nSample,nPeer,nGene,e),file=paste0(outFpI,".log.txt"),ncol=4,append=T)
}


expr = t(invnorm)
if (1 == 1){
# Peer
    # row = samples; col = gene
    
    dim(expr)
    
    model = PEER()
    
    # if returns NULL means no error
    PEER_setPhenoMean(model,as.matrix(expr))
    dim(PEER_getPhenoMean(model))
    
    # to infer K=10 hidden confounders
    PEER_setNk(model,nPeer)
    PEER_getNk(model)
    
    # Finally, the prior parameters on the noise and weight precision distributions can also be changed. 
    # As these are both gamma distributed, you can specify the a and b parameters of both:
    PEER_setPriorAlpha(model, alphaprior_a, alphaprior_b)
    PEER_setPriorEps(model, epsprior_a, epsprior_b)
    
    # As default, PEER iterates through updates of every variable 1000 times. To set it to say, 100, use
    PEER_setNmax_iterations(model, max_iter)
    
    # If there are measured experimental variables that may contribute to variability in the data, they can be included in the inference.
    # row=sample, col=covariates 
    head(as.matrix(t.cov))
    head(data.matrix(t.cov))
    save.session(file=paste0(outFpI,".",e,".",nPeer,".peer.tmp.RDa"))
    # cov
    PEER_setCovariates(model, data.matrix(t.cov))
    
    # Perform the inference: The result is the model object with posterior distributions of the variables.
    #   about 30 min - 90 min (why more time needed for PR?) 
    cat("run model: PEER_update\n")
    PEER_update(model)
    cat("PEER_update Done\n")


    # posterior mean of the inferred confounders (NxK matrix)
    # model is not restorable by saving as session
    factors = PEER_getX(model)
    weights = PEER_getW(model)
    alpha = PEER_getAlpha(model)
    R = PEER_getResiduals(model)
    
    dim(factors)
    dim(weights)
    dim(alpha)
    dim(R)
    
    head(factors,5)
    
    saveRDS(factors,file=paste0(outFpI,".",e,".",nPeer,".peer.factor.rds"))
    saveRDS(weights,file=paste0(outFpI,".",e,".",nPeer,".peer.weight.rds"))
    saveRDS(alpha,file=paste0(outFpI,".",e,".",nPeer,".peer.alpha.rds"))
    saveRDS(R,file=paste0(outFpI,".",e,".",nPeer,".peer.residual.rds"))
}

if (1 == 0){
    factors <- readRDS(file=paste0(outFpI,".",e,".",nPeer,".peer.factor.rds"))
    weights <- readRDS(file=paste0(outFpI,".",e,".",nPeer,".peer.weight.rds"))
    alpha <- readRDS(file=paste0(outFpI,".",e,".",nPeer,".peer.alpha.rds"))
    R  <- readRDS(file=paste0(outFpI,".",e,".",nPeer,".peer.residual.rds"))
}

factors <- data.frame(factors)
colnames(factors) <- c(colnames(t.cov),paste0("PEER_",1:nPeer))
rownames(factors) <- rownames(expr)
dim(factors)
nZeroN <- sum(!is.na(factors$PEER_1))

if (nZeroN == 0){
    failFp <- gsub("END","FAIL",statusFp)
    file.create(file=failFp,overwrite=T)
}else{
    # precision (inverse variance) of the weights
    colnames(weights) <- c(colnames(t.cov),paste0("PEER_",1:nPeer))
    rownames(weights) <- colnames(expr)
    dim(weights)
    
    
    bmp(file=paste0(outFpI,".",e,".",nPeer,".peer.alpha.bmp"))
    plot(alpha)
    dev.off()
    
    colnames(R) <- colnames(expr)
    rownames(R) <- rownames(expr)
    
    # For fastqtl: row covariates, col samples
    # NWD
    sampleHeader <- paste0("id\t",paste(sample_id,sep="",collapse="\t"))
    write(sampleHeader, file=paste0(outFpI,".",e,".",nPeer,".peer.covariates.txt"),append=F)
    write.table(t(factors), file=paste0(outFpI,".",e,".",nPeer,".peer.covariates.txt"), col.names=F, row.names=T, sep="\t",quote=F,append=T)
    system(paste0("module load CBI && module load htslib && bgzip -f ",outFpI,".",e,".",nPeer,".peer.covariates.txt"))
    
    
    #  ---------------------------------
    ## Part 4 : Generate covariates file
    #  ---------------------------------
    
    
    # merge phenotype with PEER factors: col var
    t.cov.str2 <- data.frame(NWD_ID=sample_id, factors)
    
    # for general
    write.table(t.cov.str2, file=paste0(outFpI,".",e,".",nPeer,".peer.covariates.t.txt"), col.names=T, row.names=F, sep="\t",quote=F)
    system(paste0("module load CBI && module load htslib && bgzip -f ",outFpI,".",e,".",nPeer,".peer.covariates.t.txt"))
    
    cat("END of script\n")
    file.create(file=statusFp,overwrite=T)
}

