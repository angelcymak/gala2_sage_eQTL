library(edgeR)
library("session")
options(warn=1)
sessionInfo()
options("scipen"=100, "digits"=4)

args <- commandArgs(trailingOnly = TRUE)
# as.numeric()
k <- as.numeric(args[1])

# in case we need differernt dir
basedir <- "/wynton/group/burchard/maka"
homedir <- "/wynton/group/burchard/maka"

source(paste0(homedir,"/data_shared/workflow.sh/bin.sub/general.qb3.wf.R"))

outdir <- paste0(basedir,"/project/topmedp4.rnaseq/gtex.preprocessing/z.9")
outFI <- "out"
outFpI <- paste0(outdir,"/",outFI)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE, mode = "0755")

v.elong <- c("African American","Mexican","Puerto Rican","All","3pop")
v.e <- c("AA","MX","PR","All","3pop")

t.phen <- read.table(file=paste0(outdir,"/t.phen.wgs.nodup.qc.relN.2733.txt"),header=T,sep="\t",stringsAsFactor=F)

rownames(t.phen) <- t.phen$TORID
table(t.phen$ethnicity.am,t.phen$asthma)


# TOPMed files
    countMatFp <- paste0(homedir,"/project/topmedp4.rnaseq/data/TOPMed_Burchard_P4.RNASeQCv2.3.4_gene_reads.gct.gz")
    tpmFp <- paste0(homedir,"/project/topmedp4.rnaseq/data/TOPMed_Burchard_P4.RNASeQCv2.3.4_gene_tpm.gct.gz")
    geneFp <- paste0(homedir,"/data_shared/GRCh38/gencode.v30.GRCh38.annotation.ERCC.gOnly.bed")

## read data
    counts.df <- read.table(gzfile(countMatFp),header=T,skip=2,sep="\t")
    head(counts.df)[,1:5]
    rownames(counts.df) <- counts.df$Name
    counts <- counts.df[,t.phen$TORID]
        
    tpm.df <- read.table(gzfile(tpmFp),header=T,skip=2,sep="\t")
    rownames(tpm.df) <- tpm.df$Name
    
    tpm <- tpm.df[,t.phen$TORID]


#  --------------------------------------------------------------------
## Last step : remove population outliers by kmenas clustering of PC1-3
#  --------------------------------------------------------------------

# see r01.burchard_endotype.log.txt > # 06/11/2021 01:26:33 AM: rerun QC for eQTL (v2021_0611)

#  ----------------------
## PART 1: gene filtering
#  ----------------------

# Based on GTEx v8 : https://gtexportal.org/home/documentationPage

# this script is part of the older version : rnaseq_2_fastqtl_dataprep_topmedf4.R

# counts.e | *.All.count.v1.rds : based on countMatFp after removing pop outliers
#                        | *.{PR.MX.AA}.count.v1.rds : subset by population from All
# tpm.e | *.All.tpm.v1.rds : based on tpmFp after removing pop outliers
#                  | *.{PR,MX,AA}.tmp.v1.rds : subset by population from All


    e <- v.e[k]
    elong <- v.elong[k]
    
    cat(e,"\n")
    
    counts.e <- NULL
    tpm.e <- NULL
    if (k <= 3){
        flag.e <- t.phen$ethnicity.am == elong
        table(flag.e,useNA="always")
        
        t.phen.e <- t.phen[flag.e,]
        counts.e <- counts[,flag.e]
        tpm.e <- tpm[,flag.e]
    }else if (k == 4){ #All
        t.phen.e <- t.phen
        counts.e <- counts
        tpm.e <- tpm
        
        # NWD TOR key
        write.table(t.phen.e[,c("NWD_ID","TORID")], file=paste0(outdir,"/",outFI,".",e,".key.txt"), row.names=F, col.names=T,sep="\t",quote=F)
    }else if (k == 5){
        flag.e <- t.phen$ethnicity.am == "African American" | t.phen$ethnicity.am == "Puerto Rican" | t.phen$ethnicity.am == "Mexican"
        table(flag.e,useNA="always")
        
        t.phen.e <- t.phen[flag.e,]
        counts.e <- counts[,flag.e]
        tpm.e <- tpm[,flag.e]
    }
    
    
    cat("Saving phenotype\n")
    saveRDS(t.phen.e,file=paste0(outdir,"/",outFI,".",e,".phen.v1.rds"))
    cat("Saving count\n")
    saveRDS(counts.e,file=paste0(outdir,"/",outFI,".",e,".count.v1.rds"))
    cat("Saving TPM\n")
    saveRDS(tpm.e,file=paste0(outdir,"/",outFI,".",e,".tpm.v1.rds"))
    
    nrow(t.phen.e)
    
    write.table(t.phen.e$TORID, file=paste0(outdir,"/",outFI,".",e,".TOR.txt"),row.names=F, col.names=F,quote=F, sep="\t")
    write.table(t.phen.e$NWD_ID, file=paste0(outdir,"/",outFI,".",e,".NWD.txt"),row.names=F, col.names=F,quote=F, sep="\t")
    
    # gene filtering
        nGene0 <- nrow(counts.e)
        print(dim(tpm.e))
        
        # filter by per population tpm
        nSample <- ncol(tpm.e)
        count.tpm <- apply(tpm.e,1,function(x)sum(x > 0.1))
        flag.tpm <- count.tpm/nSample >= 0.2 
        nGene1 <- sum(flag.tpm)
        
        # filter by per population count
        count.raw <- apply(counts.e,1,function(x)sum(x>=6))
        flag.raw <- count.raw/nSample >= 0.2
        nGene1B <- sum(flag.raw)
        
        # remove ERCC spike-in
        i.ercc <- grep("ERCC",rownames(counts.e))
        flag.ercc <- !(1:nGene0 %in% i.ercc)
        
        # filter by both
        flag.both <- flag.tpm & flag.raw & flag.ercc
        nGene2 <- sum(flag.both)
        
        counts.e.qc <- counts.e[flag.both,]
        nGene2B <- nrow(counts.e.qc)
        
        out <- c(e,nSample,nGene0,nGene1,nGene1B,nGene2,nGene2B)
        names(out) <- c("e","nSample","nGeneInput","pass.0.1.tpm.0.2.sample","pass.6.reads.0.2.sample","pass.both","final")
        if (k == 1){
            write(names(out),file=paste0(outdir,"/",outFI,".nGene.txt"),ncol=length(out),append=T,sep="\t")
        }
        write(out,file=paste0(outdir,"/",outFI,".nGene.txt"),ncol=length(out),append=T,sep="\t")
        saveRDS(counts.e.qc,file=paste0(outdir,"/",outFI,".",e,".count.geneOutliersN.cor.rds"))
        
        tpm.e.qc <- tpm.e[flag.both,]
        saveRDS(tpm.e.qc,file=paste0(outdir,"/",outFI,".",e,".tpm.geneOutliersN.cor.rds"))
        
        
    
    
    # -----------------------------------
    # Prepare the inverse normalized data
    # -----------------------------------
        cat("Inverse normalization\n")
        # create edgeR object with filtered count table in it
        # rows=genes
    
        # counts.e.qc <- readRDS(file=paste0(outdir,"/",outFI,".",e,".count.geneOutliersN.cor.rds"))
        dge <- DGEList(counts.e.qc)
        
        # Run TMM normalization with edgeR
        # If object is a DGEList, then it is returned as output with the relative normalization factors in object$samples$norm.factors
        dge2 <- calcNormFactors(dge)
        
        # Compute counts per million (CPM) or reads per kilobase per million (RPKM)
        y_cpm1 <- cpm(dge2)
        
        y_cpm <- y_cpm1
        invnorm <- y_cpm1
        
        ###Inverse normalization###
        for(i in 1:length(y_cpm[,1])){
            # https://stackoverflow.com/questions/25459864/making-non-normal-data-normal
            invnorm[i,] <- qnorm(rank(y_cpm[i,], ties.method = 'average')/(length(y_cpm[i,])+1))
        }
        saveRDS(invnorm,file=paste0(outFpI,".",e,".count_matrix.invnorm.rds"))
        saveRDS(y_cpm,file=paste0(outFpI,".",e,".count_matrix.cpm.rds"))
        
        # For data sharing
        outFp <- paste0(outFpI,".",e,".count_matrix.invnorm.txt")
        header=paste0("Gene\t",paste(colnames(invnorm),sep="",collapse="\t"))
        write(header,file=outFp,append=F)
        write.table(invnorm,file=outFp,col.names=F,row.names=T,sep="\t",quote=F,append=F)
        system(paste0("gzip ",outFp))
        
        outFp <- NULL
        outFp <- paste0(outFpI,".",e,".count_matrix.cpm.txt")
        header=paste0("Gene\t",paste(colnames(y_cpm),sep="",collapse="\t"))
        write(header,file=outFp,append=F)
        write.table(y_cpm,file=outFp,col.names=F,row.names=T,sep="\t",quote=F,append=F)
        system(paste0("gzip ",outFp))
        
    
    ## prepare expression file (gz) for fastqtl
        cat("Prepare fastqtl expression data\n")
        t.gene <- read.table(geneFp,header=1,sep="\t")
        flag.gene <- t.gene$geneId %in% rownames(counts.e.qc)
        table(flag.gene, useNA="always")
        
        # overlap of gene in count matrix with gene file, make sure that the number of genes matched
        flag.tmp <- rownames(counts.e.qc) %in% t.gene$geneId 
        table(flag.tmp, useNA="always")
        
        t.gene.filter <- t.gene[flag.gene,]
        order.gene <- match(rownames(counts.e.qc),t.gene.filter$geneId)
        t.gene.filter.order <- t.gene.filter[order.gene,]
        dim(counts.e.qc)
        dim(t.gene.filter.order)
        head(t.gene.filter.order)
        
        # merge gene coordinate with inverse normalized expression and sort genes by coordinates
        invnorm.fastqtl <- cbind(t.gene.filter.order[,c("chrf2","start","end","geneId")],invnorm)
        t.sort <- t.gene.filter.order[,c("chrf2","start","end","geneId")]
        t.sort$Chr <- gsub("chr","",t.sort$chrf2)
        order.sort <- order(t.sort$Chr,t.sort$start,t.sort$end)
        length(order.sort)
        invnorm.fastqtl.sort <- invnorm.fastqtl[order.sort,]
        
        fastqtlFp <- paste0(outFpI,".",e,".expression.bed")
        
        header=paste0("#",paste(c("Chr","Start","End","Geneid"),sep="",collapse="\t"),"\t",paste(t.phen.e$NWD_ID,sep="",collapse="\t"))
        write(header,file=fastqtlFp,append=F)
        write.table(invnorm.fastqtl.sort,file=fastqtlFp,col.names=F,row.names=F,sep="\t",quote=F,append=T)
        
        command.zip <- paste0("module load CBI && module load htslib && bgzip -f ",fastqtlFp," && tabix -f -p bed ",fastqtlFp,".gz")
        system(command.zip)
        

#  ---------
## Next step
#  ---------

    # Estimate peers, format into fastqtl input format
    
    # script : runaseq_2C_z.9_fastqtl_dataprep_topmedf4.R