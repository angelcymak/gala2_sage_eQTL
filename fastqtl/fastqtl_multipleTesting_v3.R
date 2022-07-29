library(qvalue)
library(session)

args <- commandArgs(trailingOnly = TRUE)

fdr<- as.numeric(args[1])
geneFp <- args[2]
outdir <- args[3]
step <- args[4]
varFp_raw <- args[5]
chr <- args[6]


for (i in 1:length(args)){
    cat(paste0(i," <- '",args[i],"'\n"))
}


save.session(file=paste0(outdir,'/out.All.genes.RDa'))

if (file.exists(paste0(outdir,"/log/status.ERROR.R.fastqtl_multipleTesting_v3.",step,".",chr))){
    file.remove(paste0(outdir,"/log/status.ERROR.R.fastqtl_multipleTesting_v3.",step,".",chr))
}


if (step == "fdr"){
    
    sink(paste0(outdir,"/out.All.genes.log"),append=F)
    fastqtl.df <- read.table(gzfile(geneFp),header=T,fill=T,sep="\t")
    print(head(fastqtl.df))
    
    nanrows <- is.na(fastqtl.df[, 'pval_beta'])
    fastqtl.df <- fastqtl.df[!nanrows, ]
    
    cat("  * Number of genes tested: ", nrow(fastqtl.df), " (excluding ",
        sum(nanrows), " genes w/o variants)\n", sep="")
    cat("  * Correlation between Beta-approximated and empirical p-values: ",
        round(cor(fastqtl.df[, 'pval_perm'], fastqtl.df[, 'pval_beta']), 4), "\n", sep="")
    
    Q <- qvalue(fastqtl.df[, 'pval_beta'])
    fastqtl.df$qval <- signif(Q$qvalues, 6)
    cat("  * Proportion of significant phenotypes (1-pi0): " , round((1 - Q$pi0), 2), "\n", sep="")
    cat("  * eGenes @ FDR ", fdr, ":   ", sum(fastqtl.df[, 'qval']< fdr), "\n", sep="")
    
    # determine global min(p) significance threshold and calculate nominal p-value threshold for each gene
    
    ub <- sort(fastqtl.df[fastqtl.df$qval > fdr, 'pval_beta'])[1]  # smallest p-value above FDR
    lb <- -sort(-fastqtl.df[fastqtl.df$qval <= fdr, 'pval_beta'])[1]  # largest p-value below FDR
    
    pthreshold <- (lb+ub)/2
    cat("  * min p-value threshold @ FDR ", fdr, ": ", pthreshold, "\n", sep="")
    
    fastqtl.df[, 'pval_nominal_threshold'] <- signif(qbeta(pthreshold,
        fastqtl.df[, 'beta_shape1'], fastqtl.df[, 'beta_shape2'], ncp=0, lower.tail=TRUE, log.p=FALSE), 6)
    
    write.table(fastqtl.df, gzfile(paste0(outdir,"/out.All.genes.txt.gz")), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
    
    print(summary(fastqtl.df[,"num_var"]))
    print(summary(fastqtl.df[,"tss_distance"]))
    sink()
    
    t.eGene <- fastqtl.df[which(fastqtl.df$qval < fdr),]
    write.table(t.eGene, gzfile(paste0(outdir,"/out.All.eGenes.txt.gz")), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
    
    write(c(),file=paste0(outdir,"/log/z.PASS.fdr"),ncol=1,sep="\t")
    
}else if (step == "eqtl"){
    
    
    
    t.eGene <- read.table(gzfile(paste0(outdir,"/out.All.eGenes.txt.gz")), header=TRUE, sep="\t")
    varFp <- gsub("--CHR--",chr,varFp_raw)
    
    varOutFp <- gsub(".txt",".more.txt",varFp)


    cat("chr=",chr,"\n")
    t.var <- read.table(gzfile(varFp),header=T,fill=T,sep="\t")
    dim(t.var)
    t.merge <- merge(t.var,t.eGene[,c("gene_id","qval","pval_nominal_threshold")], by="gene_id")
    dim(t.merge)
    
    t <- t.merge
    t$chr <- chr
    t$eQTL.flag <- NA
    t$eQTL.flag <- t$pval_nominal < t$pval_nominal_threshold
    write.table(t, gzfile(varOutFp), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
    
}

cat("END of script\n")