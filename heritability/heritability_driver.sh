# submit job in SGE


# e_ : AA PR MX All
# z_ : z.9


eX_=$e_
tssFlag_=.tss
nPeer=50
fastqtlTag=all
switch_=$z_..h2
if [ $switch_ == $switch_ ];
then
    echoFormat "SWITCH:$switch_"
    task_=rnaseq_heritability.scratch.res.$eX_${tssFlag_}
    outdir_=/wynton/group/burchard/maka/project/topmedp4.rnaseq/heritability/${z_}${tssFlag_}/$task_
    
    if [ 1 == 1 ];
    then
        # SNP only
        # Note that if input genotype file is not exact same number of sample in sampeFp, all variants with mac 1 will be kept
        subp_script_=$subp/z_geneExp/rnaseq_heritability.scratch.res.sh
        rFp_geneExpRes_=$subp/z_geneExp/rnaseq_heritability_res.R
        echoTask $task_
        qclogdir_=$outdir_/log
        hold_="-hold_jid a"
        qcFpStr_=$sourceFp
        
        # batch of genes as input
        indexFp_=/wynton/group/burchard/maka/project/topmedp4.rnaseq/predictdb/${z_}${tssFlag_}/all.$eX_/prepare_data/list.1.23.txt
        
        # genotype file will be subset to samples in sampleFp
        genoFp_=$zdir_/bcftools_isec.DP0.subset.mafp01.genop05.$eX_/isec.CHR-----/0000.vcf.gz
        expFp_=/wynton/group/burchard/maka/project/topmedp4.rnaseq/fastqtl/${z_}/fastqtl${tssFlag_}.${e_}.$fastqtlTag/out.-----.expression.bed.gz
        
        sampleFp_=$zdir_/out.$eX_.NWD.txt
        parFp_=$outdir_/par.txt
        # rsync every 50 lines
        save_=50
        # For gene.annotation.parsed.6.txt file
        unconstrainFlag_=1
        chrCol_=1
        geneIdCol_=2
        geneNameCol_=3
        startCol_=4
        stopCol_=5
        
        # separate into ccov (categorical) and qcov (quantitative)
        # to prepare these files: eqtl.1.compare.Rmd > ## z.9 heritability cov prep
        ccovFp_=$zdir_/peer.1000.80.$eX_/out.$eX_.$nPeer.peer.covariates.covc.txt
        qcovFp_=$zdir_/peer.1000.80.$eX_/out.$eX_.$nPeer.peer.covariates.covq.txt
        qcovFp2_=$zdir_/peer.1000.80.$eX_/out.$eX_.$nPeer.covariates.covq.txt
        peerFp_=$qcovFp_
        
        # Generate covariate file : separete into categorical and quantitative (continuous)
        if [ ! -d $outdir_ ]; then mkdir -p $outdir_;fi
        # par1: all covariates
        echo "--covar $ccovFp_ --qcovar $qcovFp_" > $parFp_
        # par2: all covariates without peer for analysis 3 and 4
        echo "--covar $ccovFp_ --qcovar $qcovFp2_" >> $parFp_
        # par3: other GCTA arguments
            
        vStr="qcFpStr=$qcFpStr_,indexFp=$indexFp_,genoFp=$genoFp_,expFp=$expFp_,parFp=$parFp_,sampleFp=$sampleFp_,save=$save_,unconstrainFlag=$unconstrainFlag_,peerFp=$peerFp_,chrCol=$chrCol_,geneIdCol=$geneIdCol_,geneNameCol=$geneNameCol_,startCol=$startCol_,stopCol=$stopCol_,nPeer=$nPeer,rFp_geneExpRes=$rFp_geneExpRes_"
        
        MEM_=30G
        HR_=336:00:00
        
        if [ ! -f $indexFp_ ];
        then
            echo Missing:$indexFp_
            kill -SIGINT $$
        fi
        
        n_=`cat $indexFp_ | wc -l`
        
        t_=1-$n_
        nT_=$n_
        
        
        
        qcStr_="END"
        
        echoFormat "Number of job to be exluded: $(wc -l $excludeSGEFp_) \n"
        vStr=$vStr",excludeSGEFp=$excludeSGEFp_"
        
        qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$qclogdir_,qctag="
        runTask $task_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" "$vStr"
        
    fi
    
fi ## Heritability using GCTA : z9 tss (unconstrained)
