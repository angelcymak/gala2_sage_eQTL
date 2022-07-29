# ------------- predictdb z.9 (tss or gene)

# prepare expression file (gene version for predictdb)
if [ 1 == 0 ];
then
    # all | case | ctrl
    fastqtlTag=all
    tssFlag_=
    
    for nPeer in 50;
    do
        eX_=$e_
        
        runInfo_="fastqtl z=$z_ eX=$eX_ tssFlag=$tssFlag_ nPeer=$nPeer"
        
        threads_=4
        task_=fastqtl${tssFlag_}.$eX_.$fastqtlTag
        # subp_script_=$subp/docker.fastqtl.v2.wf.sh
        subp_script_=$subp/docker.fastqtl.v3.scratch.sh
        outdir_=$fastqtl_basedir/$task_
        qclogdir_=$outdir_/log
        
        hold_="-hold_jid a -pe smp $threads_"
        qcFpStr_=$sourceFp
        
        chrChar_=$chrChar
        inExpFp_=""
        # Sample used in the eQTL analysis is controlled by the covFp
        #   genotype and expression file can contain more samples than needed
        inExpFpAll_=$zdir_/out.${eX_}${tssFlag_}.expression.bed.gz
        # USU DO NOT CHANGE because this is the output from step one
        inExpFp_=$outdir_/out.-----.expression.bed.gz
        
        # these settings are for each threads requested
        MEM_=4G
        HR_=72:00:00
        
        # split expression file into chromosome (only need to run once and results shared with all)
        step_=prep
        hold_="-hold_jid a"
        qcFpStr_=$sourceFp
        t_=1-23
        nT_=23
        qcStr_="END.${step_}.FINAL"
        
        
        qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$qclogdir_,qctag=.$step_"
        vStr="qcFpStr=$qcFpStr_,inExpFp=$inExpFpAll_,chrChar=$chrChar_,command=,command2=,dockerFp=,rundir=,step=$step_"
        runTask $task_.$step_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" "$vStr"
            
    # each nPeer
    done
# 2 == 0   
fi
    


# Note that for z.9: .tss is skipping one gene for every chr but not for tssFlag_ version
tssFlag_=.tss
# tssFlag_=
rundir_base=/wynton/group/burchard/maka/project/topmedp4.rnaseq/predictdb/$z_${tssFlag_}
covdir_=/wynton/group/burchard/maka/project/topmedp4.rnaseq/gtex.preprocessing/$z_
nPeer=50
fastqtlTag=all
eX_=$e_
switch_=$z_..predictdb
if [ $switch_ == 0 ];
then
    echoFormat "SWITCH:$switch_"
    
    # prepare data for predictdb using per population SNV (MAF >= 0.01)
    if [ 1 == 0 -a $eI_ -ge 4 ];
    then
        
        geneAnnoFullFp=""
        if [ "$tssFlag_" == ".tss" ];
        then
            geneAnnoFullFp=/wynton/group/burchard/maka/data_shared/GRCh38/gencode.v30.GRCh38.annotation.gOnly.tss.bed
        else
        
            geneAnnoFullFp=/wynton/group/burchard/maka/data_shared/GRCh38/gencode.v30.GRCh38.annotation.ERCC.gOnly.bed
        fi
        
        echoFormat "    geneAnnoFullFp=$geneAnnoFullFp"
        
        for fastqtlTag in all;
        do
            rundir=$rundir_base/$fastqtlTag.${eX_}
            task_=predictdb_dataprep.$fastqtlTag.${eX_}.$nPeer${tssFlag_}
            echoTask $task_
            outdir_=$rundir/z/$task_
            qclogdir_=$outdir_/log
            hold_="-hold_jid a"
            qcFpStr_=$sourceFp
            # known bug, 1 gene less from beginning of each chromosome (affect pz9 and hz9)
            subp_script_=$subp/z_geneExp/predictdb_dataprep.v2.sh
            
            covFp=$outdir_/z.cov.ori.txt.gz
            # !!!!! pick the correct number of covariates
            covFp_tmp=$covdir_/peer.1000.80.$eX_/out.$eX_.$nPeer.extract.peer.covariates.txt.gz
            
            if [ ! -d $outdir_ ]; then mkdir -p $outdir_; fi
            
            # This is probably unnessary: original version need to modify original covariate file
            cp $covFp_tmp $covFp
            
            echo $covFp
            vcf_=$zdir_/bcftools_isec.DP0.subset.mafp01.genop05.$eX_/isec.chr-----/0000.vcf.gz
            
            # fastqtl
            expFp=/wynton/group/burchard/maka/project/topmedp4.rnaseq/fastqtl/$z_/fastqtl${tssFlag_}.${eX_}.$fastqtlTag/out.-----.expression.bed.gz
            
            geneAnnoFullFp=$geneAnnoFullFp
            tissue=whole_blood
            batchSize=200
            dbsnpdir=/wynton/group/burchard/maka/data_shared/GRCh38
            dbsnpFI=All_20170710
            geneTypeColName=geneType
            
            t_=1-23
            nT_=23
            vStr="qcFpStr=$qcFpStr_,rfmix_outFp=$rfmix_outFp_,geneAnnoFullFp=$geneAnnoFullFp,tissue=$tissue,covFp=$covFp,vcfFp=${vcf_},expFp=$expFp,eGeneFp=$eGeneFp,batchSize=$batchSize,dbsnpdir=$dbsnpdir,dbsnpFI=$dbsnpFI,rundir=$rundir,geneTypeColName=$geneTypeColName"
            
            MEM_=4G
            HR_=48:00:00
            
            qcStr_="END"
            qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$qclogdir_,qctag="
          
            runTask $task_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" "$vStr"
            
        done
    fi
    
    
    
    # run predictdb_array z.9
    if [ 1 == 0 -a $eI_ -ge 4 ];
    then
        for fastqtlTag in all;
        do
            rundir=$rundir_base/$fastqtlTag.${eX_}
            
            # see below for script to set up this pipeline: ## Test script to set up predictDB_v7 pipeline ##
            task_=predictdb_array.$fastqtlTag.${eX_}
            rFp_=$rundir/model_training/scripts/gtex_tiss_chrom_training_mod.R
            rFp2_=$rundir/model_training/scripts/gtex_v7_nested_cv_elnet.mod.R
            echoTask $task_:$rFp_:$rFp2_
            outdir_=$rundir/z/$task_
            
            general_rFp="/wynton/group/burchard/maka/data_shared/workflow.sh/bin.sub/general.qb3.wf.R"
            qclogdir_=$outdir_/log
            hold_="-hold_jid $jobPre.predictdb_dataprep.$fastqtlTag.${eX_}.qc"
            qcFpStr_=$sourceFp
            subp_script_=$subp/z_geneExp/predictdb_array.wf.sh
            
            # has to be numeric
            chr=-----
            outFI=out.-----.CHRI-----.INDEX-----
            
            rundir2=$rundir
            # need to merge all array
            indexFp_=$rundir2/prepare_data/list.1.23.txt
            gene_annot_file=$rundir2/prepare_data/gene.annotation.parsed.-----.CHRI-----.txt
            genotype_file=$rundir2/prepare_data/genotype.snp.-----.txt
            snp_annot_file=$rundir2/prepare_data/genotype.snp.anno.-----.txt
            expression_file=$rundir2/prepare_data/whole_blood_Analysis.expression.-----.txt
            # covFp is the same for all chromosome unless local ancestry mode
            covFp=$rundir2/prepare_data/whole_blood_Analysis.combined_covariates.1.txt
            
            # Covariates other than these will be converted into numeric in scripts
            covcstr_="asthma:Sex"
            if [ $e_ == "All" ];
            then
                covcstr_="asthma:Sex:e.data_MX:e.data_PR:e.data_LA"
            elif [ $e_ == "3pop" ];
            then
                covcstr_="asthma:Sex:e.data_MX:e.data_PR"
            fi
            echoFormat "    covcstr=$covcstr_"
            
            a=`head -1 $covFp | cut -f2-`
            b=`head -1 $rundir2/prepare_data/whole_blood_Analysis.expression.1.txt | cut -f2-`
            c=`head -1 $rundir2/prepare_data/genotype.snp.1.txt | cut -f2-`
            
            if [ "$a" == "$b" ];
            then
                echo "PASS: sample order 1";
            else
                echo "ERROR 1"
                kill -SIGINT $$
            fi
            
            if [ "$a" == "$c" ];
            then 
                echo "PASS: sample order 2";
            else
                echoFormat "ERROR 2; a=$a b=$b c=$c"
                kill -SIGINT $$
            fi
            
            if [ ! -d $rundir/model_training/scripts ]; then mkdir -p $rundir/model_training/scripts; fi
            if [ ! -f $rundir/model_training/scripts/gtex_tiss_chrom_training_mod.R ];
            then
                ln -s /wynton/group/burchard/maka/data_shared/pipeline/PredictDB_Pipeline_GTEx_v7/model_training/scripts/* $rundir/model_training/scripts/
            fi
            rPar_="$outFI $chr $snp_annot_file $gene_annot_file $genotype_file $expression_file $covFp $rundir $covcstr_ $rFp2_ $general_rFp"
            
            if [ ! -f $indexFp_ -a 3 == 3 ];
            then
                echo "Generating merge list"
                cat $rundir/prepare_data/list.1.txt > $rundir/prepare_data/list.tmp
                for chrI in {2..23};
                do
                    cat $rundir/prepare_data/list.$chrI.txt >> $rundir/prepare_data/list.tmp
                done
                awk 'BEGIN{FS="\t";OFS="\t"}{print FNR,$0}' $rundir/prepare_data/list.tmp > $indexFp_ && rm -rf $rundir/prepare_data/list.tmp
            fi
            
            vStr="qcFpStr=$qcFpStr_,indexFp=$indexFp_,rFp=$rFp_"
            MEM_=30G
            HR_=48:00:00
            
            nLine=`cat $indexFp_ | wc -l`
            t_=1-$nLine
            nT_=$nLine
            
            # t_=1
            # nT_=0
            
            qcStr_="SUCCESS"
            qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$qclogdir_,qctag="
    
            # runTask_par $task_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" "$vStr" "$rPar_" 
            
            taskParFp_=$outdir_/z.parFp
            runTask_par2 $task_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" "$vStr" "$rPar_" $taskParFp_
        done
    fi
    
    
    ## Merge predictdb results z.9 (tss or gene)
    if [ 1 == 0 ];
    then
        # confirmed that construction was not used in the model
        #   PredictionModel.load_model : pipeline/MetaXcan/software/metax/PredictionModel.py
        
        for fastqtlTag in all;
        do
            rundir=$rundir_base/$fastqtlTag.${eX_}
            task_=predictdb_merge.$fastqtlTag.${eX_}
            echoTask $task_
            outdir_=$rundir/z/$task_
            qclogdir_=$outdir_/log
            hold_="-hold_jid $jobPre.predictdb_array.$fastqtlTag.${eX_}.qc"
            qcFpStr_=$sourceFp
            subp_script_=$subp/z_geneExp/predictdb_merge.sh
            
            dir_=$rundir
            tiss_="Whole_Blood"
            rundir2=$rundir
            indexFp_=$rundir2/prepare_data/list.1.23.txt
            
            t_=1
            nT_=1
            nChr_=23
            vStr="qcFpStr=$qcFpStr_,dir=$dir_,tissue=$tiss_,indexFp=$indexFp_,rundir2=$rundir2,nChr=$nChr_"
            
            MEM_=4G
            HR_=48:00:00
            
            qcStr_="SUCCESS"
            qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$qclogdir_,qctag="
      
            runTask $task_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" "$vStr"
            # warning: SQL statements must be issued with dbExecute() or dbSendStatement() instead of dbGetQuery() or dbSendQuery().
            # (You can safely ignore that warning, it's because the first statement returns nothing.) d
        done
    fi
    
    ## Make DB z.9 tss
    if [ 1 == 0 ];
    then
        for fastqtlTag in all;
        do
            rundir=$rundir_base/$fastqtlTag.${eX_}
            task_=runR.make_dbs_mod.$fastqtlTag.${eX_}
            rFp_=$rundir/model_training/scripts/make_dbs_mod.R
            echoTask $task_:$rFp_
            subp_script_=$subp/runR.wf.sh
            outdir_=$rundir/z/$task_
            qclogdir_=$outdir_/log
            hold_="-hold_jid $jobPre.predictdb_merge.$fastqtlTag.${eX_}.qc"
            qcFpStr_=$rundir/z/predictdb_merge.$fastqtlTag.${eX_}/log/z.PASS
            vStr="qcFpStr=$qcFpStr_"
            
            dir_=$rundir/model_training
            pop_=$eX_
            geneAnnoFp_=$rundir/prepare_data/gene.annotation.parsed.1-23.txt
            dbPrefix_=gala.sage.$eX_
            rPar_="$pop_ $geneAnnoFp_ $dir_ $dbPrefix_"
            
            t_=1
            nT_=1
            MEM_=4G
            HR_=48:00:00
            
            qcStr_=SUCCESS.run.$outFI_
            qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$qclogdir_,qctag="
            
            if [ ! -f $rFp_ -o 3 == 3 ];
            then
                ln -s -f /wynton/group/burchard/maka/data_shared/pipeline/PredictDB_Pipeline_GTEx_v7/model_training/scripts/* $rundir/model_training/scripts/
            fi
            
            runR $task_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" $qcFpStr_ $rFp_ "$rPar_"
            
        done
    fi
    
    
    # filter DB z.9 tss
    if [ 1 == 0 ];
    then
        for fastqtlTag in all;
        do
            rundir=$rundir_base/$fastqtlTag.${eX_}
            task_=runR.filter_dbs_mod.$fastqtlTag.${eX_}
            rFp_=$rundir/model_training/scripts/filter_dbs_mod.R
            echoTask $task_:$rFp_
            subp_script_=$subp/runR.wf.sh
            outdir_=$rundir/z/$task_
            qclogdir_=$outdir_/log
            hold_="-hold_jid $jobPre.runR.make_dbs_mod.$fastqtlTag.${eX_}.qc"
            qcFpStr_=$rundir/z/runR.make_dbs_mod.$fastqtlTag.${eX_}/log/z.PASS
            vStr="qcFpStr=$qcFpStr_"
            
            dir_=$rundir/model_training
            dbPrefix_=gala.sage.$eX_
            rPar_="$dir_ $dbPrefix_"
            
            t_=1
            nT_=1
            MEM_=4G
            HR_=48:00:00
            
            qcStr_="SUCCESS"
            qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$qclogdir_,qctag="
            
            if [ ! -f $rFp_ -o 3 == 3 ];
            then
                ln -s -f /wynton/group/burchard/maka/data_shared/pipeline/PredictDB_Pipeline_GTEx_v7/model_training/scripts/* $rundir/model_training/scripts/
            fi
            
            runR $task_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ $qcPar_ $qcFpStr_ "$rFp_" "$rPar_"
        done
    fi
    
    
    ## Create covariance file (needed for S-PrediXcan) z.9 tss
    if [ 1 == 0 ];
    then
        for fastqtlTag in all;
        do
            rundir=$rundir_base/$fastqtlTag.${eX_}
            task_=runCom.create_covariances_mod.$fastqtlTag.${eX_}
            echoTask $task_
            subp_script_=$subp/runCom.wf.sh
            outdir_=$rundir/z/$task_
            qclogdir_=$outdir_/log
            hold_="-hold_jid $jobPre.runR.filter_dbs_mod.$fastqtlTag.${eX_}.qc"
            qcFpStr_=$rundir/z/runR.filter_dbs_mod.$fastqtlTag.${eX_}/log/z.PASS
            vStr="qcFpStr=$qcFpStr_"
            
            dbPrefix_=gala.sage.$eX_
            
            t_=1
            nT_=1
            MEM_=4G
            HR_=48:00:00
            
                # dbPrefix_=gala.sage.AA
                # rundir=/wynton/group/burchard/maka/project/topmedp4.rnaseq/predictdb/run3B.z.8.tssY/all.AA
                # python /wynton/group/burchard/maka/data_shared/pipeline/PredictDB_Pipeline_GTEx_v7/model_training/scripts/create_covariances_mod.py -d $rundir -f $dbPrefix_
            
            command_="python $rundir/model_training/scripts/create_covariances_mod.py -d $rundir -f $dbPrefix_"
            
            qcStr_="SUCCESS"
            qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$qclogdir_,qctag="
            
            if [ ! -f $inScriptFp_ ];
            then
                ln -s -f /wynton/group/burchard/maka/data_shared/pipeline/PredictDB_Pipeline_GTEx_v7/model_training/scripts/* $rundir/model_training/scripts/
            fi
            
            runCom $task_ "$hold_" $outdir_ $MEM_ $HR_ $qcFpStr_ $subp_script_ $t_ $nT_ $qcPar_ "$command_"
        done
    fi
    # predictdb summary statistics z.9
fi # predictdb z.9 tss