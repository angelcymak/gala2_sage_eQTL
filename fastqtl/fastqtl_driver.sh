# --------- z.9: analysis grouped by self identified race/ethnicity
z_=z.9

# ---------- z.9: QC and prepare count file for each population
if [ 1 == 0 ];
then

    # generate count, tpm, inverse normalization and fastqtl gene expression file for each population
    
    if [ 1 == 0 ];
    then
        
        subp_script_=$subp/runR.v2.sh
        rFp_=$wf1/rnaseq_2B_fastqtl_dataprep_topmedf4_z.9.R
        k_=$eI_
        MEM_=8G
        HR_=72:00:00
        
        task_=z.runR.rnaseq_2B_dataprep.$e_
        echoTask $task_
        outdir_=$zdir_/$task_
        qclogdir_=$outdir_/log
        hold_="-hold_jid a"
        qcFpStr_=$sourceFp
        
        outFI_=out.$e_
        rPar_="$k_"
        
        t_=1
        nT_=1
        
        qcStr_=SUCCESS.run.$outFI_
        qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$qclogdir_,qctag="
        runRv2 $task_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" $qcFpStr_ $rFp_ $outFI_ "$rPar_"
    fi
    
    # Generate peer factor z.9
    if [ 1 == 0 -a $eI_ -ge 4 ];
    then
        # double check covariates before running
        module load CBI && module load htslib
        for max_iter_ in 1000;
        do
            for nPeer_ in 80;
            do
                # !!! SAMPLE SCRIPT
                # Make sure each run has a different output directory to a separate folder (rPar_) will be overwritten
                subp_script_=$subp/runR.v2.sh
                rFp_=$wf0/bin.standalone/rnaseq_2C_${z_}_fastqtl_dataprep_topmedf4.R
                task_=z.runR.rnaseq_2C_fastqtl_dataprep_topmedf4.$z_.$e_.$max_iter_.$nPeer_
                outdir_=$zdir_/peer.$max_iter_.$nPeer_.$e_
                outFpI_=$outdir_/out
                inFpI=$zdir_/out
                
                qclogdir_=$outdir_/log
                hold_="-hold_jid $jobPre.z.runR.rnaseq_2B_dataprep.$e_.qc"
                # qcFpStr_=$zdir_/z.runR.rnaseq_2B_dataprep.$e_/log/z.PASS
                qcFpStr_=$zdir_/z.runR.rnaseq_2B_dataprep.$e_/log/status.SUCCESS.run.out.$e_
                
                invnormFp_=$zdir_/out.${e_}.count_matrix.invnorm.rds
                logdir_=$outdir_/log
                statusFp_=$logdir_/status.END.$e_.$max_iter_.$nPeer_
                
                outFI_=$e_
                rPar_="$outFpI_ $inFpI $eI_ $max_iter_ $nPeer_ $invnormFp_ $statusFp_"
                
                MEM_=8G
                HR_=72:00:00
                t_=1
                nT_=1
                
                qcStr_=SUCCESS.run.$outFI_
                qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$qclogdir_,qctag="
                runRv2 $task_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" $qcFpStr_ $rFp_ $outFI_ "$rPar_"
                
                # QC: check covariate in Rout
                
            done
        done
    fi

fi 



#  --------------Prepare genotype file : z.9
switch_=$z_..genotype
if [ $switch_ == 0 ];
then
    echoFormat "SWITCH:$switch_"
    
    # !!!! IMPORTANT: first filter site that are 0.01 and geno 0.05 within samples
    #   To avoid missingness, use the DP10 list to filter out sites for DP0
    
    # prepare genotype file for fastqtl: subset samples only so that stats will be updated
    threads_=4
    
    if [ 1 == 0 ];
    then
        task_=bcftools_view_chr.DP10.subset.mafp01.genop05.$e_
        echoTask $task_
        subp_script_=$subp/bcftools_view_chr.wf.sh
        outdir_=$zdir_/$task_
        qclogdir_=$outdir_/log
        hold_=" -hold_jid a -pe smp $threads_"
        qcFpStr_=$sourceFp
        keepList_=$zdir_/out.$e_.NWD.txt
        inFp_=/wynton/group/burchard/maka/project/topmedf8/DP10/bcftools_dbsnpAnno_array.v3.DP10.xmono.pass.ucsf/chr-----.bcf
        # the AF filtering did not work as expected; need to set max-af to be 0.99
        # filterStr_="--samples-file $keepList_ --min-af 0.01 --include F_MISSING<0.05 -O b"
        
        # Fixed final genotype file in plink using --maf 0.01
        # filterStr_="--samples-file $keepList_ --min-af 0.01 --max-af 0.99 --include F_MISSING<0.05 -O b"
        
        vStr="qcFpStr=$qcFpStr_,inFp=$inFp_,bcftools=$bcftools,threads=$threads_,chrChar=$chrChar"
        MEM_=10G
        HR_=72:00:00
        t_=1-23
        nT_=23
        
        qcStr_="SUCCESS.FINAL"
        qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$qclogdir_,qctag="
        
        # runTask_par $task_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" "$vStr" "$filterStr_"
        
        taskParFp_=$outdir_/z.parFp
        runTask_par2 $task_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" "$vStr" "$filterStr_" $taskParFp_
    fi
    
    # ---------- Independent step
    if [ 1 == 0 ];
    then
        
        # prepare genotype file for fastqtl/predictDB (DP0)
        #   Better to use DP0 than impute
        #   Fastqtl: Missing entries (./., ./0 or ./1) are internally imputed as mean dosage at the variant site.
        #   PredictDB missing not allowed (i believe bcftools dosage will assign 0)
        #   This step can be run in parallel with DP10
        task_=bcftools_view_chr.DP0.subset.$e_
        echoTask $task_
        subp_script_=$subp/bcftools_view_chr.wf.sh
        outdir_=$zdir_/$task_
        qclogdir_=$outdir_/log
        hold_=" -hold_jid a -pe smp $threads_"
        qcFpStr_=$sourceFp
        keepList_=$zdir_/out.$e_.NWD.txt
        inFp_=/wynton/group/burchard/maka/project/topmedf8/DP0/bcftools_view_chr.v3.DP0.xmono.ucsf.mega/out.chr-----.bcf
        filterStr_="--samples-file $keepList_ -O b"
        threads_=2
        vStr="qcFpStr=$qcFpStr_,inFp=$inFp_,bcftools=$bcftools,threads=$threads_,chrChar=$chrChar"
        MEM_=10G
        HR_=72:00:00
        t_=1-23
        nT_=23
        qcStr_="SUCCESS.FINAL"
        qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$qclogdir_,qctag="
        
        # runTask_par $task_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" "$vStr" "$filterStr_"
        
        taskParFp_=$outdir_/z.parFp
        runTask_par2 $task_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" "$vStr" "$filterStr_" $taskParFp_
    fi
    
    if [ 1 == 0 ];
    then
        # Overlap DP0 and DP10
        # DEPENDENCIES:
        #   bcftools_view_chr.DP10.subset.mafp01.genop05.$e_
        #   bcftools_view_chr.DP0.subset.$e_
        task_=bcftools_isec.DP0.subset.mafp01.genop05.$e_
        echoTask $task_
        subp_script_=$subp/bcftools_isec.wf.sh
        outdir_=$zdir_/$task_
        qclogdir_=$outdir_/log
        hold_=" -hold_jid $jobPre.bcftools_view_chr.DP0.subset.$e_.qc,$jobPre.bcftools_view_chr.DP10.subset.mafp01.genop05.$e_.qc -pe smp $threads_"
        
        qcFpStr_=$zdir_/bcftools_view_chr.DP0.subset.All/log/z.PASS:$zdir_/bcftools_view_chr.DP10.subset.mafp01.genop05.$e_/log/z.PASS
        
        inFp1_=$zdir_/bcftools_view_chr.DP0.subset.$e_/out.chr-----.bcf
        inFp2_=$zdir_/bcftools_view_chr.DP10.subset.mafp01.genop05.$e_/out.chr-----.bcf
        filterStr_="$inFp1_ $inFp2_ -n =2 -O z"
        vStr="qcFpStr=$qcFpStr_,bcftools=$bcftools,threads=$threads_,chrChar=$chrChar"
        MEM_=10G
        HR_=72:00:00
        t_=1-23
        nT_=23
        qcStr_="SUCCESS.FINAL"
        qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$qclogdir_,qctag="
        
        # runTask_par $task_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" "$vStr" "$filterStr_"
        
        taskParFp_=$outdir_/z.parFp
        runTask_par2 $task_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" "$vStr" "$filterStr_" $taskParFp_
    fi
fi


#  -------- fastqtl z.9 
fastqtl_basedir=/wynton/group/burchard/maka/project/topmedp4.rnaseq/fastqtl/$z_
if [ ! -d $fastqtl_basedir ]; then mkdir -p $fastqtl_basedir;fi


# prepare TSS version of expression file
if [ 1 == 0 -a $eI_ == 1 ];
then
    dir_=$zdir_
    
    for eX in AA MX PR All;
    do
        # eX=AA
        
        outFp_=$dir_/out.$eX.tss.expression.bed
        zcat $dir_/out.$eX.expression.bed.gz | head -1 > $outFp_
        
        for chr_ in {1..22} X Y M;
        do
            echo $eX $chr_;
            tssFp=/wynton/group/burchard/maka/data_shared/GRCh38/gencode.v30.GRCh38.annotation.$chr_.gOnly.tss.bed
            awk 'BEGIN{FS="\t";OFS="\t"}{
                if (NR == FNR){
                    geneid=$5
                    start0=$2
                    start1=$3
                    array1[geneid]=start0
                    array2[geneid]=start1
                }else{
                    geneid=$4
                    if (geneid in array1){
                        $2=array1[geneid]
                        $3=array2[geneid]
                        print $0
                    }
                }
            }' $tssFp <(zcat $dir_/out.$eX.expression.bed.gz | tail -n +2 | grep "^chr$chr_\s") | sort -k2n -k3n >> $outFp_
        done
            
        bgzip -f $outFp_
        
        tabix -f -p bed $outFp_.gz
        
        zcat $dir_/out.$eX.expression.bed.gz | wc -l
        zcat $outFp_.gz | wc -l
        
    done
    
    # check difference
    if [ 2 == 0 ];
    then
        awk '{
            if (NR == FNR){
                array[$1]=1
            }else{
                if ($1 in array == 0){
                    print
                }
            }
        }' <(zcat $dir_/out.$eX.tss.expression.bed.gz | cut -f4) <(zcat $dir_/out.$eX.expression.bed.gz | cut -f4) > out.$eX.expression.diff.geneid.txt
    fi

fi

# ------------- fastqtl z.9 (tss version)


switch_=$z_..fastqtl
if [ $switch_ == 0 ];
then
    echoFormat "SWITCH:$switch_"
    

    # Generate job list
    jobFp_=$fastqtl_basedir/job.txt
    
    if [ 1 == 0 -a $eI_ == 1 ];
    then
        echo -n > $jobFp_
        
        for nPeer in 50 60 0 10 20 30 40 70 80;
        do
            for eX_ in AA MX PR All;
            do
                echo -e "$nPeer\t$eX_" >> $jobFp_
            done
        done
    fi

    nJob_=`cat $jobFp_ | wc -l`
    
    # Don't change $eI_ == 1 because specific job arrangement is available
    if [ 2 == 0 ];
    then
        minute=0
        minute_step=15
        
        niter=1000
        # all | case | ctrl
        fastqtlTag=all
        tssFlag_=.tss
        
        for i in $(seq 1 $nJob_);
        do
        
            nPeer=`sed -n ${i}p $jobFp_ | awk '{print $1}'`
            eX_=`sed -n ${i}p $jobFp_ | awk '{print $2}'`
            
            # %m = month; %M = minute
            if [ $eX_ == "AA" -a $nPeer != 50 ]; then minute=`echo $minute $minute_step | awk '{print $1+$2}'`; fi;
            submit_time=`date -d "+$minute minute" '+%F %m%d%H%M' | awk '{print $2}'`
            
            runInfo_="fastqtl z=$z_ eX=$eX_ tssFlag=$tssFlag_ nPeer=$nPeer"
            
            threads_=4
            task_=fastqtl${tssFlag_}.$eX_.$fastqtlTag
            # subp_script_=$subp/docker.fastqtl.v2.wf.sh
            subp_script_=$subp/docker.fastqtl.v3.scratch.sh
            outdir_=$fastqtl_basedir/$task_/peer.$nPeer
            outdir2_=$fastqtl_basedir/$task_
            qclogdir_=$outdir_/log
            # wait until the previous job finish, set impossible value to disable
            nPeer_pre=`echo $nPeer | awk '{print $1-90}'`
            
            hold_="-hold_jid a -pe smp $threads_"
            qcFpStr_=$sourceFp
            
            # peer factor in input peer file; use different number of peers from a single run of PEER
            input_peer=80
            chrChar_=$chrChar
            dockerFp_=/wynton/group/burchard/maka/docker/gtex_eqtl_V8.sif
            # this variable is important because where singularity is running is important. It can't access anything from the parent directory
            rundir_=/wynton/group/burchard/maka/project/topmedp4.rnaseq
            inExpFp_=""
            vcf_=""
            # Sample used in the eQTL analysis is controlled by the covFp
            #   genotype and expression file can contain more samples than needed
            covFp=""
            covFp=$zdir_/peer.$niter.$input_peer.$eX_/out.${eX_}.$nPeer.extract.peer.covariates.txt.gz
            inExpFpAll_=$zdir_/out.${eX_}${tssFlag_}.expression.bed.gz
            # USU DO NOT CHANGE because this is the output from step one
            inExpFp_=$outdir2_/out.-----.expression.bed.gz
            # use DP0 to avoid imputation of missing genotypes
            vcf_=$zdir_/bcftools_isec.DP0.subset.mafp01.genop05.$eX_/isec.chr-----/0000.vcf.gz
            
            # --seed
            
            # '--maf_threshold', default='0.0'
            # '--ma_sample_threshold', default='0'
            # gtex: --window 1e6 --ma_sample_threshold 10 --maf_threshold 0.01
            command_="/opt/fastqtl/python/run_FastQTL_threaded.py ${vcf_} $inExpFp_ $outdir_/chr----- --covariates $covFp --window 1e6 --threads $threads_"
            
            
            command2_="/opt/fastqtl/python/run_FastQTL_threaded.py ${vcf_} $inExpFp_ $outdir_/permute/chr----- --covariates $covFp --window 1e6 --threads $threads_ --permute 1000 10000"
            
            
            vStr="qcFpStr=$qcFpStr_,inExpFp=$inExpFp_,chrChar=$chrChar_,command=$command_,command2=$command2_,dockerFp=$dockerFp_,rundir=$rundir_"
            
            if [ 1 == 0 -a \( $eX_ == "3pop"  -o $eX_ == "All" \) ];
            then
                # !!!!! run once before running any steps below
                # Generate PEER covariate file for different number of PEER
                
                echo $eX_ $nPeer
                oriPeerFp=$zdir_/peer.$niter.$input_peer.$eX_/out.${eX_}.$input_peer.peer.covariates.txt.gz 
                
                nLineSkip=`echo $nPeer | awk '{print 80-$1}'`
                zcat $oriPeerFp | head -n -$nLineSkip | bgzip -c > $zdir_/peer.$niter.$input_peer.$eX_/out.${eX_}.$nPeer.extract.peer.covariates.txt.gz 
                
            fi
            
            
            # these settings are for each threads requested
            MEM_=4G
            HR_=72:00:00
            
            if [ 1 == 0 ];
            then
                
                # split expression file into chromosome (only need to run once and results shared with all)
                step_=prep
                hold_="-hold_jid a"
                qcFpStr_=$sourceFp
                t_=1-23
                nT_=23
                qcStr_="END.${step_}.FINAL"
                
                qclogdir_=$outdir2_/log
                reportdir_=$outdir2_/log
                qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$reportdir_,qctag=.$step_"
                vStr="qcFpStr=$qcFpStr_,inExpFp=$inExpFpAll_,chrChar=$chrChar_,command=$command_,command2=$command2_,dockerFp=$dockerFp_,rundir=$rundir_,step=$step_"
                runTask $task_.$step_ "$hold_" $outdir2_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" "$vStr"
            fi
                
                
            # nominal and permute steps (z9)
            # if [ 1 == 0 -a \( $eX_ == "3pop"  -o $eX_ == "All" \) ];
            if [ 1 == 0 ];
            then
                echo $i / $nJob_ e=$eX_ nPeer=$nPeer hr=$hr submit_time=$submit_time;
                echoFormat "$runInfo_"
                echoFormat "\nnominal: $command_"
                echoFormat "\npermutation: $command2_"
                
                if [ ! -d $outdir_ ]; then mkdir -p $outdir_; fi
                
                step_=nominal
                hold_="-hold_jid $jobPre.$task_.prep.qc -a $submit_time"
                qcFpStr_=$outdir2_/log/z.PASS.prep
                t_=1-23
                nT_=23
                qcStr_="END.$step_.FINAL"
                qclogdir_=$outdir_/log
                reportdir_=$outdir_/log
                qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$reportdir_,qctag=.$step_"
                vStr="qcFpStr=$qcFpStr_,inExpFp=$inExpFp_,chrChar=$chrChar_,command=$command_,command2=$command2_,dockerFp=$dockerFp_,rundir=$rundir_,step=$step_"
                # runTask $task_.$step_.$nPeer "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" "$vStr" 
            
                step_=permute
                # hold_="-hold_jid $jobPre.$task_.prep.qc -a $submit_time"
                hold_="-hold_jid $jobPre.$task_.prep.qc"
                qcFpStr_=$outdir2_/log/z.PASS.prep
                t_=1
                nT_=23
                qcStr_="END.$step_.FINAL"
                qclogdir_=$outdir_/permute/log
                reportdir_=$outdir_/log
                qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$reportdir_,qctag=.$step_"
                vStr="qcFpStr=$qcFpStr_,inExpFp=$inExpFp_,chrChar=$chrChar_,command=$command_,command2=$command2_,dockerFp=$dockerFp_,rundir=$rundir_,step=$step_"
                runTask $task_.$step_.$nPeer "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" "$vStr"
                
                
                if [ 2 == 0 ];
                then
                    ls $fz9/*/peer.*/log/*permute*
                    
                    ls */peer.*/log/*permute*
                    # remove log when rerun
                    rm -rf $fz9/*/peer.{10,20,30}/permute
                    rm -rf $fz9/*/peer.{10,20,30}/log/*permute*
                    rm -rf $fz9/*{3pop,All}*/peer.0/log/*permute*
                    
                    echo "ls /scratch/266539.21.long.q/project/topmedp4.rnaseq/fastqtl/z.9/fastqtl.tss.AA.all/peer.30/permute/21*; cat /scratch/266539.21.long.q/project/topmedp4.rnaseq/fastqtl/z.9/fastqtl.tss.AA.all/peer.30/permute/21_chunk022.log" | qsub -l h=qb3-id21
                    
                    # QC
                    # expect 18 per population
                    ls *AA*/peer*/log/z* | wc -l
                    
                    # compare which fail run was due to qvalue error
                    for fp in */peer*/log/*permute*.o*;
                    do
                        nLn_=`grep calculateSignificanceFastQTL.R $fp | wc -l`
                        if [ $nLn_ != 0 ];
                        then
                            echo $fp;
                            dir_=$(dirname "$(dirname "$fp")")
                            ls $dir_/permute/*[1,2,3,4,5,6,7,8,9,0,X].txt.gz
                        fi
                    done
                    # if fail but not see chr.txt.gz, then fail due to error that highly likely to need rerun
                    ls */peer*/log/z.FAIL.permute
                fi
            fi
            
            
            # Merge results z.9
            if [ 1 == 0 ];
            # if [ 1 == 0 -a \( $eX_ == "All" -o $eX_ == "3pop" \) ];
            then
                echoFormat "$runInfo_"
                step_=mergeResults
                hold_="-hold_jid $jobPre.$task_.nominal.$nPeer.qc,$jobPre.$task_.permute.$nPeer.qc"
                qcFpStr_=$outdir_/log/z.PASS.nominal:$outdir_/log/z.PASS.permute 
                t_=1
                nT_=1
                qcStr_="END.$step_.FINAL"
                qclogdir_=$outdir_/log
                reportdir_=$outdir_/log
                qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$reportdir_,qctag=.$step_"
                vStr="qcFpStr=$qcFpStr_,inExpFp=$inExpFp_,chrChar=$chrChar_,command=$command_,command2=$command2_,dockerFp=$dockerFp_,rundir=$rundir_,step=$step_"
                runTask $task_.$step_.$nPeer "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" "$vStr"
                
                # ls docker.fastqtl.{3pop,All}.all/peer*/log
            fi
            
            # if [ 2 == 0 -a $eX_ == "AA" -a $nPeer == 0 ];
            if [ 2 == 0 ];
            then
                echoFormat "$runInfo_"
                subp_script_=$subp/runR.v2.sh
                rFp_=$subp/fastqtl_multipleTesting_v3.R
                task_=runR.fastqtl_multipleTesting_v3.$eX_.$nPeer
                fdr_=0.05
                geneFp_=$outdir_/out.permute.txt.gz
                varFp_=$outdir_/chr-----.allpairs.txt.gz
                MEM_=8G
                HR_=72:00:00
                
                if [ 1 == 0 ];
                then
                    step_=fdr
                    
                    echoTask $task_
                    qclogdir_=$outdir_/log
                    hold_="-hold_jid $jobPre.fastqtl${tssFlag_}.$eX_.$fastqtlTag.mergeResults.$nPeer.qc"
                    qcFpStr_=$outdir_/log/z.PASS.mergeResults
                    
                    outFI_=out.$step_
                    rPar_="$fdr_ $geneFp_ $outdir_ $step_ $varFp_ dummy"
                    rParFp_=$outdir_/z.rPar.$step_.txt
                    
                    t_=1
                    nT_=1
                    
                    qcStr_=SUCCESS.run.out.$step_
                    qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$qclogdir_,qctag=.$step_"
                    # runRv3 $task_.$step_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" $qcFpStr_ $rFp_ $outFI_ "$rPar_" $rParFp_
                    
                    vStr="qcFpStr=$qcFpStr_,outFI=$outFI_"
                    runRv4 $task_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" "$vStr" $rFp_ "$rPar_" $rParFp_
                fi
                
                # merge fdr cutoff to define eqtl
                if [ 1 == 0 ];
                then
                    step_=eqtl
                    echoTask $task_
                    qclogdir_=$outdir_/log
                    hold_="-hold_jid $jobPre.$task_.fdr.qc"
                    qcFpStr_=$outdir_/log/z.PASS.fdr
                    
                    outFI_=out.$step_.-----
                    rPar_="$fdr_ $geneFp_ $outdir_ $step_ $varFp_ chr-----"
                    rParFp_=$outdir_/z.rPar.$step_.txt
                    
                    t_=1-23
                    nT_=23
                    
                    qcStr_=SUCCESS.run.out.$step_
                    qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$qclogdir_,qctag=.$step_"
                    # runRv3 $task_.$step_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" $qcFpStr_ $rFp_ $outFI_ "$rPar_" $rParFp_
                    
                    vStr="qcFpStr=$qcFpStr_,outFI=$outFI_"
                    runRv4 $task_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" "$vStr" $rFp_ "$rPar_" $rParFp_
                fi
            fi
            
        # each job
        done
    # 2 == 0   
    fi
    
    
    eX_=$e_
    
    # Counting number of gene-eQTL pair for each Peer (z.9)
    # if [ 1 == 0 -a $eX_ == "AA" ];
    # if [ 1 == 0 -a \( $eX_ == "All" -o $eX_ == "3pop" \) ];
    if [ 1 == 0 ];
    then
        
        fastqtlDir2_=$fastqtl_basedir/fastqtl${tssFlag_}.$eX_.$fastqtlTag
        outdir_=$fastqtlDir2_
        
        task_=fastqtl_count${tssFlag_}.$eX_.$fastqtlTag
        echoTask $task_
        qclogdir_=$fastqtlDir2_/log
        hold_="-hold_jid $jobPre.runR.fastqtl_multipleTesting_v3.$eX_.$nPeer2.eqtl.qc"
        subp_script_=$subp/fastqtl_count.sh
        
        if [ ! -d $fastqtlDir2_ ]; then mkdir -p $fastqtlDir2_;fi
        
        qcFpStr_=$fastqtlDir2_/peer.PEER-----/log/z.PASS.eqtl
        indexFp_=$fastqtlDir2_/z.peer.txt
        
        if [ ! -d $indexFp_ ];
        then
            echo -n > $indexFp_
            for nPeer2 in 70;
            # for nPeer2 in 0 10 20 30 40 50 60 70 80;
            do
                echo $nPeer2 >> $indexFp_
            done
        fi
        
        nLine_=`cat $indexFp_ | wc -l`
        t_=0
        # t_=1-$nLine_
        # nT_=$nLine_
        nT_=9
        vStr="qcFpStr=$qcFpStr_,pop=$eX_,fastqtlDir=$fastqtlDir2_,indexFp=$indexFp_"
        
        MEM_=4G
        HR_=48:00:00
        
        qcStr_="END.fastqtl_count"
        qcPar_="logdir=$qclogdir_,nExp=$nT_,str=$qcStr_,reportdir=$qclogdir_,qctag=.count"
        
        runTask $task_ "$hold_" $outdir_ $MEM_ $HR_ $subp_script_ $t_ $nT_ "$qcPar_" "$vStr"
        
        # ls -al docker.fastqtl.{3pop,All}.all/log/ | grep fastqtl_count
        # ls -al docker.fastqtl.{3pop,All}.all/z.gv.*
        
    fi
    
    
    
    # write all peer count into one file
    if [ 1 == 0 ];
    then
        echoFormat "Merge PEER stats into on file"
        echoFormat $eX_
        dir_=$fastqtl_basedir
        outFp=$dir_/z.gv.$eX_${tssFlag_}.txt
        
        echo -n > $outFp
        
        for nPeer2 in 0 10 20 30 40 50 60 70 80;
        do
            fp=$dir_/fastqtl${tssFlag_}.$eX_.all/z.gv.$eX_.$nPeer2.txt
            
            nLn=`cat $outFp | wc -l`
            
            if [ -f $fp -a $nLn == 0 ];
            then
                echo $fp header
                cat $fp > $outFp
            elif [ -f $fp -a $nLn -gt 0 ];
            then
                echo $fp append
                tail -n +2 $fp >> $outFp
            fi
        done
    fi
    
fi


