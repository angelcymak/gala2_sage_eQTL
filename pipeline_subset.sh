copyRunFile2(){
    # copyRunFile $parFp $projectParFp $scriptFp $outdir
    parFp_local=$1
    projectParFp_local=$2
    scriptFp_local=$3
    subp_scriptFp_local=$4
    outdir_local=$5
    
    bindir=$outdir_local/bin
    if [ ! -d $bindir ]; then mkdir -p $bindir;fi
    
    scriptF_local=$(basename "$scriptFp_local")
    
    cp $parFp_local $bindir/par_env.txt
    cp $projectParFp_local $bindir/par_project.txt
    cp $scriptFp_local $bindir/
    cp $subp_scriptFp_local $bindir/
    gzip -f $bindir/*
}

copyRunFile3(){
    # copyRunFile3 $outdir {...}
    outdir_=$1
    overwrite=$2
    # .. all files listed after outdir will be put in bin folders
    echo "    "overwrite=$overwrite
    
    bindir=$outdir_/bin
    if [ ! -d $bindir ]; then mkdir -p $bindir;fi
    
    argN=$#
    for i in $(seq 4 $argN);
    do
        fp="${!i}"
        file=$(basename $fp)
        
        fI=0
        outFp=$bindir/$file
        # echo $fp
        
        if [ $overwrite == 0 ];
        then
          while [ -f $outFp.gz ];
          do
            fI=`echo $fI | awk '{print $1+1}'`
            outFp=$bindir/$file.$fI
            # echo $fI $outFp
          done
        fi
        # echo "    "$outFp
        
        cp $fp $outFp
        
        gzip -f $outFp
    done
}

runRv5(){
    task=$1
    hold=$2
    dirStr=$3
    MEM=$4
    HR=$5
    subp_script=$6
    t=$7
    nT=$8
    qcPar=$9
    vStr=${10}
    rFp=${11}
    rPar=${12}
    rParFp=${13}
    
    # rF=$(basename "$rFp")
    # func=${FUNCNAME[0]}
    # task=${func}.$rF
    
    outdir=`echo $dirStr | cut -d ":" -f1`
    qsub_logdir=`echo $dirStr | cut -d ":" -f2`
    logdir=$outdir/log
    
    if [ -z $qsub_logdir ]; then qsub_logdir=$logdir; fi
    
    if [ ! -d $outdir ]; then mkdir -p $outdir; fi
    if [ ! -d $logdir ]; then mkdir -p $logdir; fi
    if [ ! -d $qsub_logdir ]; then mkdir -p $qsub_logdir; fi

    echo -e "\n    ===== SUBMIT TASK: $jobPre.$task"
    echo -e "    outdir=$outdir"
    
    if [ $testrun == 1 ];
    then
        echo "TEST MODE"
        echo -e "      nArg="$#
        argN=$#
        for i in $(seq 1 $argN);
        do
            echo -e "      arg $i="${!i}
        done
        return
    fi
    
    
    echo -n >>  $outdir/qsub.log
    
    if [ -z $rParFp ]; then 
      echo -e "ERROR: empty rParFp ";
      kill -SIGINT $$
    fi
    if [ -f $rParFp ]; then echo -e "WARNING:$rParFp exists";fi
    echo $rPar > $rParFp

    copyRunFile3 $outdir $overwrite $rParFp $projectParFp $scriptFp $subp_script $rFp
    
    if [ $servermode == 1 ];
    then
        echo SERVERMODE
        cat $subp_script | grep -ve "^#\!" > $outdir/bin/script.sh
        chmod 700 $outdir/bin/script.sh
        
        vStr2=`echo $vStr | sed "s/,/ /g"`
        for t_count in $(seq 1 $t);
        do
            echo $t_count of $t
            
            nohup sh -c "SGE_TASK_ID=$t_count JOB_ID=NA sourceFp=$sourceFp generalFp=$generalFp qcFpStr=$qcFpStr logdir=$logdir rFp=$rFp rParFp=$rParFp Rscript=$Rscript $vStr2 $outdir/bin/script.sh" > $logdir/servermode.log 2>&1 &
        done
        
    else
    
        if [ $t != 0 ];
        then
            echo -e "      t=$t"
            qsub $qsubHead -e $qsub_logdir -o $qsub_logdir -N $jobPre.$task $hold -l mem_free=$MEM -l h_rt=$HR -v sourceFp=$sourceFp,generalFp=$generalFp,qcFpStr=$qcFpStr,logdir=$logdir,rFp=$rFp,rParFp=$rParFp,Rscript=$Rscript,"$vStr" -t $t $subp_script >> $outdir/qsub.log
            
        fi
        
        if [ $nT -gt 0 ];
        then
            nExp=$nT
            echo -e "      nT=$nT"
            echo -e "      nExp=$nExp"
            if [ -z $qcStr ];
            then
                str=SUCCESS
            else
                str=$qcStr
            fi
            qsub $qsubHead -e $logdir -o $logdir -N $jobPre.$task.qc -hold_jid $jobPre.$task -l mem_free=1G -l h_rt=00:05:00 -v sourceFp=$sourceFp,generalFp=$generalFp,"$qcPar" $subp/qcFlag.qb3.wf.sh  >> $outdir/qsub.log
        fi
    fi
    
}

runTask(){
    task=$1
    hold=$2
    outdir=$3
    MEM=$4
    HR=$5
    subp_script=$6
    t=$7
    nT=$8
    qclogdir=$9
    vStr=${10}
    waitFp=${11}
    sleeptime=${12}
    qcStr=${13}
    
    echo -e "\n===== SUBMIT TASK: $task"
    echo -e "      outdir=$outdir"
    echo -e "      Job: $jobPre.$task"
    echo -e "      jobPre=$jobPre"
    echo -e "      task=$task"
    
    if [ $testrun == 1 ];
    then
        echo "TEST MODE"
        echo -e "      nArg="$#
        argN=$#
        for i in $(seq 1 $argN);
        do
            echo -e "      arg $i="${!i}
        done
        return
    fi
    
    if [ ! -d $outdir ]; then mkdir -p $outdir; fi
    echo  >>  $outdir/qsub.log
    
    copyRunFile2 $parFp $projectParFp $scriptFp $subp_script $outdir
    
    if [ ! -z $waitFp ];
    then
        while [ ! -f $waitFp ];
        do
          echo waiting for $waitFp
          echo sleep $sleeptime
          sleep $sleeptime
        done
    fi
    
    logdir=$outdir/log
    if [ ! -d $logdir ]; then mkdir -p $logdir; fi
    
    RETVAL=0
    if [ $t != 0 ];
    then
        echo -e "      t=$t"
        echo qsub $qsubHead -e $logdir -o $logdir -N $jobPre.$task $hold
       
        qsub $qsubHead -e $logdir -o $logdir -N $jobPre.$task $hold -l mem_free=$MEM -l h_rt=$HR -v sourceFp=$sourceFp,generalFp=$generalFp,logdir=$logdir,$vStr -t $t $subp_script >> $outdir/qsub.log
        
        RETVAL=$?
    fi
    
    if [ $nT -gt 0 -a $RETVAL == 0 ];
    then
        nExp=$nT
        echo -e "      nT=$nT"
        echo -e "      nExp=$nExp"
        if [ -z $qcStr ];
        then
            str=SUCCESS
        else
            str=$qcStr
        fi
        qsub $qsubHead -e $logdir -o $logdir -N $jobPre.$task.qc -hold_jid $jobPre.$task -l mem_free=1G -l h_rt=48:00:00 -v sourceFp=$sourceFp,generalFp=$generalFp,logdir=$qclogdir,nExp=$nExp,str=$str,reportdir=$logdir $subp/qcFlag.qb3.wf.sh  >> $outdir/qsub.log
        
    fi
    
}

runTaskv2(){
    task=$1
    hold=$2
    dirStr=$3
    MEM=$4
    HR=$5
    subp_script=$6
    t=$7
    nT=$8
    qcPar=$9
    vStr=${10}
    
    outdir=`echo $dirStr | cut -d ":" -f1`
    qsub_logdir=`echo $dirStr | cut -d ":" -f2`
    logdir=$outdir/log
    
    if [ -z $qsub_logdir ]; then qsub_logdir=$logdir; fi
    
    if [ ! -d $outdir ]; then mkdir -p $outdir; fi
    if [ ! -d $logdir ]; then mkdir -p $logdir; fi
    if [ ! -d $qsub_logdir ]; then mkdir -p $qsub_logdir; fi
    
    echo -e "\n    ===== SUBMIT TASK: $jobPre.$task"
    echo -e "    outdir=$outdir"
    
    if [ $testrun == 1 ];
    then
        echo "TEST MODE"
        echo -e "      nArg="$#
        argN=$#
        for i in $(seq 1 $argN);
        do
            echo -e "      arg $i="${!i}
        done
        return
    fi
    
    echo  >>  $outdir/qsub.log
    
    copyRunFile3 $outdir $overwrite $parFp $projectParFp $scriptFp $subp_script
    
    if [ ! -z $waitFp ];
    then
        while [ ! -f $waitFp ];
        do
          echo waiting for $waitFp
          echo sleep $sleeptime
          sleep $sleeptime
        done
    fi
    
    if [ $servermode == 1 ];
    then
        echo SERVERMODE
        cat $subp_script | grep -ve "^#\!" > $outdir/bin/script.sh
        chmod 700 $outdir/bin/script.sh
        
        vStr2=`echo $vStr | sed "s/,/ /g"`
        for t_count in $(seq 1 $t);
        do
            echo $t_count of $t
            
            nohup sh -c "SGE_TASK_ID=$t_count JOB_ID=NA sourceFp=$sourceFp generalFp=$generalFp logdir=$logdir $vStr2 $outdir/bin/script.sh" > $logdir/servermode.log 2>&1 &
        done
        
    else
    
        RETVAL=0
        if [ $t != 0 ];
        then
            echo -e "      t=$t"
           
            qsub $qsubHead -e $qsub_logdir -o $qsub_logdir -N $jobPre.$task $hold -l mem_free=$MEM -l h_rt=$HR -v sourceFp=$sourceFp,generalFp=$generalFp,logdir=$logdir,"$vStr" -t $t $subp_script >> $outdir/qsub.log
            
            # parallel: mem_free is for each core requested ($NSLOTS = $ncores)
            # qsub $qsubHead -e $logdir -o $logdir -N $jobPre.$task  -hold_jid $hold -pe smp $ncores -l mem_free=4G -l h_rt=48:00:00 -v sourceFp=$sourceFp,generalFp=$generalFp,logdir=$logdir,$vStr -t $t $subp/plink.operation.more.wf.sh >> $outdir/qsub.log
            
            RETVAL=$?
        fi
        
        
        if [ $nT -gt 0 -a $RETVAL == 0 ];
        then
            nExp=$nT
            echo -e "      nT=$nT"
            echo -e "      nExp=$nExp"
            if [ -z $qcStr ];
            then
                str=SUCCESS
            else
                str=$qcStr
            fi
            
            qsub $qsubHead -e $logdir -o $logdir -N $jobPre.$task.qc -hold_jid $jobPre.$task -l mem_free=1G -l h_rt=48:00:00 -v sourceFp=$sourceFp,generalFp=$generalFp,"$qcPar" $subp/qcFlag.qb3.wf.sh  >> $outdir/qsub.log
            
        fi
    fi
    
}


