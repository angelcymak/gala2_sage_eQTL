#!/bin/bash                        #-- what is the language of this shell
#$ -S /bin/bash                    #-- the shell for the job

date
hostname
ulimit -c 0 
echo $JOB_ID

host=`hostname`

##    it to local /scratch, if it exists, otherwise to /tmp
if [[ -z "$TMPDIR" ]]; then
  if [[ -d /scratch ]]; then TMPDIR=/scratch/$USER; else TMPDIR=/tmp/$USER; fi
  mkdir -p "$TMPDIR"
  export TMPDIR
  echo TMPDIR=$TMPDIR
fi


logdir=${logdir}
sourceFp=${sourceFp}
generalFp=${generalFp}
source $sourceFp
source $generalFp

command_raw=${command}
command2_raw=${command2}
inExpFp=${inExpFp}
chrChar=${chrChar}
# indexFp is the actual index to be given to output file for each SGE_TASK_ID
indexFp=${indexFp}
dockerFp=${dockerFp}
# this variable is important because where singularity is running is important. It can't access anything from the parent directory
rundir_raw=${rundir}
step=${step}



if [ ! -z $indexFp ];
then
    SGE_TASK_ID_i=`sed -n ${SGE_TASK_ID}p $indexFp`
else
    SGE_TASK_ID_i=$SGE_TASK_ID
fi
echo "SGE_TASK_ID=$SGE_TASK_ID SGE_TASK_ID_i=$SGE_TASK_ID_i"
chr="${chrList[$SGE_TASK_ID_i]}"

qcFpStr_raw=${qcFpStr}
qcFpStr_raw2=`echo $qcFpStr_raw | sed "s|chr-----|$chr|g"`
qcFpStr=`echo $qcFpStr_raw2 | sed "s|-----|$SGE_TASK_ID_i|g"`
qcFpStr=$qcFpStr
if [ $step == 'prep' ];
then
    qcFpStr=$qcFpStr:$inExpFp
fi
qc_tag $qcFpStr $logdir $SGE_TASK_ID_i

RETVAL=$?
if [ $RETVAL -ne 0 ];
then
  kill -SIGINT $$
fi
rm -rf $logdir/status.HALT.qc.$SGE_TASK_ID_i


logdirH=$logdir
outdirH=$(dirname "$logdirH")
outdirHPar=$(dirname "$outdirH")
if [ ! -d $outdirH ];then mkdir -p $outdirH; fi
if [ ! -d $logdirH ];then mkdir -p $logdirH; fi

outdir=`echo $outdirH | sed "s|/wynton/group/burchard/maka|$TMPDIR|g"`
if [ ! -d $outdir ];then mkdir -p $outdir; fi

logdir=`echo $logdirH | sed "s|/wynton/group/burchard/maka|$TMPDIR|g"`
if [ ! -d $logdir ];then mkdir -p $logdir; fi

rundir=`echo $rundir_raw | sed "s|/wynton/group/burchard/maka|$TMPDIR|g"`
if [ ! -d $logdir ];then mkdir -p $logdir; fi
echo "rundir=$rundir"
cd $rundir

outdir0=$outdir


outFpI=$outdir/out.$SGE_TASK_ID_i

start_date0=`date +"%Y_%m%d"`
start_time0=`date +"%T"`


# expFp=$outFpI.expression.bed.gz

command_raw2=`echo $command_raw | sed "s|/wynton/group/burchard/maka|$TMPDIR|g"`
command_raw3=`echo $command_raw2 | sed "s|chr-----|$chr|g"`
command=`echo $command_raw3 | sed "s|-----|$SGE_TASK_ID_i|g"`

command2_raw2=`echo $command2_raw | sed "s|/wynton/group/burchard/maka|$TMPDIR|g"`
command2_raw3=`echo $command2_raw2 | sed "s|chr-----|$chr|g"`
command2=`echo $command2_raw3 | sed "s|-----|$SGE_TASK_ID_i|g"`




if [ $step == 'prep' ];
then
    logSuf=$step
    
    if [ ! -f $inExpFp ];
    then
      echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time ERROR $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.ERROR.$logSuf.$SGE_TASK_ID_i
      kill -SIGINT $$
    fi
    
    # subset expression file
    zcat $inExpFp | head -1  > $outFpI.expression.bed
    tabix $inExpFp $chrChar$chr >> $outFpI.expression.bed
    bgzip -f $outFpI.expression.bed
    tabix -f -p bed $outFpI.expression.bed.gz


elif [ $step == 'nominal' ];
then
    
    vcfFp=`echo $command | awk '{print $2}'`
    expFp=`echo $command | awk '{print $3}'`
    covFp=`echo $command | awk '{print $6}'`
    
    vcfdir=$(dirname "$vcfFp")
    expdir=$(dirname "$expFp")
    covdir=$(dirname "$covFp")
    
    if [ ! -d $vcfdir ];then mkdir -p $vcfdir; fi
    if [ ! -d $expdir ];then mkdir -p $expdir; fi
    if [ ! -d $covdir ];then mkdir -p $covdir; fi
    
    vcfFpH=`echo $vcfFp | sed "s|$TMPDIR|/wynton/group/burchard/maka|g"`
    expFpH=`echo $expFp | sed "s|$TMPDIR|/wynton/group/burchard/maka|g"`
    covFpH=`echo $covFp | sed "s|$TMPDIR|/wynton/group/burchard/maka|g"`
    
    rsync -arv $vcfFpH* $vcfdir/
    rsync -arv $expFpH* $expdir/
    rsync -arv $covFpH $covdir/
    
    ls $vcfdir/*
    ls $expdir/*
    ls $covdir/*
    
    # fastqtl nominal pass
    start_date=`date +"%Y_%m%d"`
    start_time=`date +"%T"`
    logSuf=$step
    
    echo -n "$logSuf $JOB_ID $start_date $start_time $host $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.RUN.$logSuf.$SGE_TASK_ID_i
    rm -rf $logdir/status.ERROR.$logSuf.$SGE_TASK_ID_i
    rm -rf $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID_i
    
    echo -e "\n\n\ncommand1=$command"
    singularity exec --bind $(pwd) $dockerFp $command
    
    RETVAL=$?
    
    ls $outdir
    echo outdirHPar=$outdirHPar
    
    end_date=`date +"%Y_%m%d"`
    end_time=`date +"%T"`
    
    if [ $RETVAL -ne 0 ];
    then
      echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time ERROR $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.ERROR.$logSuf.$SGE_TASK_ID_i
      kill -SIGINT $$
    else
      echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time SUCCESS $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID_i
    fi
    
elif [ $step == 'permute' ];
then
    
    vcfFp=`echo $command | awk '{print $2}'`
    expFp=`echo $command | awk '{print $3}'`
    covFp=`echo $command | awk '{print $6}'`
    
    vcfdir=$(dirname "$vcfFp")
    expdir=$(dirname "$expFp")
    covdir=$(dirname "$covFp")
    
    if [ ! -d $vcfdir ];then mkdir -p $vcfdir; fi
    if [ ! -d $expdir ];then mkdir -p $expdir; fi
    if [ ! -d $covdir ];then mkdir -p $covdir; fi
    
    vcfFpH=`echo $vcfFp | sed "s|$TMPDIR|/wynton/group/burchard/maka|g"`
    expFpH=`echo $expFp | sed "s|$TMPDIR|/wynton/group/burchard/maka|g"`
    covFpH=`echo $covFp | sed "s|$TMPDIR|/wynton/group/burchard/maka|g"`
    
    rsync -arv $vcfFpH* $vcfdir/
    rsync -arv $expFpH* $expdir/
    rsync -arv $covFpH $covdir/
    
    ls $vcfdir/*
    ls $expdir/*
    ls $covdir/*
    
    
    outdir2=$outdir/permute
    logdir2=$outdir2/log
    outdir=$outdir2
    logdir=$logdir2
    
    if [ ! -d $logdir2 ];then mkdir -p $logdir2; fi
    
    start_date=`date +"%Y_%m%d"`
    start_time=`date +"%T"`
    logSuf=$step
    
    echo -n "$logSuf $JOB_ID $start_date $start_time $host $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.RUN.$logSuf.$SGE_TASK_ID_i
    rm -rf $logdir/status.ERROR.$logSuf.$SGE_TASK_ID_i
    rm -rf $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID_i
    
    echo -e "\n\n\ncommand2=$command2"
    singularity exec --bind $(pwd) $dockerFp $command2
    
    RETVAL=$?
    
    rsync -arv $outdir0 $outdirHPar/
    echo -e "\n\n cat $outdir/*_chunk022.log"
    cat $outdir/*_chunk010.log
    echo -e "\n\n"
    
    end_date=`date +"%Y_%m%d"`
    end_time=`date +"%T"`
    
    if [ $RETVAL -ne 0 ];
    then
      echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time ERROR $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.ERROR.$logSuf.$SGE_TASK_ID_i
      kill -SIGINT $$
    else
      echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time SUCCESS $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID_i
    fi
    
    
elif [ $step == 'mergeResults' ];
then
    
    # merge all genes for fdr
    logSuf=$step
    
    permFp=$outdir/out.permute.txt
    
    rm -rf $outdirH/permute/*chunk*
    numF_qvalueN=`ls $outdirH/permute/*[0-9X].txt.gz | wc -l`
    numF_qvalueY=`ls $outdirH/permute/*[0-9X].genes.txt.gz | wc -l`
    touch $outdirH/permute/z.numF_qvalueN.$numF_qvalueN
    touch $outdirH/permute/z.numF_qvalueY.$numF_qvalueY
    
    for i in {1..22} X;
    do
        fp=$outdirH/permute/$i.genes.txt.gz
        
        fp2=$outdirH/permute/$i.txt.gz
        
        rm -rf $outdirH/permute/$i.lite.txt.gz
        
        if [ -f $fp ];
        then
            zcat $fp | cut -f1-17 | gzip -c > $outdirH/permute/$i.lite.txt.gz
        fi
        
        if [ -f $fp2 ];
        then
            ln -s -f $fp2 $outdirH/permute/$i.lite.txt.gz
        fi
    done
    
    
    zcat $outdirH/permute/1.lite.txt.gz | head -1 > $permFp
    zcat $outdirH/permute/*.lite.txt.gz | grep -ve "^gene_id"  >> $permFp
    bgzip -f $permFp
    
fi


end_date0=`date +"%Y_%m%d"`
end_time0=`date +"%T"`

echo "FINAL $JOB_ID $start_date0 $start_time0 $host $end_date0 $end_time0 END $SGE_TASK_ID" > $logdir/status.END.$step.FINAL.$SGE_TASK_ID

rsync -arv $outdir0 $outdirHPar/

# put at end of script
jobsum $JOB_ID | grep usage