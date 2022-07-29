#!/bin/bash                        #-- what is the language of this shell
#$ -S /bin/bash                    #-- the shell for the job

date
hostname
ulimit -c 0 
echo $JOB_ID

host=`hostname`

logdir=${logdir}
sourceFp=${sourceFp}
generalFp=${generalFp}
source $sourceFp
source $generalFp

chr="${chrList[$SGE_TASK_ID]}"

inFp_raw=${inFp}
bcftools=${bcftools}
chrChar=${chrChar}
threads=${threads}
filterStrFp=${parFp}
filterStr=`sed -n 1p $filterStrFp`

inFp_raw2=`echo $inFp_raw | sed "s|chr-----|$chr|g"`
inFp=`echo $inFp_raw2 | sed "s|-----|$SGE_TASK_ID|g"`

qcFpStr_raw=${qcFpStr}
qcFpStr_raw2=`echo $qcFpStr_raw | sed s/chr-----/$chr/g`
qcFpStr=`echo $qcFpStr_raw2 | sed s/-----/$SGE_TASK_ID/g`
qcFpStr=$qcFpStr:$inFp
qc_tag $qcFpStr $logdir $SGE_TASK_ID

RETVAL=$?
if [ $RETVAL -ne 0 ];
then
  kill -SIGINT $$
fi

rm -rf $logdir/status.HALT.qc.$SGE_TASK_ID

# inFp=/hernandez/netapp/maka/project/topmedf8/DP10/rfmixv2/phasing_BEAGLE.refS/out.22.vcf.gz
# filterStr="-S $sampleFp_$sampleFp_local -O z"

suf=`echo $filterStr | awk '{if ($0 ~ /-O z/){print ".vcf.gz"}else if ($0 ~ /-O b/){print ".bcf"}else if ($0 ~ /-O v/){print ".vcf"} }'`
echo filterStr=$filterStr suf=$suf

outdir=$(dirname $logdir)
outFpI=$outdir/out.$chr
outFp=$outFpI${suf}

if [ ! -d $logdir ]; then mkdir -p $logdir; fi



if [ ! -f $inFp.tbi -a ! -f $inFp.csi ];
then
    start_date=`date +"%Y_%m%d"`
    start_time=`date +"%T"`
    logSuf=indexIn

    echo -n "$logSuf $JOB_ID $start_date $start_time $host $SGE_TASK_ID" > $logdir/status.RUN.$logSuf.$SGE_TASK_ID
    rm -rf $logdir/status.ERROR.$logSuf.$SGE_TASK_ID
    rm -rf $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID

    # -O z: compressed VCF
    $bcftools index $inFp

    RETVAL=$?

    end_date=`date +"%Y_%m%d"`
    end_time=`date +"%T"`

    indexFlag=0
    if [ $RETVAL -ne 0 ];
    then
        echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time ERROR $SGE_TASK_ID" > $logdir/status.ERROR.$logSuf.$SGE_TASK_ID
        kill -SIGINT $$
    else
        echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time SUCCESS $SGE_TASK_ID" > $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID
    fi
fi



start_date=`date +"%Y_%m%d"`
start_time=`date +"%T"`
logSuf=bcftools

echo -n "$logSuf $JOB_ID $start_date $start_time $host $SGE_TASK_ID" > $logdir/status.RUN.$logSuf.$SGE_TASK_ID
rm -rf $logdir/status.ERROR.$logSuf.$SGE_TASK_ID
rm -rf $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID


# -O z: compressed VCF
echo $bcftools
echo $inFp
$bcftools view $filterStr \
  -o $outFp \
  --threads $threads \
  -r $chrChar$chr \
  $inFp

RETVAL=$?

end_date=`date +"%Y_%m%d"`
end_time=`date +"%T"`

indexFlag=0
if [ $RETVAL -ne 0 ];
then
    echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time ERROR $SGE_TASK_ID" > $logdir/status.ERROR.$logSuf.$SGE_TASK_ID
    indexFlag=0
    kill -SIGINT $$
else
    
    echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time SUCCESS $SGE_TASK_ID" > $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID
    
    start_date=`date +"%Y_%m%d"`
    start_time=`date +"%T"`
    nIn=`$bcftools view -H -G $inFp | wc -l`
    touch $outFpI.$nIn.nIn
    
    n=`$bcftools view -H -G $outFp | wc -l`
    touch $outFpI.$n.nOut
    
    end_date=`date +"%Y_%m%d"`
    end_time=`date +"%T"`

    if [ $n -gt 0 ];
    then
        echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time NONEMPTY $SGE_TASK_ID" > $logdir/status.NONEMPTY.$logSuf.$SGE_TASK_ID
        indexFlag=1
    else
        echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time SUCCESS $SGE_TASK_ID" > $logdir/status.SUCCESS.FINAL.$SGE_TASK_ID
        echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time EMPTY $SGE_TASK_ID" > $logdir/status.EMPTY.$logSuf.$SGE_TASK_ID
        indexFlag=0
    fi
fi


if [ $indexFlag == 1 ];
then
    
    indexPar=""
    
    if [ "$suf" == ".vcf" -o "$suf" == ".vcf.gz" ];
    then
        indexPar="-t"
    fi
    
    echo indexPar=$indexPar
    
    start_date=`date +"%Y_%m%d"`
    start_time=`date +"%T"`
    logSuf=FINAL
    
    echo -n "$logSuf $JOB_ID $start_date $start_time $host $SGE_TASK_ID" > $logdir/status.RUN.$logSuf.$SGE_TASK_ID
    rm -rf $logdir/status.ERROR.$logSuf.$SGE_TASK_ID
    rm -rf $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID
    
    # index 
    # -----
    $bcftools index $indexPar $outFp

    RETVAL=$?
    
    end_date=`date +"%Y_%m%d"`
    end_time=`date +"%T"`
    if [ $RETVAL -ne 0 ];
    then
      echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time ERROR $SGE_TASK_ID" > $logdir/status.ERROR.$logSuf.$SGE_TASK_ID
    else
      echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time SUCCESS $SGE_TASK_ID" > $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID
    fi
fi

# put at end of script
jobsum $JOB_ID | grep usage
