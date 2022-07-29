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

bcftools=${bcftools}
chrChar=${chrChar}
threads=${threads}
filterStrFp=${parFp}
filterStr_raw=`sed -n 1p $filterStrFp`

filterStr_raw2=`echo $filterStr_raw | sed s/chr-----/$chr/g`
filterStr=`echo $filterStr_raw2 | sed s/-----/$SGE_TASK_ID/g`

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

echo filterStr=$filterStr

outdir=$(dirname $logdir)
outFpI=$outdir/out.$chr
outFp=$outFpI${suf}

if [ ! -d $logdir ]; then mkdir -p $logdir; fi

start_date=`date +"%Y_%m%d"`
start_time=`date +"%T"`
logSuf=bcftools

echo -n "$logSuf $JOB_ID $start_date $start_time $host $SGE_TASK_ID" > $logdir/status.RUN.$logSuf.$SGE_TASK_ID
rm -rf $logdir/status.ERROR.$logSuf.$SGE_TASK_ID
rm -rf $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID


# -O z: compressed VCF
$bcftools isec -p $outdir/isec.$chr $filterStr \
  --threads $threads \
  -r $chrChar$chr

RETVAL=$?

end_date=`date +"%Y_%m%d"`
end_time=`date +"%T"`

rm -rf $logdir/status.ERROR.$logSuf.$SGE_TASK_ID
rm -rf $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID

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
    
    f=`ls $outdir/isec.$chr/0000.* | head -1 `
    # will not work if use 0000.* 
    n=`$bcftools view -H -G $f | wc -l`
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
    
    indexPar="-t"
    
    echo indexPar=$indexPar
    
    start_date=`date +"%Y_%m%d"`
    start_time=`date +"%T"`
    logSuf=FINAL
    
    echo -n "$logSuf $JOB_ID $start_date $start_time $host $SGE_TASK_ID" > $logdir/status.RUN.$logSuf.$SGE_TASK_ID
    
    rm -rf $logdir/status.ERROR.$logSuf.$SGE_TASK_ID
    rm -rf $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID
    
    # index 
    # -----
    
    listFp2=$outdir/isec.$chr.list.txt
    ls $outdir/isec.$chr/0* > $listFp2
    
    nList=`cat $listFp2  | wc -l`
    errorFlag=0
    for k in $(seq 1 $nList);
    do
        fp=`sed -n ${k}p $listFp2`
        suf=`echo $fp | awk '{n=split($1,array,"."); print array[n]}'
        if [ $suf == "vcf" ];
        then
            bgzip $fp
            fp=$fp.gz
        fi
        $bcftools index -t $fp
        
        RETVAL=$?
        if [ $RETVAL -ne 0 ];
        then
            errorFlag=1
        fi
    done

    end_date=`date +"%Y_%m%d"`
    end_time=`date +"%T"`
    if [ $errorFlag -ne 0 ];
    then
      echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time ERROR $SGE_TASK_ID" > $logdir/status.ERROR.$logSuf.$SGE_TASK_ID
      kill -SIGINT $$
    else
      echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time SUCCESS $SGE_TASK_ID" > $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID
    fi
fi

# put at end of script
jobsum $JOB_ID | grep usage
