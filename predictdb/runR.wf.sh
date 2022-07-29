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

rFp=${rFp}
rParFp=${rParFp}
Rscript=${Rscript}

chr="${chrList[$SGE_TASK_ID]}"

qcFpStr_raw=${qcFpStr}
qcFpStr_raw2=`echo $qcFpStr_raw | sed "s/chr-----/$chr/g"`
qcFpStr=`echo $qcFpStr_raw2 | sed "s/-----/$SGE_TASK_ID/g"`
qcFpStr=$qcFpStr
qc_tag $qcFpStr $logdir $SGE_TASK_ID

RETVAL=$?
if [ $RETVAL -ne 0 ];
then
  touch $logdir/status.NOTRUN.qc.$SGE_TASK_ID
  kill -SIGINT $$
fi

echo "checkpoint"
outdir=$(dirname "$logdir")
if [ ! -d $logdir ];then mkdir -p $logdir; fi
outFI=out.$SGE_TASK_ID

rPar_raw=`sed -n 1p $rParFp`
rPar_raw2=`echo $rPar_raw | sed "s|chr-----|$chr|g"`
rPar=`echo $rPar_raw2 | sed "s|-----|$SGE_TASK_ID|g"`
echo $rParFp $rPar_raw $rPar

start_date=`date +"%Y_%m%d"`
start_time=`date +"%T"`
logSuf=run

echo -n "$logSuf $JOB_ID $start_date $start_time $host $SGE_TASK_ID" > $logdir/status.RUN.$logSuf.$SGE_TASK_ID

rm -rf $logdir/status.ERROR.$logSuf.$SGE_TASK_ID
rm -rf $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID


if [[ -z "$TMPDIR" ]]; then
  if [[ -d /scratch ]]; 
  then 
    TMPDIR=/scratch/$USER; 
  else 
    TMPDIR=/tmp/$USER; 
  fi
  mkdir -p "$TMPDIR"
  export TMPDIR
fi

## 1. Use a temporary working directory
# cd "$TMPDIR"


echo "Run Rscript $rFp $rPar"
$Rscript --verbose ${rFp} $rPar > $logdir/$outFI.Rout 2>&1

RETVAL=$?

end_date=`date +"%Y_%m%d"`
end_time=`date +"%T"`

if [ $RETVAL -ne 0 ];
then
  echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time ERROR $SGE_TASK_ID" > $logdir/status.ERROR.$logSuf.$SGE_TASK_ID
  kill -SIGINT $$
else
  echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time SUCCESS $SGE_TASK_ID" > $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID
fi

# put at end of script
jobsum $JOB_ID | grep usage