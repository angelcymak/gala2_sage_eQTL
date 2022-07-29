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

commandFp_raw=${commandFp}

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

outdir=$(dirname "$logdir")
if [ ! -d $logdir ];then mkdir -p $logdir; fi
outFI=out.$SGE_TASK_ID
outFpI=$outdir/$outFI

echo outdir=$outdir

commandFp=$commandFp.$SGE_TASK_ID.final.sh
cat $commandFp_raw | sed "s|chr-----|$chr|g" > ${commandFp_raw}.$SGE_TASK_ID.2.sh
cat ${commandFp_raw}.$SGE_TASK_ID.2.sh | sed "s|outdir-----|$outdir|g" > ${commandFp_raw}.$SGE_TASK_ID.3.sh
cat ${commandFp_raw}.$SGE_TASK_ID.3.sh | sed "s|ls-----|$SGE_TASK_ID|g" > $commandFp

echo -e "\ncommandFp=$commandFp\n"
cat $commandFp
echo
chmod 750 $commandFp

start_date=`date +"%Y_%m%d"`
start_time=`date +"%T"`
logSuf=run

echo -n "$logSuf $JOB_ID $start_date $start_time $host $SGE_TASK_ID" > $logdir/status.RUN.$logSuf.$SGE_TASK_ID

rm -rf $logdir/status.ERROR.$logSuf.$SGE_TASK_ID
rm -rf $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID

$commandFp

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