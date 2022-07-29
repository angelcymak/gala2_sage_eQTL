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

# indexFp is the actual index to be given to output file for each SGE_TASK_ID
indexFp=${indexFp}
parFp=${parFp}

echo $indexFp
echo $parFp

if [ ! -z $indexFp ];
then
    SGE_TASK_ID_i=`sed -n ${SGE_TASK_ID}p $indexFp | awk '{print $1}'`
    chr=`sed -n ${SGE_TASK_ID}p $indexFp | awk '{print $5}'`
    chrI=`sed -n ${SGE_TASK_ID}p $indexFp | awk '{print $2}'`
else
    touch $logdir/status.ERROR.index.$SGE_TASK_ID_i
    echo "Missing $indexFp"
    kill -SIGINT $$
fi
echo "SGE_TASK_ID=$SGE_TASK_ID SGE_TASK_ID_i=$SGE_TASK_ID_i chr=$chr"

parStr_raw=`sed -n 1p $parFp`

parStr_raw2=`echo $parStr_raw | sed "s|INDEX-----|$SGE_TASK_ID_i|g"`
parStr_raw3=`echo $parStr_raw2 | sed "s|CHRI-----|$chrI|g"`
parStr=`echo $parStr_raw3 | sed "s|-----|$chr|g"`

qcFpStr_raw=${qcFpStr}
qcFpStr_raw2=`echo $qcFpStr_raw | sed "s|INDEX-----|$SGE_TASK_ID_i|g"`
qcFpStr_raw3=`echo $qcFpStr_raw2 | sed "s|CHRI-----|$CHRI|g"`
qcFpStr=`echo $qcFpStr_raw3 | sed "s|-----|$chr|g"`
qcFpStr=$qcFpStr:$inFp
qc_tag $qcFpStr $logdir $SGE_TASK_ID_i

RETVAL=$?
if [ $RETVAL -ne 0 ];
then
  kill -SIGINT $$
fi
rm -rf $logdir/status.HALT.qc.$SGE_TASK_ID_i

outdir=$(dirname "$logdir")
if [ ! -d $logdir ];then mkdir -p $logdir; fi
outFI=out.$chr.$chrI.$SGE_TASK_ID_i
outFpI=$outdir/$outFI


start_time=`date`
logSuf=run

echo -n "$logSuf $JOB_ID $start_time $host $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.RUN.$logSuf.$SGE_TASK_ID_i
rm -rf $logdir/status.ERROR.$logSuf.$SGE_TASK_ID_i
rm -rf $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID_i

echo parStr=$parStr
dummy=1
rm -rf $logdir/$outFI.Rout
Rscript --verbose ${rFp} $dummy $parStr > $logdir/$outFI.Rout 2>&1

RETVAL=$?

end_time=`date`
runtime_m=$(( ($(date -d "$end_time" "+%s") - $(date -d "$start_time" "+%s"))/60 ))

if [ $RETVAL -ne 0 ];
then
  echo "$logSuf $JOB_ID $start_time $host $end_time $runtime_m ERROR $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.ERROR.$logSuf.$SGE_TASK_ID_i
  kill -SIGINT $$
else
  echo "$logSuf $JOB_ID $start_time $host $end_time $runtime_m SUCCESS $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID_i
fi

# put at end of script
jobsum $JOB_ID | grep usage