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

dir=${dir}
rundir2=${rundir2}
tissue=${tissue}
indexFp=${indexFp}
nChr=${nChr}

if [ 1 == 0 ];
then
    dir=/wynton/group/burchard/maka/project/topmedp4.rnaseq/predictdb/run1/all.AA
    tissue="Whole_Blood"
    indexFp=$dir/prepare_data/list.1.23.txt
    logdir=$dir/z/predictdb_merge.all.AA
fi

qcFpStr=${qcFpStr}
qcFpStr=$qcFpStr
qc $qcFpStr $logdir

RETVAL=$?
if [ $RETVAL -ne 0 ];
then
  kill -SIGINT $$
fi
rm -rf $logdir/status.HALT.qc

outdir=$(dirname "$logdir")
if [ ! -d $logdir ];then mkdir -p $logdir; fi
outFpI=$outdir/out


start_time=`date`
logSuf=run

echo -n "$logSuf $JOB_ID $start_time $host" > $logdir/status.RUN.$logSuf

rm -rf $logdir/status.ERROR.$logSuf
rm -rf $logdir/status.SUCCESS.$logSuf

geneAnnoFp_raw=$rundir2/prepare_data/gene.annotation.parsed.-----.txt
geneAnnoFp=`echo $geneAnnoFp_raw | sed "s|-----|1|g" `
geneAnnoAllFp=$dir/prepare_data/gene.annotation.parsed.1-23.txt
if [ ! -d $dir/prepare_data ]; then mkdir -p $dir/prepare_data;fi
cat $geneAnnoFp > $geneAnnoAllFp
nGeneAnnoSource=`cat $geneAnnoFp | wc -l`
for chr in $(seq 2 $nChr);
do
    geneAnnoFp=`echo $geneAnnoFp_raw | sed "s|-----|$chr|g" `
    tail -n +2 $geneAnnoFp >> $geneAnnoAllFp
    nTmp=`tail -n +2 $geneAnnoFp | wc -l`
    nGeneAnnoSource=`echo -e "$nGeneAnnoSource\t$nTmp" | awk '{print $1+$2}'`
done
nGeneAnnoFinal=`cat $geneAnnoAllFp | wc -l`

n=`cat $indexFp | wc -l`

line=`sed -n 1p $indexFp`
chr=`echo $line | awk '{print $4}'`
chrI=`echo $line | awk '{print $2}'`
t=`echo $line | awk '{print $1}'`
inFI=out.$chr.$chrI.${t}_nested_cv_chr$chr
outFI=${tiss}_nested_cv_chr$chr

modSumFp=${inFI}_model_summaries.txt
tisSumFp=${inFI}_tiss_chr_summary.txt
weightFp=${inFI}_weights.txt
covarFp=${inFI}_covariances.txt

modSumOutFp=${outFI}_model_summaries.txt
tisSumOutFp=${outFI}_tiss_chr_summary.txt
weightOutFp=${outFI}_weights.txt
covarOutFp=${outFI}_covariances.txt

nModSumSource=0
nTisSumSource=0
nWeightSource=0
nCovarSource=0
for j in $(seq 1 $n);
do
    echo $j
    line=`sed -n ${j}p $indexFp`
    chr=`echo $line | awk '{print $5}'`
    chrI=`echo $line | awk '{print $2}'`
    t=`echo $line | awk '{print $1}'`
    inFI=out.$chr.$chrI.${t}_nested_cv_chr$chr
    outFI=${tissue}_nested_cv_chr$chr
    
    modSumFp=$dir/model_training/summary/${inFI}_model_summaries.txt
    tisSumFp=$dir/model_training/summary/${inFI}_tiss_chr_summary.txt
    weightFp=$dir/model_training/weights/${inFI}_weights.txt
    covarFp=$dir/model_training/covariances/${inFI}_covariances.txt
    
    modSumOutFp=$dir/model_training/summary/${outFI}_model_summaries.txt
    tisSumOutFp=$dir/model_training/summary/${outFI}_tiss_chr_summary.txt
    weightOutFp=$dir/model_training/weights/${outFI}_weights.txt
    covarOutFp=$dir/model_training/covariances/${outFI}_covariances.txt
    
    if [ $chrI == 0 ];
    then
        cat $modSumFp > $modSumOutFp
        cat $tisSumFp > $tisSumOutFp
        cat $weightFp > $weightOutFp
        cat $covarFp > $covarOutFp
        nModSumSource=`tail -n +2 $modSumFp | wc -l | xargs echo $nModSumSource | awk '{print $1+$2}'`
        nTisSumSource=`tail -n +2 $tisSumFp | wc -l | xargs echo $nTisSumSource | awk '{print $1+$2}'`
        nWeightSource=`tail -n +2 $weightFp | wc -l | xargs echo $nWeightSource | awk '{print $1+$2}'`
        nCovarSource=`tail -n +2 $covarFp | wc -l | xargs echo $nCovarSource | awk '{print $1+$2}'`
    else
        tail -n +2 $modSumFp >> $modSumOutFp
        tail -n +2 $tisSumFp >> $tisSumOutFp
        tail -n +2 $weightFp >> $weightOutFp
        tail -n +2 $covarFp >> $covarOutFp
        nModSumSource=`tail -n +2 $modSumFp | wc -l | xargs echo $nModSumSource | awk '{print $1+$2}'`
        nTisSumSource=`tail -n +2 $tisSumFp | wc -l | xargs echo $nTisSumSource | awk '{print $1+$2}'`
        nWeightSource=`tail -n +2 $weightFp | wc -l | xargs echo $nWeightSource | awk '{print $1+$2}'`
        nCovarSource=`tail -n +2 $covarFp | wc -l | xargs echo $nCovarSource | awk '{print $1+$2}'`
    fi
done

nModSumFinal=`cat $dir/model_training/summary/${tissue}_nested_cv_chr*_model_summaries.txt | wc -l | awk -v nChr=$nChr '{ print $1 - nChr}'`

nTisSumFinal=`cat $dir/model_training/summary/${tissue}_nested_cv_chr*_tiss_chr_summary.txt | wc -l | awk -v nChr=$nChr '{ print $1 - nChr}'`

nWeightFinal=`cat $dir/model_training/weights/${tissue}_nested_cv_chr*_weights.txt | wc -l | awk -v nChr=$nChr '{print $1- nChr}'`

nCovarFinal=`cat $dir/model_training/covariances/${tissue}_nested_cv_chr*_covariances.txt | wc -l | awk -v nChr=$nChr '{ print $1- nChr}'`


echo -e "geneAnno $nGeneAnnoSource $nGeneAnnoFinal"
echo -e "modSum $nModSumSource $nModSumFinal"
echo -e "tisSum $nTisSumSource $nTisSumFinal"
echo -e "weight $nWeightSource $nWeightFinal"
echo -e "covar $nCovarSource $nCovarFinal"

end_time=`date`
runtime_m=$(( ($(date -d "$end_time" "+%s") - $(date -d "$start_time" "+%s"))/60 ))

if [ $nGeneAnnoSource == $nGeneAnnoFinal -a $nModSumSource == $nModSumFinal -a $nTisSumSource == $nTisSumFinal -a $nWeightSource == $nWeightFinal -a $nCovarSource == $nCovarFinal ];
then
    echo "$logSuf $JOB_ID $start_time $host $end_time $runtime_m SUCCESS" > $logdir/status.SUCCESS.$logSuf
else
    echo "$logSuf $JOB_ID $start_time $host $end_time $runtime_m ERROR" > $logdir/status.ERROR.$logSuf
    echo "ERROR one of the above number did not match"
fi

# put at end of script
qstat -j $JOB_ID | grep usage
