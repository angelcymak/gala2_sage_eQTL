#!/bin/bash                        #-- what is the language of this shell
#$ -S /bin/bash                    #-- the shell for the job

# NOTES ABOUT IMPLEMENTATION
#   - SNP only
#   - Note that if input genotype file is not exact same number of sample in 
#   sampeFp, all variants with mac 1 will be kept



# 03/10/2021 06:20:22 PM
#   Use gene start +- 1MB flanking

date
hostname
ulimit -c 0 
echo $JOB_ID

host=`hostname`

# Add a new mechanism for rerun
# contain a list of SGE_TASK_ID that should be excluded
excludeSGEFp=${excludeSGEFp}

excludeFlag=`cat $excludeSGEFp | grep ^$SGE_TASK_ID$ | wc -l`

echo excludeSGEFp=$excludeSGEFp
echo excludeFlag=$excludeFlag

if [ $excludeFlag != 0 ];
then
    echo "THIS JOB IS SKIPPED"
    kill -SIGINT $$
fi



## 0. In case TMPDIR is not set, e.g. on development nodes, set
##    it to local /scratch, if it exists, otherwise to /tmp
if [[ -z "$TMPDIR" ]]; then
  if [[ -d /scratch ]]; then TMPDIR=/scratch/$USER; else TMPDIR=/tmp/$USER; fi
  mkdir -p "$TMPDIR"
  export TMPDIR
fi


logdir=${logdir}

if [ ! -d $logdir ];then mkdir -p $logdir; fi

sourceFp=${sourceFp}
generalFp=${generalFp}
source $sourceFp
source $generalFp

# indexFp is the actual index to be given to output file for each SGE_TASK_ID
indexFp_raw=${indexFp}
genoFp_raw=${genoFp}
expFp_raw=${expFp}
parFp=${parFp}
sampleFp_raw=${sampleFp}
save=${save}
unconstrainFlag=${unconstrainFlag}
peerFp=${peerFp}
chrCol=${chrCol}
geneIdCol=${geneIdCol}
geneNameCol=${geneNameCol}
startCol=${startCol}
stopCol=${stopCol}
bcfChrPre=chr
nPeer=${nPeer}
rFp_geneExpRes=${rFp_geneExpRes}

indexFp_lvl1=`echo $indexFp_raw | sed "s|-----|$SGE_TASK_ID|g"`
echo indexFp_lvl1:
head -2 $indexFp_lvl1

echo indexFp_lvl1=$indexFp_lvl1
echo genoFp_raw=$genoFp_raw
echo expFp_raw=$expFp_raw
echo parFp=$parFp
echo chrCol=$chrCol
echo geneIdCol=$geneIdCol
echo geneNameCol=$geneNameCol
echo startCol=$startCol
echo stopCol=$startCol
echo peerFp=$peerFp
echo rFp_geneExpRes=$rFp_geneExpRes

SGE_TASK_ID_i=$SGE_TASK_ID
indexFp=`sed -n ${SGE_TASK_ID}p $indexFp_lvl1 | awk '{print $6}'`
chrNum=`sed -n ${SGE_TASK_ID}p $indexFp_lvl1 | awk '{print $5}'`
chrStr="${chrList[$chrNum]}"
echo "SGE_TASK_ID=$SGE_TASK_ID SGE_TASK_ID_i=$SGE_TASK_ID_i chr=$chr"

qcFpStr_raw=${qcFpStr}
qcFpStr=`echo $qcFpStr_raw | sed "s|INDEX-----|$SGE_TASK_ID_i|g"`
qcFpStr=`echo $qcFpStr_raw3 | sed "s|-----|$chr|g"`


genoFp_raw2=`echo $genoFp_raw | sed "s|CHR-----|$chrStr|g"`
genoFp=`echo $genoFp_raw2 | sed "s|-----|$chrNum|g"`
expFp=`echo $expFp_raw | sed "s|-----|$chrNum|g"`
echo genoFp=$genoFp
echo expFp=$expFp

qcFpStr=$qcFpStr:$genoFp:$expFp
qc_tag $qcFpStr $logdir $SGE_TASK_ID_i

RETVAL=$?
if [ $RETVAL -ne 0 ];
then
  kill -SIGINT $$
fi
rm -rf $logdirH/status.HALT.qc.$SGE_TASK_ID_i


logdirH=$logdir
outdirH=$(dirname "$logdirH")
outdirHPar=$(dirname "$outdirH")
outfolder=$SGE_TASK_ID_i
outdirH2=$outdirH/$outfolder
logdirH2=$logdirH/$outfolder

echo TMPDIR=$TMPDIR
echo logdirH=$logdirH
echo outdirH=$outdirH
echo outdirHPar=$outdirHPar
echo outdirH2=$outdirH2
echo logdirH2=$logdirH2

if [ ! -d $logdirH ];then mkdir -p $logdirH; fi
if [ ! -d $outdirH ];then mkdir -p $outdirH; fi
if [ ! -d $outdirH2 ];then mkdir -p $outdirH2; fi
if [ ! -d $logdirH2 ];then mkdir -p $logdirH2; fi


logdir=`echo $logdirH | sed "s|/wynton/group/burchard/maka|$TMPDIR|g"`
outdir=`echo $outdirH2 | sed "s|/wynton/group/burchard/maka|$TMPDIR|g"`


genoFpH=$genoFp
genoFp=`echo $genoFpH | sed "s|/wynton/group/burchard/maka|$TMPDIR|g"`
genoFp_dir=$(dirname "$genoFp")

echo logdir=$logdir
echo outdir=$outdir
echo genoFpH=$genoFpH
echo genoFp=$genoFp
if [ ! -d $logdir ];then mkdir -p $logdir; fi
if [ ! -d $outdir ];then mkdir -p $outdir; fi
if [ ! -d $genoFp_dir ];then mkdir -p $genoFp_dir; fi

rsync -arv $genoFpH* $genoFp_dir/
if [ ! -f $genoFp ];
then
    touch $logdirH/status.ERROR.genoFp.$SGE_TASK_ID_i
    kill -SIGINT $$
fi


par1=`sed -n 1p $parFp`
par2=`sed -n 2p $parFp`
par3=

if [ $unconstrainFlag == 1 ];
then
    par3='--reml-no-constrain'
fi

echo par1=$par1
echo par2=$par2
echo par3=$par3


module load plink
module load gcta

# --covar : categorical
# --qcovar : quantitative

rm -rf $outdirH2/out.log.txt
rm -rf $outdir/out.converged.allVar.txt
rm -rf $outdir/out.converged.allVar.cov.txt
rm -rf $outdir/out.converged.allVar.cov2.txt
rm -rf $outdir/out.converged.allVar.res.cov2.txt
rm -rf $outdir/out.converged.txt
rm -rf $outdir/out.converged.cov.txt
rm -rf $outdir/out.converged.cov2.txt
rm -rf $outdir/out.converged.res.cov2.txt


meltResults(){
    # $outFpI.gcta.reml.hsq
    resultFp=$1
    anno=$2
    outFp2=$3
    awk -v anno=$anno 'BEGIN{FS="\t";OFS="\t"}{
        if (FNR == 1){
            printf $0"\t"anno
        }else{
            printf "\t"$0
        }
    }END{printf "\n"}' $resultFp >> $outFp2
}

start_date0=`date +"%Y_%m%d"`
start_time0=`date +"%T"`

touch $logdirH/status.RUN.$outfolder
rm -rf $logdirH/status.END.$outfolder


nLn=`cat $indexFp | wc -l`
# For each gene
saveCounter=0
for i in $(seq 2 $nLn);
do
    # i=5
    echo -e "\n\n $i / $nLn"
    
    echo line:
    sed -n ${i}p $indexFp

    saveCounter=`echo $saveCounter | awk '{print $1+1}'`
    chr=`sed -n ${i}p $indexFp | awk -v col=$chrCol '{print $col}'`
    chrStr2="${chrList[$chr]}"
    
    # now use tss
    start=`sed -n ${i}p $indexFp | awk -v col=$startCol '{pos=$col-1000000; if (pos < 0){print 0}else{print pos}}'`
    stop=`sed -n ${i}p $indexFp | awk -v col=$stopCol '{print $col+1000000}'`
    geneid=`sed -n ${i}p $indexFp | awk -v col=$geneIdCol '{print $col}'`
    geneName=`sed -n ${i}p $indexFp | awk -v col=$geneNameCol '{print $col}'`
    anno=$i":"$geneid":"$geneName":"$chr":"$start":"$stop
    echo anno=$anno
    
    sampleFp_raw2=`echo $sampleFp_raw | sed "s|GENE-----|$geneid|g"`
    sampleFp=`echo $sampleFp_raw2 | sed "s|-----|$SGE_TASK_ID_i|g"`
    
    # extract variants in 1Mb flanking region
    outFI=gene.$i
    outFpI=$outdirH2/$outFI
    echo outFpI=$outFpI
    
    
    # no ld score
    # 1= no cov; 2=all cov; 3= cov without peers
    nFlag=0
    passFlag01=NR
    passFlag02=NR
    passFlag03=NR
    passFlag04=NR
    # ld score (for WGS)
    passFlag1=NR
    passFlag2=NR
    passFlag3=NR
    passFlag4=NR
    failFlag=0
    
    geneFlag=`zcat $expFp | grep $geneid | wc -l`
    echo geneFlag=$geneFlag
    
    if [ $geneFlag -gt 0 ];
    then
        nFlag=1
        
        echo -e "cp $sampleFp $outFpI.sample.txt"
        cp $sampleFp $outFpI.sample.txt
        
        # Gerenate phenFp
        #   FID IID phen
        #   NA = missing
        
        awk 'BEGIN{FS="\t";OFS="\t"}{
            if (FNR == 1){
                fI++
            }
            
            if (fI == 1){
                nid=$1
                array[nid]=1
            }else if (fI == 2){
                for (i=5;i<=NF;i++){
                    nid=$i
                    nidArr[i]=nid
                }
            }else{
                for (i=5;i<=NF;i++){
                    nid=nidArr[i]
                    val=$i
                    if (nid in array){
                        print nid,nid,val
                    }
                }
            }
            
            
        }' $outFpI.sample.txt <(zcat $expFp | head -1) <(zcat $expFp | grep $geneid) > $outFpI.gcta.phen.txt
        
        
        nCheck=`cat $outFpI.gcta.phen.txt | wc -l`
        
        if [ $nCheck == 0 ];
        then
            echo empty $outFpI.gcta.phen.txt
            echo "$logSuf $JOB_ID ERROR $SGE_TASK_ID_i $SGE_TASK_ID" > $logdirH/status.ERROR.emptyFile.$SGE_TASK_ID_i.$i
            kill -SIGINT $$
        fi
        
        
        # echo "----->$outFpI.sample.txt"
        # head $outFpI.sample.txt
        
        # echo "----->$expFp"
        # zcat $expFp | head  | cut -f1-10
        
        awk 'BEGIN{FS="\t";OFS="\t"}{
            if(NR == FNR){
                if (FNR == 1){
                    header=$0
                }else{
                    array[$1]=$0
                }
                
            }else{
                if (FNR == 1){
                    print "FID\tIID\tphen\t"header
                }
                if ($1 in array){
                    print $0,array[$1]
                }
            }
        
        }' $peerFp $outFpI.gcta.phen.txt > $outFpI.gcta.phen.plus.txt
        
        # echo "-----> head peerFp"
        # head $peerFp | cut -f1-10
        
        
        # echo "-----> $outFpI.gcta.phen.txt"
        # head $outFpI.gcta.phen.txt
        
        echo peerFp=$peerFp
        
        echo "R"
        
        logSuf=R.res
        rFp_=$rFp_geneExpRes
        inFp_=$outFpI.gcta.phen.plus.txt
        statusFp_=$logdirH2/status.SUCCESS.R.res.$SGE_TASK_ID_i.$i
        outFp_=$outFpI.gcta.phen.res.txt
        covStr_="PEER_1"
        for m in $(seq 2 $nPeer);
        do
            covStr_=$covStr_"+PEER_"$m
        done
        
        if [ ! -f $rFp_ ];
        then
            echo Missing rFp: $rFp_
            echo "$logSuf $JOB_ID ERROR $SGE_TASK_ID_i $SGE_TASK_ID" > $logdirH/status.ERROR.rFpMissing.$SGE_TASK_ID_i.$i
            kill -SIGINT $$
        
        fi
        
        Rscript --verbose ${rFp_} $inFp_ $statusFp_ $outFp_ $covStr_
        
        # Rscript --verbose ${rFp_} $inFp_ $statusFp_ $outFp_ $covStr_  > $logdirH2/out.$SGE_TASK_ID_i.$i.Rout 2>&1
        
        RETVAL=$?
        if [ $RETVAL -ne 0 ];
        then
            echo "$logSuf $JOB_ID ERROR $SGE_TASK_ID_i $SGE_TASK_ID"
            echo "$logSuf $JOB_ID ERROR $SGE_TASK_ID_i $SGE_TASK_ID" > $logdirH/status.ERROR.$logSuf.$SGE_TASK_ID_i.$i
            kill -SIGINT $$
        fi
        
        # QC to make sure chr format matched
        chr_bcf=`bcftools view -H $genoFp | head -1 | cut -f1`
        chr_script=$bcfChrPre$chrStr2
        
        if [ $chr_bcf != $chr_script ];
        then
            echo "$logSuf $JOB_ID ERROR $SGE_TASK_ID_i $SGE_TASK_ID" > $logdirH/status.ERROR.$logSuf.$SGE_TASK_ID_i.$i
            echo "Mismatch chr format"
            kill -SIGINT $$
        fi
        
        
        # Subset genotype file
        
        logSuf=0.bcftools
        bcfFp=$outFpI.bcf
        # bcftools view --min-ac 1 $genoFp chr$chr:$start-$stop  -O b -o $bcfFp
        # run3B.z.4
        echo bcftools view --min-ac 1 --types snps --samples-file $outFpI.sample.txt $genoFp $bcfChrPre$chrStr2:$start-$stop  -O b -o $bcfFp
        
        bcftools view --min-ac 1 --types snps --samples-file $outFpI.sample.txt $genoFp $bcfChrPre$chrStr2:$start-$stop  -O b -o $bcfFp
        
        RETVAL=$?
        if [ $RETVAL -ne 0 ];
        then
            echo "$logSuf $JOB_ID ERROR $SGE_TASK_ID_i $SGE_TASK_ID" > $logdirH/status.ERROR.$logSuf.$SGE_TASK_ID_i.$i
            kill -SIGINT $$
        fi
        bcftools index $bcfFp
        
        nVar=`bcftools view -H $bcfFp | wc -l`
        echo nVar=$nVar
        
        if [ $nVar == 0 ];
        then
            echo -e "$anno\tFAIL:noVar\t$nLn\t$nVar\t$n1\t$n2\t\t\t\t\t\t\t\t\t" >> $outdirH2/out.log.txt
        else
        
            logSuf=1.plink
            rm -rf $logdirH/status.ERROR.$logSuf.$SGE_TASK_ID_i.$i
            plink --bcf $outFpI.bcf --keep-allele-order --allow-extra-chr --make-bed --out $outFpI >/dev/null
            RETVAL=$?
            
            if [ $RETVAL -ne 0 ];
            then
              echo "$logSuf $JOB_ID ERROR $SGE_TASK_ID_i $SGE_TASK_ID" > $logdirH/status.ERROR.$logSuf.$SGE_TASK_ID_i.$i
              kill -SIGINT $$
            fi
            
            
            # phen = res; covariates without peer
                if [ -f $statusFp_ ];
                then
                    gcta64 --grm $outFpI.gcta.allVar --pheno $outFpI.gcta.phen.res.txt $par2 --reml $par3 --out $outFpI.gcta.reml.allVar.res.cov2 >/dev/null
                    
                    RETVAL=$?
                    if [ $RETVAL == 0 ];
                    then
                        meltResults $outFpI.gcta.reml.allVar.res.cov2.hsq "$anno" $outdir/$outFI.converged.allVar.res.cov2.txt
                        cat $outdir/$outFI.converged.allVar.res.cov2.txt >> $outdir/out.converged.allVar.res.cov2.txt
                        passFlag04=1
                    else
                        failFlag=1
                    fi
                else
                    failFlag=1
                fi
            
        # $nVar == 0, else
        fi
        
        rm -rf $outFpI.*
        
        echo -e "\n\n\n\n\n"
        
        if [ $saveCounter == $save ];
        then
            echo from=$outdir
            echo to=$outdirH/
            rsync -ar $outdir $outdirH/
            saveCounter=0
        fi
    fi
done


end_date0=`date +"%Y_%m%d"`
end_time0=`date +"%T"`

echo "FINAL $JOB_ID $start_date0 $start_time0 $host $end_date0 $end_time0 END $SGE_TASK_ID" > $logdirH/status.END.$outfolder

rsync -ar $outdir $outdirH/

# put at end of script
jobsum $JOB_ID | grep usage