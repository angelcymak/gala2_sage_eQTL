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

gtfFp=${gtfFp}
tissue=${tissue}
# eGeneFp=${eGeneFp}
batchSize=${batchSize}
dbsnpdir=${dbsnpdir}
dbsnpFI=${dbsnpFI}
covFp_raw=${covFp}

vcfFp_raw=${vcfFp}
expFp_raw=${expFp}
rundir=${rundir}
geneTypeColName=${geneTypeColName}


if [ 1 == 0 ];
then
    e_=MX
    # Prepare gene annotation file
    chr=22
    tmpdir=/wynton/group/burchard/maka/project/topmedp4.rnaseq/gtex.preprocessing
    fastqtlTag=all                
    covFp=$tmpdir/run1.${e_}.cov.str.txt.gz
    vcfFp_raw=$tmpdir/bcftools_isec.DP0.subset.mafp01.genop05.${e_}/isec.chr-----/0000.vcf.gz
    
    geneAnnoFullFp=/wynton/group/burchard/maka/data_shared/GRCh38/gencode.v30.GRCh38.annotation.ERCC.gOnly.bed
    tissue=whole_blood
    outdir=/wynton/group/burchard/maka/project/topmedp4.rnaseq/predictdb/run1.predictdb_dataprep.all.MX
    expFp_raw=/wynton/group/burchard/maka/project/topmedp4.rnaseq/fastqtl/all.docker.fastqtl.MX/out.$chr.expression.bed.gz
    # eGeneFp=/wynton/home/burchard/maka/project/topmed.rnaseqb1-2/fastqtl/docker.fastqtl.MX/permute/allChr.eGenes.txt.gz
    batchSize=200
    dbsnpdir=/wynton/group/burchard/maka/data_shared/GRCh38
    dbsnpFI=All_20170710
fi

if [ ! -z $indexFp ];
then
    SGE_TASK_ID_i=`sed -n ${SGE_TASK_ID}p $indexFp`
else
    SGE_TASK_ID_i=$SGE_TASK_ID
fi
echo "SGE_TASK_ID=$SGE_TASK_ID SGE_TASK_ID_i=$SGE_TASK_ID_i"
chr="${chrList[$SGE_TASK_ID_i]}"

covFp_raw2=`echo $covFp_raw | sed "s|chr-----|$chr|g"`
covFp=`echo $covFp_raw2 | sed "s|-----|$SGE_TASK_ID_i|g"`

vcfFp_raw2=`echo $vcfFp_raw | sed "s|chr-----|$chr|g"`
vcfFp=`echo $vcfFp_raw2 | sed "s|-----|$SGE_TASK_ID_i|g"`

expFp_raw2=`echo $expFp_raw | sed "s|chr-----|$chr|g"`
expFp=`echo $expFp_raw2 | sed "s|-----|$SGE_TASK_ID_i|g"`

qcFpStr_raw=${qcFpStr}
qcFpStr_raw2=`echo $qcFpStr_raw | sed "s|chr-----|$chr|g"`
qcFpStr=`echo $qcFpStr_raw2 | sed "s|-----|$SGE_TASK_ID_i|g"`
qcFpStr=$qcFpStr:$inFp:$covFp
qc_tag $qcFpStr $logdir $SGE_TASK_ID_i

RETVAL=$?
if [ $RETVAL -ne 0 ];
then
  kill -SIGINT $$
fi
rm -rf $logdir/status.HALT.qc.$SGE_TASK_ID_i

outdir=$(dirname "$logdir")
if [ ! -d $logdir ];then mkdir -p $logdir; fi
if [ ! -d $rundir/prepare_data ]; then mkdir -p $rundir/prepare_data; fi


start_date0=`date +"%Y_%m%d"`
start_time0=`date +"%T"`
logSuf=run

echo -n "$logSuf $JOB_ID $start_date0 $start_time0 $host $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.RUN.$logSuf.$SGE_TASK_ID_i
rm -rf $logdir/status.END.$logSuf.$SGE_TASK_ID_i


chrChar=$chr
chr=$SGE_TASK_ID

#   it can contain all gene because what got analyzed was determined by the gene annotation file

if [ 1 == 0  ];
then
    # use eGene only
    zcat $expFp | head -1 > $outdir/tmp.$SGE_TASK_ID_i.txt
    awk 'BEGIN{FS="\t";OFS="\t"}{
        if (NR == FNR){
            gene=$1
            array[gene]=1
        }else{
            gene=$4
            if (gene in array){
                print
            }
        }
    }'  <(zcat $eGeneFp) <(zcat $expFp | tail -n +2) >> $outdir/tmp.$SGE_TASK_ID_i.txt
    
    n=`cat $outdir/tmp.$SGE_TASK_ID_i.txt | wc -l`
    echo $chr $n
    if [ $n -gt 1 ];
    then
        awk '{
            if (FNR == 1){
                $4="NAME"
            }
            printf $4
            for (i=5;i<=NF;i++){
                printf "\t"$i
            }
            printf "\n"
        }' $outdir/tmp.$SGE_TASK_ID_i.txt > $rundir/prepare_data/${tissue}_Analysis.expression.$chr.txt
        rm -rf $outdir/tmp.$SGE_TASK_ID_i.txt
    fi
fi

awk '{
    if (FNR == 1){
        $4="NAME"
    }
    printf $4
    for (i=5;i<=NF;i++){
        printf "\t"$i
    }
    printf "\n"
}' <(zcat $expFp) > $rundir/prepare_data/${tissue}_Analysis.expression.$chr.txt


start_date=`date +"%Y_%m%d"`
start_time=`date +"%T"`
logSuf=exp.final

echo -n "$logSuf $JOB_ID $start_date $start_time $host $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.RUN.$logSuf.$SGE_TASK_ID_i
rm -rf $logdir/status.ERROR.$logSuf.$SGE_TASK_ID_i
rm -rf $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID_i

# gene annotation file
# gene annotation file determined what gets analysed
# convert form this
# chrf2  start   end      geneName        geneId  txid    feature geneType
# chr gene_id           gene_name     start end    gene_type
geneAnnoFp=$rundir/prepare_data/gene.annotation.parsed.$chr.txt

geneTypeCol=`head -1 $geneAnnoFullFp | awk -v colName=$geneTypeColName 'BEGIN{FS="\t"}{
  for (i=1;i<=NF;i++){
    if ($i == colName){
      print i
    }
  }
}'`
        
echo -e "chr\tgene_id\tgene_name\tstart\tend\tgene_type" > $geneAnnoFp
tail -n +2 $geneAnnoFullFp | sed s/^chr$chrChar/$chr/g | awk -v geneTypeCol=$geneTypeCol 'BEGIN{FS="\t";OFS="\t"}{
    if (NR == FNR){
        gene=$1
        array[gene]=1
    }else{
        chr=$1
        start=$2
        end=$3
        geneName=$4
        gene=$5
        # txid=$6
        # feature=$7
        geneType=$geneTypeCol
        if (gene in array){
            print chr,gene,geneName,start,end,geneType
        }
    }
}' $rundir/prepare_data/${tissue}_Analysis.expression.$chr.txt - >> $geneAnnoFp




n1=`cat $rundir/prepare_data/${tissue}_Analysis.expression.$chr.txt | wc -l `
n2=`cat $geneAnnoFp | wc -l `
if [ $n1 != $n2 ];
then
    echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time ERROR $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.ERROR.$logSuf.$SGE_TASK_ID_i
    kill -SIGINT $$
else

    echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time SUCCESS $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID_i
fi


# Further split gene annotation file (100 gene per file)
nGene=`grep -ve "^#" $geneAnnoFp | wc -l`
nBatch=`echo $nGene $batchSize | awk '{nGene=$1; batchSize=$2; nBatch=int(nGene/batchSize); if (nGene % batchSize > 0){nBatch=nBatch+1}; print nBatch}'`
maxBatchI=`echo $nBatch | awk '{print $1-1}'`
echo -n > $rundir/prepare_data/list.$chr.txt
for batchI in $(seq 0 $maxBatchI);
do
    geneAnnoFpTmp=$rundir/prepare_data/gene.annotation.parsed.$chr.$batchI.txt
    head -1 $geneAnnoFp > $geneAnnoFpTmp
    
    startLine=`echo $batchI $batchSize | awk '{batchI=$1; batchSize=$2; print batchI*batchSize+1}'`
    startLine2=`echo $batchI $batchSize | awk '{batchI=$1; batchSize=$2; print batchI*batchSize+2}'`
    endLine=`echo $batchI $batchSize | awk '{batchI=$1; batchSize=$2; print (batchI+1)*batchSize}'`
    if [ $endLine -gt $nGene ];
    then
        endLine=$nGene
    fi
    echo -e "$batchI\t$startLine\t$endLine\t$chr\t$geneAnnoFpTmp" >> $rundir/prepare_data/list.$chr.txt
    echo -e "$batchI\t$startLine\t$endLine\t$chr\t$nGene" 
    
    # This version has a known bug that the first gene of $geneAnnoFp will be excluded
    # tail -n +2 $geneAnnoFp | tail -n +$startLine2 | head -$batchSize >> $geneAnnoFpTmp
    
    # This is the correct version
    cat $geneAnnoFp | tail -n +$startLine | head -$batchSize | grep -ve "gene_name" >> $geneAnnoFpTmp
done



# prepare covariate file (no missingness allowed)
#   needed for each chr because of different local ancestry by gene
logSuf=cov
echo -n "$logSuf $JOB_ID NA NA $host $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.RUN.$logSuf.$SGE_TASK_ID_i
rm -rf $logdir/status.ERROR.$logSuf.$SGE_TASK_ID_i
rm -rf $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID_i

if [ 1 == 0 ];
then
    covFp=/wynton/group/burchard/maka/project/topmedp4.rnaseq/gtex.preprocessing/gtex.random/rand1.AA.cov.str.txt.gz
    rundir=/wynton/group/burchard/maka/project/topmedp4.rnaseq/predictdb/run2.gtex670/gtexCaCt.AA
    tissue=test
    chr=22
fi



# must not have NA in covariates
zcat $covFp > $rundir/prepare_data/${tissue}_Analysis.combined_covariates.$chr.txt

len1=`zcat $covFp | wc -l`

len2=`cat $rundir/prepare_data/${tissue}_Analysis.combined_covariates.$chr.txt | wc -l `


if [ $len1 != $len2 ];
then
  echo "$logSuf $JOB_ID ERROR NA NA NA NA NA $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.ERROR.$logSuf.$SGE_TASK_ID_i
  kill -SIGINT $$
else
  echo "$logSuf $JOB_ID SUCCESS NA NA NA NA NA $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID_i
fi



# Process dbsnp data (see par_topmedp4_rnaseq.sh) run once only


# Genotype data (independent of eQTLs analysis, include all snps)
# create a snp only verson with genotype

start_date=`date +"%Y_%m%d"`
start_time=`date +"%T"`
logSuf=bcftools.1

echo -n "$logSuf $JOB_ID $start_date $start_time $host $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.RUN.$logSuf.$SGE_TASK_ID_i
rm -rf $logdir/status.ERROR.$logSuf.$SGE_TASK_ID_i
rm -rf $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID_i


bcftools view --types snps -O b -o $outdir/tmp1.$SGE_TASK_ID_i.bcf $vcfFp
RETVAL=$?

end_date=`date +"%Y_%m%d"`
end_time=`date +"%T"`

if [ $RETVAL -ne 0 ];
then
  echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time ERROR $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.ERROR.$logSuf.$SGE_TASK_ID_i
  kill -SIGINT $$
else
  echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time SUCCESS $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID_i
fi

bcftools index $outdir/tmp1.$SGE_TASK_ID_i.bcf




start_date=`date +"%Y_%m%d"`
start_time=`date +"%T"`
logSuf=bcftools.2

echo -n "$logSuf $JOB_ID $start_date $start_time $host $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.RUN.$logSuf.$SGE_TASK_ID_i
rm -rf $logdir/status.ERROR.$logSuf.$SGE_TASK_ID_i
rm -rf $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID_i

# intersect geno-snpOnly and dbsnp-snpOnly and write intersect in geno-snpOnly (-w1: output intersection in file 1)
bcftools isec -n =2 -p $outdir/tmp.$SGE_TASK_ID_i $outdir/tmp1.$SGE_TASK_ID_i.bcf $dbsnpdir/$dbsnpFI.$chrChar.2.bcf

RETVAL=$?

end_date=`date +"%Y_%m%d"`
end_time=`date +"%T"`

if [ $RETVAL -ne 0 ];
then
  echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time ERROR $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.ERROR.$logSuf.$SGE_TASK_ID_i
  kill -SIGINT $$
else
  echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time SUCCESS $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID_i
fi


# QC check
# bcftools view -H -G $outdir/tmp.$SGE_TASK_ID_i/0000.vcf | wc -l
# bcftools view -H -G $outdir/tmp.$SGE_TASK_ID_i/0001.vcf | wc -l

cat $outdir/tmp.$SGE_TASK_ID_i/0000.vcf | grep "#CHROM" | cut -f10- > $outdir/tmp1.$SGE_TASK_ID_i.sample.txt


start_date=`date +"%Y_%m%d"`
start_time=`date +"%T"`
logSuf=bcftools.3

echo -n "$logSuf $JOB_ID $start_date $start_time $host $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.RUN.$logSuf.$SGE_TASK_ID_i
rm -rf $logdir/status.ERROR.$logSuf.$SGE_TASK_ID_i
rm -rf $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID_i

bcftools +dosage $outdir/tmp.$SGE_TASK_ID_i/0000.vcf  > $outdir/tmp1.$SGE_TASK_ID_i.dosage.v1.txt

RETVAL=$?

end_date=`date +"%Y_%m%d"`
end_time=`date +"%T"`

if [ $RETVAL -ne 0 ];
then
  echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time ERROR $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.ERROR.$logSuf.$SGE_TASK_ID_i
  kill -SIGINT $$
else
  echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time SUCCESS $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID_i
fi


cat $outdir/tmp1.$SGE_TASK_ID_i.dosage.v1.txt | awk '{print NF}' | head

cat $outdir/tmp1.$SGE_TASK_ID_i.sample.txt > $outdir/tmp1.$SGE_TASK_ID_i.dosage.v2.txt
tail -n +2 $outdir/tmp1.$SGE_TASK_ID_i.dosage.v1.txt |  sed 's/\.0//g' | cut -f5- >> $outdir/tmp1.$SGE_TASK_ID_i.dosage.v2.txt
# head $outdir/tmp1.$SGE_TASK_ID_i.dosage.v2.txt | cut -f1-10
cat $outdir/tmp1.$SGE_TASK_ID_i.dosage.v2.txt | awk '{print NF}' | uniq -c


start_date=`date +"%Y_%m%d"`
start_time=`date +"%T"`
logSuf=bcftools.final

echo -n "$logSuf $JOB_ID $start_date $start_time $host $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.RUN.$logSuf.$SGE_TASK_ID_i
rm -rf $logdir/status.ERROR.$logSuf.$SGE_TASK_ID_i
rm -rf $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID_i


echo "varID" > $outdir/tmp1.$SGE_TASK_ID_i.varID.txt
bcftools query -f '%ID\n' $outdir/tmp.$SGE_TASK_ID_i/0000.vcf >> $outdir/tmp1.$SGE_TASK_ID_i.varID.txt
n1=`cat $outdir/tmp1.$SGE_TASK_ID_i.dosage.v2.txt | wc -l`
n2=`cat $outdir/tmp1.$SGE_TASK_ID_i.varID.txt | wc -l`
if [ $n1 != $n2 ];
then
    echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time ERROR $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.ERROR.$logSuf.$SGE_TASK_ID_i
    kill -SIGINT $$
else
    echo "$logSuf $JOB_ID $start_date $start_time $host $end_date $end_time SUCCESS $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.SUCCESS.$logSuf.$SGE_TASK_ID_i
    paste $outdir/tmp1.$SGE_TASK_ID_i.varID.txt $outdir/tmp1.$SGE_TASK_ID_i.dosage.v2.txt > $rundir/prepare_data/genotype.snp.$chr.txt
fi


# e.g. pos 1077948
# snp annotation file (I didn't see the rsID was used in a useful way in the script gtex_v7_nested_cv_elnet.R, may be rsID was not important as long as it's not empty; Another information that is not used was R2)
logSuf=FINAL

echo -e "chromosome\tpos\tvarID\tref_vcf\talt_vcf\trsid_dbSNP150" > $rundir/prepare_data/genotype.snp.anno.$chr.txt
awk 'BEGIN{FS="\t";OFS="\t"}{
    if (NR == FNR){
        chr=$1
        pos=$2
        id=$3
        ref=$4
        alt=$5
        key=chr"\t"pos"\t"ref"\t"alt
        array[key]=id
    }else{
        chr=$1
        pos=$2
        id=$3
        ref=$4
        alt=$5
        key=chr"\t"pos"\t"ref"\t"alt
        if (key in array){
            rsid=array[key]
            print chr,pos,id,ref,alt,rsid
        }else{
            exit 1
        }
        
    }
}' <(bcftools view -H $outdir/tmp.$SGE_TASK_ID_i/0001.vcf) <(bcftools view -H $outdir/tmp.$SGE_TASK_ID_i/0000.vcf) | sed s/^chr$chrChar/$chr/g >> $rundir/prepare_data/genotype.snp.anno.$chr.txt



echo "$logSuf $JOB_ID $start_date0 $start_time0 $host $end_date $end_time END $SGE_TASK_ID_i $SGE_TASK_ID" > $logdir/status.END.$logSuf.$SGE_TASK_ID_i

rm -rf $outdir/tmp1.$SGE_TASK_ID_i.*
rm -rf $outdir/tmp.$SGE_TASK_ID_i

# put at end of script
jobsum $JOB_ID | grep usage