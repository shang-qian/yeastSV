#!/bin/bash

######## input SRP ########

#srp=SRP066294
#organism=arabidopsis
srp=$1
organism=$2


reference_genome=../$organism/reference_genome/$organism.fa

####  sraTofastq   ####
find ../$organism/$srp/* | grep sra > ../$organism/$organism.$srp\.fastq.fofn

while read line
do
    ../sratoolkit.2.9.4-centos_linux64/bin/fastq-dump --table SEQUENCE $line -O ../$organism/$srp/ >> ../$organism/dump_result.txt
    echo $line
done < ../$organism/$organism.$srp.fastq.fofn

find ../$organism/$srp/* | grep fastq > ../$organism/$organism.$srp.fofn

####################################################################################################################################################################################
#-----------------------------------------------------------------------Alignment--------------------------------------------------------------------------------------------------#
####################################################################################################################################################################################

bamprefix=../$organism/`echo ${srp/SRP/srp}`
bamprefix=`echo ${bamprefix/ERP/erp}`
bamprefix=`echo ${bamprefix/DRP/drp}`


#######################################################################################################################################################################################
#######################################################################################################################################################################################
#######################################################################################################################################################################################





##################################################################  minimap2 alignment   ###################################################

mmtime=../$organism/mm2time.$srp.txt

echo '###### statistics run time ######'> $mmtime

while read line
do
    srrtmp=`echo $line|cut -d "/" -f7`
    srr=`echo $srrtmp|cut -d "." -f1`
    echo 'alignment start:    ' $srr `    date "+%m_%d %H:%M:%S"`>> $mmtime
    
    minimap2 -t 24 --MD -x map-pb -R "@RG\tID:default\tSM:$srr" $reference_genome $line -a |samtools view -bS > ../$organism/$srp/$srr.mm2.bam
    echo 'alignment end:    ' $srr `    date "+%m_%d %H:%M:%S"`>> $mmtime
     
done < ../$organism/$organism.$srp\.fofn

###-------------------------------------------------------------------------------------------------------------------------------------------------
echo '####### SV calling time #######'>> $mmtime

outputfile=../$organism/$srp.mm2list.fofn
bamlist=`echo ${outputfile/SRP/srp}`
bamlist=`echo ${bamlist/ERP/erp}`
bamlist=`echo ${bamlist/DRP/drp}`
find ../$organism/$srp/* |grep RR | grep mm2 |grep bam |grep -v bai > $bamlist 
samtools merge -f -@ 24 -b $bamlist $bamprefix.merge.mm2.bam
if [ $? -eq 0 ];then
    rm -f ../$organism/$srp/*.mm2.bam
else
    exit
fi

samtools sort -@ 24 -m 2G $bamprefix.merge.mm2.bam -o $bamprefix.merge.sort.mm2.bam
if [ $? -eq 0 ];then
    rm -f $bamprefix.merge.mm2.bam
else
    exit
fi

samtools index -@ 24 $bamprefix.merge.sort.mm2.bam




########################################################################  ngmlr alignment   #######################################################################

ngmlrtime=../$organism/ngmlrtime.$srp.txt

echo '###### statistics run time ######'> $ngmlrtime

while read line
do
    srrtmp=`echo $line|cut -d "/" -f7`
    srr=`echo $srrtmp|cut -d "." -f1`
    echo 'alignment start:    ' $srr `    date "+%m_%d %H:%M:%S"`>> $ngmlrtime
     
    ngmlr -t 24 -r $reference_genome -q $line -o ../$organism/$srp/$srr.ngmlr.bam 
    echo 'alignment end:    ' $srr `    date "+%m_%d %H:%M:%S"`>> $ngmlrtime
     
done < ../$organism/$organism.$srp\.fofn

##-------------------------------------------------------------------------------------------------------------------------------------------------
echo '####### SV calling time #######'>> $ngmlrtime

outputfile=../$organism/$srp.ngmlrlist.fofn
bamlist=`echo ${outputfile/SRP/srp}`
bamlist=`echo ${bamlist/ERP/erp}`
bamlist=`echo ${bamlist/DRP/drp}`
find ../$organism/$srp/* |grep RR | grep ngmlr |grep bam |grep -v bai > $bamlist 
samtools merge -f -@ 24 -b $bamlist $bamprefix.merge.ngmlr.bam
if [ $? -eq 0 ];then
    rm -f ../$organism/$srp/*.ngmlr.bam
else
    exit
fi

samtools sort -@ 24 -m 2G $bamprefix.merge.ngmlr.bam -o $bamprefix.merge.sort.ngmlr.bam
if [ $? -eq 0 ];then
    rm -f $bamprefix.merge.ngmlr.bam
else
    exit
fi

samtools index -@ 24 $bamprefix.merge.sort.ngmlr.bam



#################################################################################################################################################################################
#################################################################################################################################################################################
#################################################################################################################################################################################









####################################################################################################################################################################################
#-----------------------------------------------------------------------SV calling-------------------------------------------------------------------------------------------------#
####################################################################################################################################################################################



###################################################################  sniffles ################################################################################################
 
echo 'sniffles start:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $mmtime

sniffles -t 24 --tmp_file ../tmpdir/temp --genotype --skip_parameter_estimation --min_support 10 -m $bamprefix.merge.sort.mm2.bam -v $bamprefix.mm2.sniffles.vcf

echo 'sniffles end:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $mmtime

echo 'sniffles start:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $ngmlrtime

sniffles -t 24 --tmp_file ../tmpdir/temp --genotype --skip_parameter_estimation --min_support 10 -m $bamprefix.merge.sort.ngmlr.bam -v $bamprefix.ngmlr.sniffles.vcf

echo 'sniffles end:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $ngmlrtime

################################################################  Picky #####################################################################################################
Picky=../Picky/src/picky.pl

echo 'sniffles start:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $mmtime

samtools sort -n $bamprefix.merge.sort.mm2.bam | samtools view -h | $Picky sam2align | $Picky callSV --oprefix $bamprefix.mm2. # --genome $reference_genome

$Picky xls2vcf --xls $bamprefix.mm2.profile.TTLC.xls --xls $bamprefix.mm2.profile.TDSR.xls --xls $bamprefix.mm2.profile.TDC.xls --xls $bamprefix.mm2.profile.INV.xls --xls $bamprefix.mm2.profile.INS.xls --xls $bamprefix.mm2.profile.INDEL.xls --xls $bamprefix.mm2.profile.DEL.xls > $bamprefix.mm2.picky.vcf

rm -f $bamprefix.mm2.profile*

echo 'sniffles end:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $mmtime


echo 'sniffles start:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $ngmlrtime

samtools sort -n $bamprefix.merge.sort.ngmlr.bam | samtools view -h | $Picky sam2align | $Picky callSV --oprefix $bamprefix.ngmlr. # --genome $reference_genome

$Picky xls2vcf --xls $bamprefix.ngmlr.profile.TTLC.xls --xls $bamprefix.ngmlr.profile.TDSR.xls --xls $bamprefix.ngmlr.profile.TDC.xls --xls $bamprefix.ngmlr.profile.INV.xls --xls $bamprefix.ngmlr.profile.INS.xls --xls $bamprefix.ngmlr.profile.INDEL.xls --xls $bamprefix.ngmlr.profile.DEL.xls > $bamprefix.ngmlr.picky.vcf

rm -f $bamprefix.ngmlr.profile*

echo 'sniffles end:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $ngmlrtime


################################################################  smartie-SV ###############################################################################################

echo 'smartie-SV start:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $mmtime

samtools view $bamprefix.merge.sort.mm2.bam | ../smartie-sv/bin/printgaps $reference_genome $bamprefix.mm2.smartiesv
rm $bamprefix.mm2.smartiesv.snp.bed
rm $bamprefix.mm2.smartiesv.indel.bed
rm $bamprefix.mm2.smartiesv.aln.coords.bed

echo 'smartie-SV end:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $mmtime

echo 'smartie-SV start:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $ngmlrtime

samtools view $bamprefix.merge.sort.ngmlr.bam | ../smartie-sv/bin/printgaps $reference_genome $bamprefix.ngmlr.smartiesv
rm $bamprefix.ngmlr.smartiesv.snp.bed
rm $bamprefix.ngmlr.smartiesv.indel.bed
rm $bamprefix.ngmlr.smartiesv.aln.coords.bed

echo 'smartie-SV end:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $ngmlrtime


###############################################################  PBHoney ####################################################################################################

source ../luanmw/resource/PBSuite/setup.sh

echo 'PBHoney start:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $mmtime

Honey.py spots -n 24 $bamprefix.merge.sort.mm2.bam --reference $reference_genome --consensus None -o $bamprefix.mm2.hon

echo 'PBHoney end:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $mmtime

echo 'PBHoney start:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $ngmlrtime

Honey.py spots -n 24 $bamprefix.merge.sort.ngmlr.bam --reference $reference_genome --consensus None -o $bamprefix.ngmlr.hon

echo 'PBHoney end:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $ngmlrtime

deactivate

#################################################################  nanoSV ###################################################################################################

source ../snakemake/bin/activate

echo 'nanoSV start:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $mmtime

NanoSV -o $bamprefix.mm2.nanosv.vcf $bamprefix.merge.sort.mm2.bam -b ../$organism/Sacc.bed  -s samtools

echo 'nanoSV end:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $mmtime

echo 'nanoSV start:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $ngmlrtime

NanoSV -o $bamprefix.ngmlr.nanosv.vcf $bamprefix.merge.sort.ngmlr.bam -b ../$organism/Sacc.bed  -s samtools

echo 'nanoSV end:    ' $srp `    date "+%m_%d %H:%M:%S"`>> $ngmlrtime

deactivate


