#/usr/bin/bash

########################################################################################################################################################################################################


#sample="SAMD00082707"
sample=$1
##############################   minimap2    ################

##############  Sniffles   #########################

sort -k1,1V -k2,2n -o $sample.mm2.sniffles.sort.bed $sample.mm2.sniffles.bed 
bedtools getfasta -fi ../yeast.fa -bed $sample.mm2.sniffles.sort.bed -fo $sample.mm2.sniffles.fa 

minimap2 -t 24 --MD -x map-pb -R "@RG\tID:default\tSM:$sample" $sample.mm2.sniffles.fa ../DRP003723/DRR095880.fastq -a | samtools view -@ 24 -bS > $sample.mm2.sniffles.bam

samtools sort -@ 24 -o $sample.mm2.sniffles.sort.bam $sample.mm2.sniffles.bam
bedtools genomecov -ibam $sample.mm2.sniffles.sort.bam -bga > $sample.mm2.sniffles.stat
##############  NanoSV   #########################
sort -k1,1V -k2,2n -o $sample.mm2.nanosv.sort.bed $sample.mm2.nanosv.bed 
bedtools getfasta -fi ../yeast.fa -bed $sample.mm2.nanosv.sort.bed -fo $sample.mm2.nanosv.fa 

minimap2 -t 24 --MD -x map-pb -R "@RG\tID:default\tSM:$sample" $sample.mm2.nanosv.fa ../DRP003723/DRR095880.fastq -a |samtools view -@ 24 -bS > $sample.mm2.nanosv.bam

samtools sort -@ 24 -o $sample.mm2.nanosv.sort.bam $sample.mm2.nanosv.bam
bedtools genomecov -ibam $sample.mm2.nanosv.sort.bam -bga > $sample.mm2.nanosv.stat

##############  PBHoney   #########################
sort -k1,1V -k2,2n -o $sample.mm2.hon.sort.bed $sample.mm2.hon.bed 
bedtools getfasta -fi ../yeast.fa -bed $sample.mm2.hon.sort.bed -fo $sample.mm2.hon.fa 

minimap2 -t 24 --MD -x map-pb -R "@RG\tID:default\tSM:$sample" $sample.mm2.hon.fa ../DRP003723/DRR095880.fastq -a | samtools view -@ 24 -bS > $sample.mm2.hon.bam

samtools sort -@ 24 -o $sample.mm2.hon.sort.bam $sample.mm2.hon.bam
bedtools genomecov -ibam $sample.mm2.hon.sort.bam -bga > $sample.mm2.hon.stat

##############  smartie-sv   #########################
sort -k1,1V -k2,2n -o $sample.mm2.smartiesv.sort.bed $sample.mm2.smartiesv.bed 
bedtools getfasta -fi ../yeast.fa -bed $sample.mm2.smartiesv.sort.bed -fo $sample.mm2.smartiesv.fa 

minimap2 -t 24 --MD -x map-pb -R "@RG\tID:default\tSM:$sample" $sample.mm2.smartiesv.fa ../DRP003723/DRR095880.fastq -a | samtools view -@ 24 -bS >$sample.mm2.smartiesv.bam

samtools sort -@ 24 -o $sample.mm2.smartiesv.sort.bam $sample.mm2.smartiesv.bam
bedtools genomecov -ibam $sample.mm2.smartiesv.sort.bam -bga > $sample.mm2.smartiesv.stat

##############  Pickey   #########################
sort -k1,1V -k2,2n -o $sample.mm2.picky.sort.bed $sample.mm2.picky.bed 
bedtools getfasta -fi ../yeast.fa -bed $sample.mm2.picky.sort.bed -fo $sample.mm2.picky.fa 

minimap2 -t 24 --MD -x map-pb -R "@RG\tID:default\tSM:$sample" $sample.mm2.picky.fa ../DRP003723/DRR095880.fastq -a | samtools view -@ 24 -bS > $sample.mm2.picky.bam

samtools sort -@ 24 -o $sample.mm2.picky.sort.bam $sample.mm2.picky.bam
bedtools genomecov -ibam $sample.mm2.picky.sort.bam -bga > $sample.mm2.picky.stat


##############################   NGMLR    ################

##############  Sniffles   #########################

sort -k1,1V -k2,2n -o $sample.ngmlr.sniffles.sort.bed $sample.ngmlr.sniffles.bed 
bedtools getfasta -fi ../yeast.fa -bed $sample.ngmlr.sniffles.sort.bed -fo $sample.ngmlr.sniffles.fa 

minimap2 -t 24 --MD -x map-pb -R "@RG\tID:default\tSM:$sample" $sample.ngmlr.sniffles.fa ../DRP003723/DRR095880.fastq -a | samtools view -@ 24 -bS > $sample.ngmlr.sniffles.bam

samtools sort -@ 24 -o $sample.ngmlr.sniffles.sort.bam $sample.ngmlr.sniffles.bam
bedtools genomecov -ibam $sample.ngmlr.sniffles.sort.bam -bga > $sample.ngmlr.sniffles.stat



##############  NanoSV   #########################
sort -k1,1V -k2,2n -o $sample.ngmlr.nanosv.sort.bed $sample.ngmlr.nanosv.bed 
bedtools getfasta -fi ../yeast.fa -bed $sample.ngmlr.nanosv.sort.bed -fo $sample.ngmlr.nanosv.fa 

minimap2 -t 24 --MD -x map-pb -R "@RG\tID:default\tSM:$sample" $sample.ngmlr.nanosv.fa ../DRP003723/DRR095880.fastq -a |samtools view -@ 24 -bS > $sample.ngmlr.nanosv.bam

samtools sort -@ 24 -o $sample.ngmlr.nanosv.sort.bam $sample.ngmlr.nanosv.bam
bedtools genomecov -ibam $sample.ngmlr.nanosv.sort.bam -bga > $sample.ngmlr.nanosv.stat

##############  PBHoney   #########################
sort -k1,1V -k2,2n -o $sample.ngmlr.hon.sort.bed $sample.ngmlr.hon.bed 
bedtools getfasta -fi ../yeast.fa -bed $sample.ngmlr.hon.sort.bed -fo $sample.ngmlr.hon.fa 

minimap2 -t 24 --MD -x map-pb -R "@RG\tID:default\tSM:$sample" $sample.ngmlr.hon.fa ../DRP003723/DRR095880.fastq -a | samtools view -@ 24 -bS > $sample.ngmlr.hon.bam

samtools sort -@ 24 -o $sample.ngmlr.hon.sort.bam $sample.ngmlr.hon.bam
bedtools genomecov -ibam $sample.ngmlr.hon.sort.bam -bga > $sample.ngmlr.hon.stat

##############  smartie-sv   #########################
sort -k1,1V -k2,2n -o $sample.ngmlr.smartiesv.sort.bed $sample.ngmlr.smartiesv.bed 
bedtools getfasta -fi ../yeast.fa -bed $sample.ngmlr.smartiesv.sort.bed -fo $sample.ngmlr.smartiesv.fa 

minimap2 -t 24 --MD -x map-pb -R "@RG\tID:default\tSM:$sample" $sample.ngmlr.smartiesv.fa ../DRP003723/DRR095880.fastq -a | samtools view -@ 24 -bS > $sample.ngmlr.smartiesv.bam

samtools sort -@ 24 -o $sample.ngmlr.smartiesv.sort.bam $sample.ngmlr.smartiesv.bam
bedtools genomecov -ibam $sample.ngmlr.smartiesv.sort.bam -bga > $sample.ngmlr.smartiesv.stat

##############  Pickey   #########################
sort -k1,1V -k2,2n -o $sample.ngmlr.picky.sort.bed $sample.ngmlr.picky.bed 
bedtools getfasta -fi ../yeast.fa -bed $sample.ngmlr.picky.sort.bed -fo $sample.ngmlr.picky.fa 

minimap2 -t 24 --MD -x map-pb -R "@RG\tID:default\tSM:$sample" $sample.ngmlr.picky.fa ../DRP003723/DRR095880.fastq -a | samtools view -@ 24 -bS > $sample.ngmlr.picky.bam

samtools sort -@ 24 -o $sample.ngmlr.picky.sort.bam $sample.ngmlr.picky.bam
bedtools genomecov -ibam $sample.ngmlr.picky.sort.bam -bga > $sample.ngmlr.picky.stat
