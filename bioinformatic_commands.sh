# The following is the commands used for all bioinformatic and phylogenetic analyses in Osborne et al 2020 "Sympatric speciation in Mountain Roses (Metrosideros) on an oceanic island" Phil. Trans. Roy. Soc. B. For analyses that were run on multiple individuals separately, only one example is shown. Simplified file names and paths are used.

# Rcorrector v. 1.0.2
perl run_rcorrector.pl -k 32 -t 20 -od /output/path/ -1 rawreads_1.fq.gz -2 rawreads_2.fq.gz

# Trimmomatic v. 0.36
java -jar trimmomatic-0.36.jar PE -threads 20 -phred33 -trimlog log.txt corrected_reads_1.fq.gz corrected_reads_2.fq.gz trimmed_reads_FP.fq.gz trimmed_reads_FU.fq.gz trimmed_reads_RP.fq.gz trimmed_reads_RU.fq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:75 TOPHRED33

# BinPacker v.1.0. $il is the mean internal fragment size.
for K in 19 25 32 ; do 
BinPacker -s fq -p pair -l trimmed_reads_FP.fq.gz -r trimmed_reads_RP.fq.gz -o /output/path/k_${K} -m FR -g ${il} -k ${K}
done

# Bridger v.2014-12-01. $il is the mean internal fragment size.
for K in 19 25 32 ; do 
	perl Bridger.pl --seqType fq --left trimmed_reads_FP.fq.gz --right trimmed_reads_RP.fq.gz --CPU 20 --pair_gap_length ${il} --kmer_length ${K} --output /output/path/k_${K}
done

# IDBA-tran v.1.1.0. Forward and reverse reads were merged before assembly 
idba_tran -r trimmed_reads_merged.fa --num_threads = 20 --mink 21 --maxk 71 --step 10 -o /output/path/

# Velvet v.1.2.10_b and Oases v.0.2.08. velveth is run first, then velvetg then oases. $il is the mean internal fragment size and $ilsd is the standard deviation of internal insert size.
for K in 21 31 41 51 61 71 ; do 
	velveth /working/dir/path/${K} ${K} -shortPaired -fastq.gz -separate trimmed_reads_FP.fq.gz trimmed_reads_RP.fq.gz
	velvetg /working/dir/path/${K} -ins_length ${il} -ins_length_sd ${ilsd} -unused_reads yes -read_trkg yes -amos_file yes
	oases /working/dir/path/${K} -ins_length ${il} -ins_length_sd ${ilsd} -unused_reads yes
done

# Shannon v.0.0.2
for K in 21 31 41 51 61 71 ; do 
	python shannon.py -o /output/path/K_${K} --left trimmed_reads_FP.cor.fa --right trimmed_reads_RP.cor.fa -p 20 -K ${K} 
done

# SOAPdenovo-trans. Contents of the config file are pasted below the command line. $il is the mean internal fragment size.
for K in 21 31 41 51 61 71 ; do 
	SOAPdenovo-Trans-127mer all -s /path/to/config/file.txt -p 20 -K ${K} -o /output/path/k_${K}
done
## CONFIG FILE:
#max_rd_len=100
#[LIB]
#avg_ins=$il
#asm_flags=3
#reverse_seq=0
#q1=trimmed_reads_FP.cor.fq.gz
#q2=trimmed_reads_RP.cor.fq.gz

# TransABYSS v.1.5.5
for K in 21 31 41 51 61 71 ; do 
	transabyss --name k_${K} --outdir /output/path/k_${K} --kmer ${K} --threads 20 --pe trimmed_reads_FP.cor.fq.gz trimmed_reads_RP.cor.fq.gz
done

# Trinity v. 2.4.0
for K in 19 25 32 ; do 
Trinity --seqType fq --max_memory 92G --left trimmed_reads_FP.fq.gz --right trimmed_reads_RP.fq.gz --CPU 20 --KMER_SIZE ${K} --output /output/path/k_${K}
done

# TransDecoder v.3.0.1
TransDecoder.LongOrfs -t contigs.fasta
TransDecoder.Predict -t contigs.fasta --cpu 50 --single_best_orf 

# CD-HIT-EST v.4.6.1 
cd-hit-est -T 0 -i contigs.cds.fasta -c 0.99 -o cdhit.output.fasta -M 0 -aL 0.005 -aS 1 -G 0 -d 0 


# STAR v.2.5.3a  
STAR --runThreadN 50 --runMode genomeGenerate --genomeDir /path/reference/ --genomeFastaFiles cdhit.output.fasta --genomeChrBinNbits 13 
STAR --runThreadN 20 --limitBAMsortRAM 100000000000 --genomeDir /path/reference/ --readFilesIn trimmed_reads_FP.fq.gz trimmed_reads_RP.fq.gz  --outFilterMultimapNmax 10000000 --readFilesCommand zcat --outSAMtype BAM Sorted  --outFileNamePrefix /output/path/mapping.bam

# Corset v.1.06
corset -d 0.3 -m 0 -g MN,MN,MN,MN,MN,MN,MN,MN,MN,MN,MS,MS,MS,MS,MS,MS,MS,MS,MS,MS MN01_mappings.bam MN02_mappings.bam MN03_mappings.bam MN04_mappings.bam MN05_mappings.bam MN06_mappings.bam MN07_mappings.bam MN08_mappings.bam MN09_mappings.bam MN10_mappings.bam MS01_mappings.bam MS02_mappings.bam MS03_mappings.bam MS04_mappings.bam MS05_mappings.bam MS06_mappings.bam MS07_mappings.bam MS08_mappings.bam MS09_mappings.bam MS10_mappings.bam

# Picard tools v.2.6.0  
picard AddOrReplaceReadGroups I=mappings.bam O=mappings_rg.bam RGLB=library RGPL=platform RGPU=machine RGSM=sample
picard MarkDuplicates I=mappings_rg.bam O=mappings_rg_rd.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics

# GATK SplitNCigarReads v. 3.7-0-gcfedb67
java -jar GenomeAnalysisTK.jar  -T SplitNCigarReads -R cdhit.output.fasta -I mappings_rg_rd.bam -o mappings_rg_rd_s.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

# samtools-bcftools v. 1.3.1 
samtools mpileup -Q 20 -t DP -d 8000 -vuf cdhit.output.fasta mappings_rg_rd_s.bam | bcftools call -M -f GQ -mg 3 -Ou | bcftools filter -g 3 -e '(GT="0/1" & ((DP4[0]+DP4[1]) < 2 || (DP4[2]+DP4[3]) < 2)) || (DP4[0]+DP4[1]+DP4[2]+DP4[3]) < 3 || GQ < 20 || %QUAL < 20' > snps.vcf

# vcf2fas (available from https://github.com/brunonevado/vcf2fas). vcflist is a text file with vcfs from each individual on a separate line.
vcf2fas_linux -reference cdhit.output.fasta -vcfs vcf_list.txt

# BLASTx v2.2.25 
blastx -query cdhit.output.fasta -db blastdb -out blastout.csv -outfmt "10 qseqid sseqid sgi sacc saccver qlen slen length qstart qend sstart send pident nident evalue bitscore" -max_target_seqs 1 -max_hsps 1 -num_threads 32 

# RAxML v. 8.1.17 
/raxmlHPC-PTHREADS-AVX -s concat_alignments.fasta -T 32 -f a -x 12345 -# 100 -p 12345 -m GTRGAMMA -n output_tree

# PLINK v. 1.9 
plink --bfile SNPs --read-genome plink.out.genome --allow-extra-chr --cluster --mds-plot 2 


# fastSTRUCTURE  v.1.0
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ; do
	structure.py -K ${K} --cv=5 --input=SNPs --output FS_output1 --full --seed=12345 --prior=simple
done
# following identification of 2 as most likely number of clusters:
structure.py -K 2 --cv=5 --input=SNPs --output=FS_output2 --full --seed=12345 --prior=logistic









