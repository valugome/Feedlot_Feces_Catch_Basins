# BIOINFORMATIC PIPELINE 
## Trimmomatic 
Trimmomatic (version 0.39-Java-11) was run within the AMR++ v3.0 pipeline 

Installation followed the tutorial: https://github.com/Microbial-Ecology-Group/AMRplusplus/blob/master/docs/installation.md

The ‘trim_qc’ subworkflow from the ‘main_AMR++.nf’ pipeline was run:
```
nextflow run main_AMR++.nf --pipeline trim_qc 
```
Only paired end reads for which both mates passed quality control were included in downstream analyses.


## SortMeRNA
The SortMeRNA database was downloaded from https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz

Input reads were quality-controlled (QC) paired mates (i.e., paired-end reads output by Trimmomatic for which both mates passed QC)

```
sortmerna \
	--ref $DB \
	--reads ${qc_trimmed_reads}/${sample}*.1P.fastq.gz \
	--reads ${qc_trimmed_reads}/${sample}*.2P.fastq.gz \
	--fastx \
	--blast '1 cigar qcov' \
	--other ${sample_outdir}/no_rRNA \
	--aligned ${sample_outdir}/rRNA \
	--paired_in \
	--out2 \
	--log 
```

## FLASH Read Merging
Overlapping quality-controlled (QC) metagenomic and rRNA-filtered metatranscriptomic paired-end reads were merged using FLASH (version 2.2.0).

```
flash \
	${trimmed_reads_dir}/${readID}.1P.fastq.gz \
	${trimmed_reads_dir}/${readID}.2P.fastq.gz  \
	-M 150 -d $output_dir  -o $sampleID \
	--interleaved-output \
	--compress-prog pigz \
	--suffix=gz \
	--compress-prog-args '-p 8'
```

## Host Removal
Merged and unmerged reads were aligned separately to the UMD_3.1.1 Bos taurus reference genome (NCBI RefSeq accession number: GCF_000003055.6) using BWA-MEM (version 0.7.18).

For merged reads: 
```
###Align to host genome
bwa mem ${index_file} ${merged_reads_dir}/${sample_id}.extendedFrags.fastq.gz -t ${task_cpus} \
	| samtools sort -@ ${task_cpus} -o ${sample_id}.merged.host.sorted.bam

###Index BAM
samtools index ${sample_id}.merged.host.sorted.bam

##Get idxstats
samtools idxstats ${sample_id}.merged.host.sorted.bam > ${sample_id}.merged.samtools.idxstats

# Extract unmapped reads (non-host) and compress
samtools view -b -f 4 -F 256 ${sample_id}.merged.host.sorted.bam \
	| samtools fastq -@ ${task_cpus} -c 6 - \
	| pigz -p ${task_cpus} -c > ${sample_id}.merged.non.host.fastq.gz
```

For unmerged reads:
```
# Align to host genome (interleaved paired-end)
bwa mem -p ${index_file} ${unmerged_reads_dir}/${sample_id}.notCombined.fastq.gz -t ${task_cpus} \
	| samtools sort -@ ${task_cpus} -o ${sample_id}.unmerged.host.sorted.bam

# Index BAM
samtools index ${sample_id}.unmerged.host.sorted.bam

# Get idxstats
samtools idxstats ${sample_id}.unmerged.host.sorted.bam > ${sample_id}.unmerged.samtools.idxstats

# Collate -- grouping reads by readID.
tmpdir=$(mktemp -d -t ${sample_id}_collate_XXXX) # Make temp directory for collate to avoid collisions
samtools collate -@ ${task_cpus} -T ${tmpdir} ${sample_id}.unmerged.host.sorted.bam -o ${sample_id}.unmerged.host.collated.bam
#rm -rf ${tmpdir}

# Extract unmapped reads (non-host) and compress. samtools fastq -n makes sure the reads are kept interleaved. -f 12 keeps reads for which both mates are unmapped.
samtools view -b -f 12 ${sample_id}.unmerged.host.collated.bam \
	| samtools fastq -@ ${task_cpus} -c 6 -n - \
	| pigz -p ${task_cpus} -c > ${sample_id}.unmerged.non.host.fastq.gz

#Since the version of Kraken2 run did not accept interleaved reads, had to go back and extract forward and reverse unmerged reads in separate files: 
samtools view -h -f 12 -b ${sample_id}.unmerged.host.sorted.bam | \
	samtools sort -n -@ ${task_cpus} - | \
	samtools fastq -@ ${task_cpus} -0 /dev/null -s /dev/null -n \
	-1 ${sample_id}.unmerged.non.host.R1.fastq.gz \
	-2 ${sample_id}.unmerged.non.host.R2.fastq.gz
```

## Kraken2 Taxonomic Classification
Non-host reads were classified taxonomically using Kraken2 (version 2.1.2) with a confidence threshold of 0 and a prebuilt ‘core_nt’ Kraken2 index (comprised of NCBI’s Core Nucleotide database – index built 10/15/2025), which was downloaded from https://benlangmead.github.io/aws-indexes/k2. 

Merged reads were classified in single-end mode:
```
kraken2 $reads \
	--db ${kraken2_db} \
	--confidence 0.0 \
	--report ${report_dir}/${sample_id}.merged.kraken.report \
	--output ${output_dir}/${sample_id}.merged.kraken.result \
	--use-names 
```

Unmerged reads were classified in paired-end mode:
```
kraken2 --paired \
	${unmerged_reads_dir}/${sample_id}.unmerged.non.host.R1.fastq.gz \
	${unmerged_reads_dir}/${sample_id}.unmerged.non.host.R2.fastq.gz \
	--db ${kraken2_db} \
	--confidence 0.0 \
	--report ${report_dir}/${sample_id}.unmerged.kraken.report \
	--output ${output_dir}/${sample_id}.unmerged.kraken.result \
	--use-names 

```

## Resistome profiling using AMR++

These steps were performed using the AMR++ v 3.0 ‘merged_resistome’ pipeline under default settings.

The pipeline was run to include single-nucleotide polymorphism (SNP) confirmation (‘-SNP Y’), but not deduplication. 
```
nextflow run main_AMR++.nf -profile local \
	--pipeline merged_resistome \
	--snp Y \
	--deduped N
```

## HUMAnN Functional Profiling 
The HUMAnN 4.0.0.alpha.1 version was installed following the tutorial:
https://docs.google.com/document/d/1rCx5JkuO7wCKWrL8_-UJx_FkopJAfcDFtZktgPspak0/edit?tab=t.0#heading=h.i5hn0zprhwld

First, merged and unmerged host-free reads were concatenated into one file per sample: 
```
concatenate_reads() {
    sample=$1
    output_file="${output_dir}/${sample}_concatenated.fastq.gz"

    # Find all matching files (merged and unmerged)
    files=(${input_dir}/${sample}.*.non.host.fastq.gz)

    #Concatenate them
    cat ${files[@]} > $output_file
}
```

Here, the ‘vOct22_CHOCOPhlAnSGB_202403’ database was used. 
The MetaPhlan4 database was downloaded from https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/
The Bowtie index for MetaPhlan4 was downloaded from http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/

DNA reads were processed first: 
```
humann \
        --input "$reads" \
        --output "$output_dir" \
        --output-basename "$sample_id" \
        --metaphlan-options  "-t rel_ab_w_read_stats --unclassified_estimation --index mpa_vOct22_CHOCOPhlAnSGB_202403" \
        --prescreen-threshold 0.5 \
        --threads "$threads" \
        --verbose
```

cDNA reads were processed next: 
```
humann \
        --input "$reads" \
        --output "$output_dir" \
        --output-basename "$rna_sample_id" \
        --taxonomic-profile "${output_dir}/${dna_sample_id}_1_metaphlan_profile.tsv" \
        --prescreen-threshold 0.5 \
        --threads "${threads}" \
        --verbose
```

The per-sample tables were joined into a single output using the 'humann_join_tables' utility tool: 
```
mkdir -p Final_count_tables

#Gene families
humann_join_tables \
	--input Functional_Profiling/ \
	--output Final_count_tables/humann_genefamilies.tsv \
	--file_name genefamilies.tsv \
	--verbose

##Pathway abundances
humann_join_tables \
	--input Functional_Profiling/ \
	--output Final_count_tables/humann_pathabundance.tsv \
	--file_name pathabundance.tsv \
	--verbose
```

UniRef90 gene families from each sample were grouped into KO groups:
```
humann_convert() {
    sample=$1
    humann_regroup_table \
        -i Functional_Profiling/${sample}_2_genefamilies.tsv \
        -c '${utility_mapping_humman4_update}/map_ko_uniref90.txt.gz'  \
        -o Functional_Profiling/${sample}_humann_KO_group_counts.tsv
}

#Once grouped, merge per sample tables
humann_join_tables \
	--input Functional_Profiling/ \
	--output Final_count_tables/humann_KO_group_counts.tsv \
	--file_name humann_KO_group_counts.tsv
```


The KO group and pathway count tables were subset for only the community-level abundances:
```
# Pathways
grep -v '|' Final_count_tables/humann_pathabundance.tsv \
	> Final_count_tables/HUMAnN_Pathways_Counts_Feedlot_CatchBasins.tsv
#Gene groups
grep -v '|' Final_count_tables/humann_KO_group_counts.tsv \
	> Final_count_tables/humann_KO_groupcounts_community.tsv
```

More informative names were given to the KO groups: 
```
humann_rename_table \
	-i Final_count_tables/humann_KO_groupcounts_community.tsv \	
	-c '${utility_mapping_humman4_update}/map_ko_name.txt.gz' \
	-o Final_count_tables/HUMAnN_KO_GroupCounts_Feedlot_CatchBasins.tsv
```

