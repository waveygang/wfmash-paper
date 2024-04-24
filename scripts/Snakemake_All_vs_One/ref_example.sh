Building DAG of jobs...
Job stats:
job                   count    min threads    max threads
------------------  -------  -------------  -------------
AnchorWave_align          2              1              1
AnchorWave_maf2paf        2              1              1
AnchorWave_pre            2              1              1
all                       1              1              1
dv                        2              1              1
dv_prep                   2              1              1
faindex                   2              1              1
hifi_mapping              2              1              1
minimap2                  2              1              1
nucmer                    2              1              1
pafcall                   8              1              1
sv_genotype               8              1              1
wfmash                    2              1              1
total                    37              1              1


[Wed Apr 24 22:01:32 2024]
rule nucmer:
    input: fasta/Slycopersicum/SL5.fa.gz, fasta/Slycopersicum/TS421.fa.gz
    output: results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.paf
    jobid: 36
    benchmark: benchmark/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.txt
    reason: Missing output files: results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.paf
    wildcards: species=Slycopersicum, sample=TS421, params=b500c500l100
    resources: tmpdir=/tmp


        nucmer --maxmatch -c 500 -b 500 -l 100 -p results/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100 -t 1 fasta/Slycopersicum/SL5.fa.gz fasta/Slycopersicum/TS421.fa.gz
        

[Wed Apr 24 22:01:32 2024]
rule dv_prep:
    input: fasta/Slycopersicum/SL5.fa.gz
    output: results/02.reads_mapping/Slycopersicum/SL5.fa, results/02.reads_mapping/Slycopersicum/SL5.fa.fai
    jobid: 22
    reason: Missing output files: results/02.reads_mapping/Slycopersicum/SL5.fa
    wildcards: species=Slycopersicum, ref=SL5
    resources: tmpdir=/tmp


        zcat fasta/Slycopersicum/SL5.fa.gz > results/02.reads_mapping/Slycopersicum/SL5.fa
        samtools faidx results/02.reads_mapping/Slycopersicum/SL5.fa
        

[Wed Apr 24 22:01:32 2024]
rule minimap2:
    input: fasta/Slycopersicum/SL5.fa.gz, fasta/Slycopersicum/TS421.fa.gz
    output: results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.paf
    jobid: 23
    benchmark: benchmark/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.txt
    reason: Missing output files: results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.paf
    wildcards: species=Slycopersicum, sample=TS421, params=asm20
    resources: tmpdir=/tmp


        minimap2 -x asm20 -c --eqx -t 1 fasta/Slycopersicum/SL5.fa.gz fasta/Slycopersicum/TS421.fa.gz > results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.paf
        

[Wed Apr 24 22:01:32 2024]
rule AnchorWave_pre:
    input: fasta/Athaliana/Col-CC.fa.gz, ref/Athaliana/Col-CC.gene.gff3
    output: ref/Athaliana/Col-CC.filter.cds.fa, ref/Athaliana/Col-CC.cds.sam
    jobid: 9
    benchmark: benchmark/Athaliana/Athaliana_Col-CC_AnchorWave_pre.txt
    reason: Missing output files: ref/Athaliana/Col-CC.filter.cds.fa, ref/Athaliana/Col-CC.cds.sam
    wildcards: species=Athaliana, ref=Col-CC
    resources: tmpdir=/tmp


        anchorwave gff2seq -i ref/Athaliana/Col-CC.gene.gff3 -r fasta/Athaliana/Col-CC.fa.gz -o ref/Athaliana/Col-CC.filter.cds.fa
        minimap2 -x splice -t 1 -k 12 -a -p 0.4 -N 20 fasta/Athaliana/Col-CC.fa.gz ref/Athaliana/Col-CC.filter.cds.fa > ref/Athaliana/Col-CC.cds.sam
        

[Wed Apr 24 22:01:32 2024]
rule hifi_mapping:
    input: fasta/Athaliana/Col-CC.fa.gz, fastq/Athaliana/Ler-0.hifi.fastq.gz
    output: results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam
    jobid: 1
    reason: Missing output files: results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam
    wildcards: species=Athaliana, sample=Ler-0
    resources: tmpdir=/tmp


        minimap2 -ax map-hifi -t 1 fasta/Athaliana/Col-CC.fa.gz fastq/Athaliana/Ler-0.hifi.fastq.gz|samtools sort -@ 1 -O bam -o results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam
        

[Wed Apr 24 22:01:32 2024]
rule nucmer:
    input: fasta/Athaliana/Col-CC.fa.gz, fasta/Athaliana/Ler-0.fa.gz
    output: results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.paf
    jobid: 17
    benchmark: benchmark/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.txt
    reason: Missing output files: results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.paf
    wildcards: species=Athaliana, sample=Ler-0, params=b500c500l100
    resources: tmpdir=/tmp


        nucmer --maxmatch -c 500 -b 500 -l 100 -p results/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100 -t 1 fasta/Athaliana/Col-CC.fa.gz fasta/Athaliana/Ler-0.fa.gz
        

[Wed Apr 24 22:01:32 2024]
rule faindex:
    input: fasta/Slycopersicum/TS421.fa.gz
    output: fasta/Slycopersicum/TS421.fa.gz.gzi, fasta/Slycopersicum/TS421.fa.gz.fai
    jobid: 33
    reason: Missing output files: fasta/Slycopersicum/TS421.fa.gz.fai
    wildcards: species=Slycopersicum, sample=TS421
    resources: tmpdir=/tmp


        samtools faidx fasta/Slycopersicum/TS421.fa.gz 
        

[Wed Apr 24 22:01:32 2024]
rule faindex:
    input: fasta/Slycopersicum/SL5.fa.gz
    output: fasta/Slycopersicum/SL5.fa.gz.gzi, fasta/Slycopersicum/SL5.fa.gz.fai
    jobid: 32
    reason: Missing output files: fasta/Slycopersicum/SL5.fa.gz.fai
    wildcards: species=Slycopersicum, sample=SL5
    resources: tmpdir=/tmp


        samtools faidx fasta/Slycopersicum/SL5.fa.gz 
        

[Wed Apr 24 22:01:32 2024]
rule wfmash:
    input: fasta/Athaliana/Col-CC.fa.gz, fasta/Athaliana/Col-CC.fa.gz.fai, fasta/Athaliana/Ler-0.fa.gz, fasta/Athaliana/Ler-0.fa.gz.fai
    output: results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.paf
    jobid: 12
    benchmark: benchmark/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.txt
    reason: Missing output files: results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.paf
    wildcards: species=Athaliana, sample=Ler-0, params=p80s10000c20000k19
    resources: tmpdir=/tmp


        wfmash -p 80 -s 10000 -c 20000 -k 19 -n 1 -H 0.001 -t 1 --hg-filter-ani-diff 30 fasta/Athaliana/Col-CC.fa.gz fasta/Athaliana/Ler-0.fa.gz > results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.paf
        

[Wed Apr 24 22:01:32 2024]
rule AnchorWave_pre:
    input: fasta/Slycopersicum/SL5.fa.gz, ref/Slycopersicum/SL5.gene.gff3
    output: ref/Slycopersicum/SL5.filter.cds.fa, ref/Slycopersicum/SL5.cds.sam
    jobid: 28
    benchmark: benchmark/Slycopersicum/Slycopersicum_SL5_AnchorWave_pre.txt
    reason: Missing output files: ref/Slycopersicum/SL5.cds.sam, ref/Slycopersicum/SL5.filter.cds.fa
    wildcards: species=Slycopersicum, ref=SL5
    resources: tmpdir=/tmp


        anchorwave gff2seq -i ref/Slycopersicum/SL5.gene.gff3 -r fasta/Slycopersicum/SL5.fa.gz -o ref/Slycopersicum/SL5.filter.cds.fa
        minimap2 -x splice -t 1 -k 12 -a -p 0.4 -N 20 fasta/Slycopersicum/SL5.fa.gz ref/Slycopersicum/SL5.filter.cds.fa > ref/Slycopersicum/SL5.cds.sam
        

[Wed Apr 24 22:01:32 2024]
rule dv_prep:
    input: fasta/Athaliana/Col-CC.fa.gz
    output: results/02.reads_mapping/Athaliana/Col-CC.fa, results/02.reads_mapping/Athaliana/Col-CC.fa.fai
    jobid: 3
    reason: Missing output files: results/02.reads_mapping/Athaliana/Col-CC.fa
    wildcards: species=Athaliana, ref=Col-CC
    resources: tmpdir=/tmp


        zcat fasta/Athaliana/Col-CC.fa.gz > results/02.reads_mapping/Athaliana/Col-CC.fa
        samtools faidx results/02.reads_mapping/Athaliana/Col-CC.fa
        

[Wed Apr 24 22:01:32 2024]
rule minimap2:
    input: fasta/Athaliana/Col-CC.fa.gz, fasta/Athaliana/Ler-0.fa.gz
    output: results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.paf
    jobid: 4
    benchmark: benchmark/Athaliana/Athaliana_Ler-0_minimap2_asm20.txt
    reason: Missing output files: results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.paf
    wildcards: species=Athaliana, sample=Ler-0, params=asm20
    resources: tmpdir=/tmp


        minimap2 -x asm20 -c --eqx -t 1 fasta/Athaliana/Col-CC.fa.gz fasta/Athaliana/Ler-0.fa.gz > results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.paf
        

[Wed Apr 24 22:01:32 2024]
rule hifi_mapping:
    input: fasta/Slycopersicum/SL5.fa.gz, fastq/Slycopersicum/TS421.hifi.fastq.gz
    output: results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam
    jobid: 20
    reason: Missing output files: results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam
    wildcards: species=Slycopersicum, sample=TS421
    resources: tmpdir=/tmp


        minimap2 -ax map-hifi -t 1 fasta/Slycopersicum/SL5.fa.gz fastq/Slycopersicum/TS421.hifi.fastq.gz|samtools sort -@ 1 -O bam -o results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam
        

[Wed Apr 24 22:01:32 2024]
rule AnchorWave_align:
    input: fasta/Athaliana/Col-CC.fa.gz, fasta/Athaliana/Ler-0.fa.gz, ref/Athaliana/Col-CC.filter.cds.fa, ref/Athaliana/Col-CC.cds.sam, ref/Athaliana/Col-CC.gene.gff3
    output: results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.cds.sam, results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.maf, results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.anchors
    jobid: 8
    benchmark: benchmark/Athaliana/Athaliana_Ler-0_AnchorWave_default.txt
    reason: Missing output files: results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.maf; Input files updated by another job: ref/Athaliana/Col-CC.filter.cds.fa, ref/Athaliana/Col-CC.cds.sam
    wildcards: species=Athaliana, sample=Ler-0, params=default
    resources: tmpdir=/tmp


        minimap2 -x splice -t 1 -k 12 -a -p 0.4 -N 20 fasta/Athaliana/Ler-0.fa.gz ref/Athaliana/Col-CC.filter.cds.fa > results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.cds.sam
        anchorwave genoali -IV -t 1 -i ref/Athaliana/Col-CC.gene.gff3 -as ref/Athaliana/Col-CC.filter.cds.fa -r fasta/Athaliana/Col-CC.fa.gz -a results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.cds.sam -ar ref/Athaliana/Col-CC.cds.sam -s fasta/Athaliana/Ler-0.fa.gz -n results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.anchors -o results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.maf -f results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.m.maf
        

[Wed Apr 24 22:01:32 2024]
rule pafcall:
    input: fasta/Slycopersicum/SL5.fa.gz, fasta/Slycopersicum/TS421.fa.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.paf
    output: results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.all_variants.filter.svs.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.all_variants.filter.short.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.all_variants.raw.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.SiteDepth.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.depth_eq1.bed, results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.filter.50kb.paf, results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.filter.50kb.maf, results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.filter.50kb.rename.maf
    jobid: 37
    reason: Missing output files: results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.all_variants.filter.svs.vcf.gz; Input files updated by another job: results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.paf
    wildcards: species=Slycopersicum, sample=TS421, aligner=nucmer, params=b500c500l100
    resources: tmpdir=/tmp


        ~/software/wgatools/wgatools filter -f paf -a 50000 results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.paf|sort -k1,1V -k3,3n > results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.filter.50kb.paf

        ~/software/PanDepth/bin/pandepth -i results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.filter.50kb.paf -a -o results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100
        zcat results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.SiteDepth.gz|awk '$3<=1'|awk '{print $1"\t"$2"\t"$2+1}'|bedtools merge -d 1 > results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.depth_eq1.bed
        
        ~/software/wgatools/wgatools paf2maf -g fasta/Slycopersicum/SL5.fa.gz -q fasta/Slycopersicum/TS421.fa.gz results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.filter.50kb.paf > results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.filter.50kb.maf

        REF_NAME=$(basename fasta/Slycopersicum/SL5.fa.gz .fa.gz)
        ~/software/wgatools/wgatools rename --prefixs "${REF_NAME}#1#,TS421#1#" results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.filter.50kb.maf > results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.filter.50kb.rename.maf
        ~/software/wgatools/wgatools mi results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.filter.50kb.rename.maf
        ~/software/wgatools/wgatools call -s -l 1 -n TS421 results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.filter.50kb.rename.maf|sed "s/${REF_NAME}#1#//g"|bgzip -@ 1 -c > results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.all_variants.raw.vcf.gz
        bcftools view -T results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.depth_eq1.bed results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.all_variants.raw.vcf.gz|bcftools view -i 'SVLEN>=50' -O z -o results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.all_variants.filter.svs.vcf.gz 
        bcftools view -T results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.depth_eq1.bed results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.all_variants.raw.vcf.gz|bcftools view -e 'SVLEN>=50' -O z -o results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.all_variants.filter.short.vcf.gz 
        
        

[Wed Apr 24 22:01:32 2024]
rule AnchorWave_align:
    input: fasta/Slycopersicum/SL5.fa.gz, fasta/Slycopersicum/TS421.fa.gz, ref/Slycopersicum/SL5.filter.cds.fa, ref/Slycopersicum/SL5.cds.sam, ref/Slycopersicum/SL5.gene.gff3
    output: results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.cds.sam, results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.maf, results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.anchors
    jobid: 27
    benchmark: benchmark/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.txt
    reason: Missing output files: results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.maf; Input files updated by another job: ref/Slycopersicum/SL5.cds.sam, ref/Slycopersicum/SL5.filter.cds.fa
    wildcards: species=Slycopersicum, sample=TS421, params=default
    resources: tmpdir=/tmp


        minimap2 -x splice -t 1 -k 12 -a -p 0.4 -N 20 fasta/Slycopersicum/TS421.fa.gz ref/Slycopersicum/SL5.filter.cds.fa > results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.cds.sam
        anchorwave genoali -IV -t 1 -i ref/Slycopersicum/SL5.gene.gff3 -as ref/Slycopersicum/SL5.filter.cds.fa -r fasta/Slycopersicum/SL5.fa.gz -a results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.cds.sam -ar ref/Slycopersicum/SL5.cds.sam -s fasta/Slycopersicum/TS421.fa.gz -n results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.anchors -o results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.maf -f results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.m.maf
        

[Wed Apr 24 22:01:32 2024]
rule pafcall:
    input: fasta/Athaliana/Col-CC.fa.gz, fasta/Athaliana/Ler-0.fa.gz, results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.paf
    output: results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.all_variants.filter.svs.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.all_variants.filter.short.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.all_variants.raw.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.SiteDepth.gz, results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.depth_eq1.bed, results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.filter.50kb.paf, results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.filter.50kb.maf, results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.filter.50kb.rename.maf
    jobid: 15
    reason: Missing output files: results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.all_variants.filter.svs.vcf.gz; Input files updated by another job: results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.paf
    wildcards: species=Athaliana, sample=Ler-0, aligner=wfmash, params=p80s10000c20000k19
    resources: tmpdir=/tmp


        ~/software/wgatools/wgatools filter -f paf -a 50000 results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.paf|sort -k1,1V -k3,3n > results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.filter.50kb.paf

        ~/software/PanDepth/bin/pandepth -i results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.filter.50kb.paf -a -o results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19
        zcat results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.SiteDepth.gz|awk '$3<=1'|awk '{print $1"\t"$2"\t"$2+1}'|bedtools merge -d 1 > results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.depth_eq1.bed
        
        ~/software/wgatools/wgatools paf2maf -g fasta/Athaliana/Col-CC.fa.gz -q fasta/Athaliana/Ler-0.fa.gz results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.filter.50kb.paf > results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.filter.50kb.maf

        REF_NAME=$(basename fasta/Athaliana/Col-CC.fa.gz .fa.gz)
        ~/software/wgatools/wgatools rename --prefixs "${REF_NAME}#1#,Ler-0#1#" results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.filter.50kb.maf > results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.filter.50kb.rename.maf
        ~/software/wgatools/wgatools mi results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.filter.50kb.rename.maf
        ~/software/wgatools/wgatools call -s -l 1 -n Ler-0 results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.filter.50kb.rename.maf|sed "s/${REF_NAME}#1#//g"|bgzip -@ 1 -c > results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.all_variants.raw.vcf.gz
        bcftools view -T results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.depth_eq1.bed results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.all_variants.raw.vcf.gz|bcftools view -i 'SVLEN>=50' -O z -o results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.all_variants.filter.svs.vcf.gz 
        bcftools view -T results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.depth_eq1.bed results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.all_variants.raw.vcf.gz|bcftools view -e 'SVLEN>=50' -O z -o results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.all_variants.filter.short.vcf.gz 
        
        

[Wed Apr 24 22:01:32 2024]
rule dv:
    input: results/02.reads_mapping/Athaliana/Col-CC.fa, results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam
    output: results/03.variants/dv/Athaliana/Ler-0.g.vcf.gz, results/03.variants/dv/Athaliana/Ler-0.vcf.gz
    jobid: 2
    reason: Missing output files: results/03.variants/dv/Athaliana/Ler-0.vcf.gz; Input files updated by another job: results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam, results/02.reads_mapping/Athaliana/Col-CC.fa
    wildcards: species=Athaliana, sample=Ler-0
    resources: tmpdir=/tmp


        sample=Ler-0
        REF_NAME=$(basename results/02.reads_mapping/Athaliana/Col-CC.fa)
        singularity exec -B /ebio/abt6_projects/AtGraph/tmp/wfmash/Snakemake-workflow/results/02.reads_mapping/Athaliana:/input -B /ebio/abt6_projects/AtGraph/tmp/wfmash/Snakemake-workflow/results/03.variants/dv/Athaliana:/output /ebio/abt6_projects7/small_projects/zbao/software/DeepVariant/dv-v1.6.0.sif /bin/bash -c "/opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref /input/${REF_NAME} --reads /input/Ler-0.hifi.sorted.bam --output_vcf=/output/Ler-0.vcf.gz --output_gvcf=/output/Ler-0.g.vcf.gz --intermediate_results_dir=/output/Ler-0_tmp --num_shards=1 --sample_name=Ler-0"
        rm -rf /ebio/abt6_projects/AtGraph/tmp/wfmash/Snakemake-workflow/results/03.variants/dv/Athaliana/Ler-0_tmp
        

[Wed Apr 24 22:01:32 2024]
rule wfmash:
    input: fasta/Slycopersicum/SL5.fa.gz, fasta/Slycopersicum/SL5.fa.gz.fai, fasta/Slycopersicum/TS421.fa.gz, fasta/Slycopersicum/TS421.fa.gz.fai
    output: results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.paf
    jobid: 31
    benchmark: benchmark/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.txt
    reason: Missing output files: results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.paf; Input files updated by another job: fasta/Slycopersicum/SL5.fa.gz.fai, fasta/Slycopersicum/TS421.fa.gz.fai
    wildcards: species=Slycopersicum, sample=TS421, params=p80s10000c20000k19
    resources: tmpdir=/tmp


        wfmash -p 80 -s 10000 -c 20000 -k 19 -n 1 -H 0.001 -t 1 --hg-filter-ani-diff 30 fasta/Slycopersicum/SL5.fa.gz fasta/Slycopersicum/TS421.fa.gz > results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.paf
        

[Wed Apr 24 22:01:32 2024]
rule pafcall:
    input: fasta/Athaliana/Col-CC.fa.gz, fasta/Athaliana/Ler-0.fa.gz, results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.paf
    output: results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.all_variants.filter.svs.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.all_variants.filter.short.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.all_variants.raw.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.SiteDepth.gz, results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.depth_eq1.bed, results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.filter.50kb.paf, results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.filter.50kb.maf, results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.filter.50kb.rename.maf
    jobid: 5
    reason: Missing output files: results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.all_variants.filter.svs.vcf.gz; Input files updated by another job: results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.paf
    wildcards: species=Athaliana, sample=Ler-0, aligner=minimap2, params=asm20
    resources: tmpdir=/tmp


        ~/software/wgatools/wgatools filter -f paf -a 50000 results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.paf|sort -k1,1V -k3,3n > results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.filter.50kb.paf

        ~/software/PanDepth/bin/pandepth -i results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.filter.50kb.paf -a -o results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20
        zcat results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.SiteDepth.gz|awk '$3<=1'|awk '{print $1"\t"$2"\t"$2+1}'|bedtools merge -d 1 > results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.depth_eq1.bed
        
        ~/software/wgatools/wgatools paf2maf -g fasta/Athaliana/Col-CC.fa.gz -q fasta/Athaliana/Ler-0.fa.gz results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.filter.50kb.paf > results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.filter.50kb.maf

        REF_NAME=$(basename fasta/Athaliana/Col-CC.fa.gz .fa.gz)
        ~/software/wgatools/wgatools rename --prefixs "${REF_NAME}#1#,Ler-0#1#" results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.filter.50kb.maf > results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.filter.50kb.rename.maf
        ~/software/wgatools/wgatools mi results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.filter.50kb.rename.maf
        ~/software/wgatools/wgatools call -s -l 1 -n Ler-0 results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.filter.50kb.rename.maf|sed "s/${REF_NAME}#1#//g"|bgzip -@ 1 -c > results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.all_variants.raw.vcf.gz
        bcftools view -T results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.depth_eq1.bed results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.all_variants.raw.vcf.gz|bcftools view -i 'SVLEN>=50' -O z -o results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.all_variants.filter.svs.vcf.gz 
        bcftools view -T results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.depth_eq1.bed results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.all_variants.raw.vcf.gz|bcftools view -e 'SVLEN>=50' -O z -o results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.all_variants.filter.short.vcf.gz 
        
        

[Wed Apr 24 22:01:32 2024]
rule dv:
    input: results/02.reads_mapping/Slycopersicum/SL5.fa, results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam
    output: results/03.variants/dv/Slycopersicum/TS421.g.vcf.gz, results/03.variants/dv/Slycopersicum/TS421.vcf.gz
    jobid: 21
    reason: Missing output files: results/03.variants/dv/Slycopersicum/TS421.vcf.gz; Input files updated by another job: results/02.reads_mapping/Slycopersicum/SL5.fa, results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam
    wildcards: species=Slycopersicum, sample=TS421
    resources: tmpdir=/tmp


        sample=TS421
        REF_NAME=$(basename results/02.reads_mapping/Slycopersicum/SL5.fa)
        singularity exec -B /ebio/abt6_projects/AtGraph/tmp/wfmash/Snakemake-workflow/results/02.reads_mapping/Slycopersicum:/input -B /ebio/abt6_projects/AtGraph/tmp/wfmash/Snakemake-workflow/results/03.variants/dv/Slycopersicum:/output /ebio/abt6_projects7/small_projects/zbao/software/DeepVariant/dv-v1.6.0.sif /bin/bash -c "/opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref /input/${REF_NAME} --reads /input/TS421.hifi.sorted.bam --output_vcf=/output/TS421.vcf.gz --output_gvcf=/output/TS421.g.vcf.gz --intermediate_results_dir=/output/TS421_tmp --num_shards=1 --sample_name=TS421"
        rm -rf /ebio/abt6_projects/AtGraph/tmp/wfmash/Snakemake-workflow/results/03.variants/dv/Slycopersicum/TS421_tmp
        

[Wed Apr 24 22:01:32 2024]
rule pafcall:
    input: fasta/Slycopersicum/SL5.fa.gz, fasta/Slycopersicum/TS421.fa.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.paf
    output: results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.all_variants.filter.svs.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.all_variants.filter.short.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.all_variants.raw.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.SiteDepth.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.depth_eq1.bed, results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.filter.50kb.paf, results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.filter.50kb.maf, results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.filter.50kb.rename.maf
    jobid: 24
    reason: Missing output files: results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.all_variants.filter.svs.vcf.gz; Input files updated by another job: results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.paf
    wildcards: species=Slycopersicum, sample=TS421, aligner=minimap2, params=asm20
    resources: tmpdir=/tmp


        ~/software/wgatools/wgatools filter -f paf -a 50000 results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.paf|sort -k1,1V -k3,3n > results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.filter.50kb.paf

        ~/software/PanDepth/bin/pandepth -i results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.filter.50kb.paf -a -o results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20
        zcat results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.SiteDepth.gz|awk '$3<=1'|awk '{print $1"\t"$2"\t"$2+1}'|bedtools merge -d 1 > results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.depth_eq1.bed
        
        ~/software/wgatools/wgatools paf2maf -g fasta/Slycopersicum/SL5.fa.gz -q fasta/Slycopersicum/TS421.fa.gz results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.filter.50kb.paf > results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.filter.50kb.maf

        REF_NAME=$(basename fasta/Slycopersicum/SL5.fa.gz .fa.gz)
        ~/software/wgatools/wgatools rename --prefixs "${REF_NAME}#1#,TS421#1#" results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.filter.50kb.maf > results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.filter.50kb.rename.maf
        ~/software/wgatools/wgatools mi results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.filter.50kb.rename.maf
        ~/software/wgatools/wgatools call -s -l 1 -n TS421 results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.filter.50kb.rename.maf|sed "s/${REF_NAME}#1#//g"|bgzip -@ 1 -c > results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.all_variants.raw.vcf.gz
        bcftools view -T results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.depth_eq1.bed results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.all_variants.raw.vcf.gz|bcftools view -i 'SVLEN>=50' -O z -o results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.all_variants.filter.svs.vcf.gz 
        bcftools view -T results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.depth_eq1.bed results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.all_variants.raw.vcf.gz|bcftools view -e 'SVLEN>=50' -O z -o results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.all_variants.filter.short.vcf.gz 
        
        

[Wed Apr 24 22:01:32 2024]
rule pafcall:
    input: fasta/Athaliana/Col-CC.fa.gz, fasta/Athaliana/Ler-0.fa.gz, results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.paf
    output: results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.all_variants.filter.svs.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.all_variants.filter.short.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.all_variants.raw.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.SiteDepth.gz, results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.depth_eq1.bed, results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.filter.50kb.paf, results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.filter.50kb.maf, results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.filter.50kb.rename.maf
    jobid: 18
    reason: Missing output files: results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.all_variants.filter.svs.vcf.gz; Input files updated by another job: results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.paf
    wildcards: species=Athaliana, sample=Ler-0, aligner=nucmer, params=b500c500l100
    resources: tmpdir=/tmp


        ~/software/wgatools/wgatools filter -f paf -a 50000 results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.paf|sort -k1,1V -k3,3n > results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.filter.50kb.paf

        ~/software/PanDepth/bin/pandepth -i results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.filter.50kb.paf -a -o results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100
        zcat results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.SiteDepth.gz|awk '$3<=1'|awk '{print $1"\t"$2"\t"$2+1}'|bedtools merge -d 1 > results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.depth_eq1.bed
        
        ~/software/wgatools/wgatools paf2maf -g fasta/Athaliana/Col-CC.fa.gz -q fasta/Athaliana/Ler-0.fa.gz results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.filter.50kb.paf > results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.filter.50kb.maf

        REF_NAME=$(basename fasta/Athaliana/Col-CC.fa.gz .fa.gz)
        ~/software/wgatools/wgatools rename --prefixs "${REF_NAME}#1#,Ler-0#1#" results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.filter.50kb.maf > results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.filter.50kb.rename.maf
        ~/software/wgatools/wgatools mi results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.filter.50kb.rename.maf
        ~/software/wgatools/wgatools call -s -l 1 -n Ler-0 results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.filter.50kb.rename.maf|sed "s/${REF_NAME}#1#//g"|bgzip -@ 1 -c > results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.all_variants.raw.vcf.gz
        bcftools view -T results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.depth_eq1.bed results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.all_variants.raw.vcf.gz|bcftools view -i 'SVLEN>=50' -O z -o results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.all_variants.filter.svs.vcf.gz 
        bcftools view -T results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.depth_eq1.bed results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.all_variants.raw.vcf.gz|bcftools view -e 'SVLEN>=50' -O z -o results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.all_variants.filter.short.vcf.gz 
        
        
Would remove temporary output results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.cds.sam
Would remove temporary output results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.anchors
Would remove temporary output results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.SiteDepth.gz
Would remove temporary output results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.filter.50kb.maf
Would remove temporary output results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.cds.sam
Would remove temporary output results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.anchors
Would remove temporary output results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.SiteDepth.gz
Would remove temporary output results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.filter.50kb.maf
Would remove temporary output results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.SiteDepth.gz
Would remove temporary output results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.filter.50kb.maf
Would remove temporary output results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.SiteDepth.gz
Would remove temporary output results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.filter.50kb.maf
Would remove temporary output results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.SiteDepth.gz
Would remove temporary output results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.filter.50kb.maf

[Wed Apr 24 22:01:32 2024]
rule sv_genotype:
    input: fasta/Slycopersicum/SL5.fa.gz, results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam, results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.all_variants.filter.svs.vcf.gz
    output: results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.cutesv.support.vcf.gz
    jobid: 25
    reason: Missing output files: results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.cutesv.support.vcf.gz; Input files updated by another job: results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.all_variants.filter.svs.vcf.gz, results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam
    wildcards: species=Slycopersicum, sample=TS421, aligner=minimap2, params=asm20
    resources: tmpdir=/tmp


        mkdir -p results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_minimap2_asm20_tmp
        cuteSV results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam fasta/Slycopersicum/SL5.fa.gz results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.cutesv.support.vcf results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_minimap2_asm20_tmp -s 1 --genotype -mi 500 -md 500 --min_mapq 0 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --min_size 50 -Ivcf results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.all_variants.filter.svs.vcf.gz
        bgzip -@ 1 results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.cutesv.support.vcf
        rm -rf results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_minimap2_asm20_tmp
        

[Wed Apr 24 22:01:32 2024]
rule sv_genotype:
    input: fasta/Athaliana/Col-CC.fa.gz, results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam, results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.all_variants.filter.svs.vcf.gz
    output: results/03.variants/SVs/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.cutesv.support.vcf.gz
    jobid: 19
    reason: Missing output files: results/03.variants/SVs/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.cutesv.support.vcf.gz; Input files updated by another job: results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam, results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.all_variants.filter.svs.vcf.gz
    wildcards: species=Athaliana, sample=Ler-0, aligner=nucmer, params=b500c500l100
    resources: tmpdir=/tmp


        mkdir -p results/03.variants/SVs/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100_tmp
        cuteSV results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam fasta/Athaliana/Col-CC.fa.gz results/03.variants/SVs/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.cutesv.support.vcf results/03.variants/SVs/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100_tmp -s 1 --genotype -mi 500 -md 500 --min_mapq 0 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --min_size 50 -Ivcf results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.all_variants.filter.svs.vcf.gz
        bgzip -@ 1 results/03.variants/SVs/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.cutesv.support.vcf
        rm -rf results/03.variants/SVs/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100_tmp
        

[Wed Apr 24 22:01:32 2024]
rule AnchorWave_maf2paf:
    input: results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.maf
    output: results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.paf
    jobid: 7
    reason: Missing output files: results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.paf; Input files updated by another job: results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.maf
    wildcards: species=Athaliana, sample=Ler-0, params=default
    resources: tmpdir=/tmp


        wgatools maf2paf results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.maf > results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.paf 
        

[Wed Apr 24 22:01:32 2024]
rule sv_genotype:
    input: fasta/Slycopersicum/SL5.fa.gz, results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam, results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.all_variants.filter.svs.vcf.gz
    output: results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.cutesv.support.vcf.gz
    jobid: 38
    reason: Missing output files: results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.cutesv.support.vcf.gz; Input files updated by another job: results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.all_variants.filter.svs.vcf.gz, results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam
    wildcards: species=Slycopersicum, sample=TS421, aligner=nucmer, params=b500c500l100
    resources: tmpdir=/tmp


        mkdir -p results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100_tmp
        cuteSV results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam fasta/Slycopersicum/SL5.fa.gz results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.cutesv.support.vcf results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100_tmp -s 1 --genotype -mi 500 -md 500 --min_mapq 0 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --min_size 50 -Ivcf results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.all_variants.filter.svs.vcf.gz
        bgzip -@ 1 results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.cutesv.support.vcf
        rm -rf results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100_tmp
        

[Wed Apr 24 22:01:32 2024]
rule AnchorWave_maf2paf:
    input: results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.maf
    output: results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.paf
    jobid: 26
    reason: Missing output files: results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.paf; Input files updated by another job: results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.maf
    wildcards: species=Slycopersicum, sample=TS421, params=default
    resources: tmpdir=/tmp


        wgatools maf2paf results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.maf > results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.paf 
        

[Wed Apr 24 22:01:32 2024]
rule sv_genotype:
    input: fasta/Athaliana/Col-CC.fa.gz, results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam, results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.all_variants.filter.svs.vcf.gz
    output: results/03.variants/SVs/Athaliana/Athaliana_Ler-0_minimap2_asm20.cutesv.support.vcf.gz
    jobid: 6
    reason: Missing output files: results/03.variants/SVs/Athaliana/Athaliana_Ler-0_minimap2_asm20.cutesv.support.vcf.gz; Input files updated by another job: results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam, results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.all_variants.filter.svs.vcf.gz
    wildcards: species=Athaliana, sample=Ler-0, aligner=minimap2, params=asm20
    resources: tmpdir=/tmp


        mkdir -p results/03.variants/SVs/Athaliana/Athaliana_Ler-0_minimap2_asm20_tmp
        cuteSV results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam fasta/Athaliana/Col-CC.fa.gz results/03.variants/SVs/Athaliana/Athaliana_Ler-0_minimap2_asm20.cutesv.support.vcf results/03.variants/SVs/Athaliana/Athaliana_Ler-0_minimap2_asm20_tmp -s 1 --genotype -mi 500 -md 500 --min_mapq 0 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --min_size 50 -Ivcf results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.all_variants.filter.svs.vcf.gz
        bgzip -@ 1 results/03.variants/SVs/Athaliana/Athaliana_Ler-0_minimap2_asm20.cutesv.support.vcf
        rm -rf results/03.variants/SVs/Athaliana/Athaliana_Ler-0_minimap2_asm20_tmp
        

[Wed Apr 24 22:01:32 2024]
rule sv_genotype:
    input: fasta/Athaliana/Col-CC.fa.gz, results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam, results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.all_variants.filter.svs.vcf.gz
    output: results/03.variants/SVs/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.cutesv.support.vcf.gz
    jobid: 16
    reason: Missing output files: results/03.variants/SVs/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.cutesv.support.vcf.gz; Input files updated by another job: results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam, results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.all_variants.filter.svs.vcf.gz
    wildcards: species=Athaliana, sample=Ler-0, aligner=wfmash, params=p80s10000c20000k19
    resources: tmpdir=/tmp


        mkdir -p results/03.variants/SVs/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19_tmp
        cuteSV results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam fasta/Athaliana/Col-CC.fa.gz results/03.variants/SVs/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.cutesv.support.vcf results/03.variants/SVs/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19_tmp -s 1 --genotype -mi 500 -md 500 --min_mapq 0 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --min_size 50 -Ivcf results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.all_variants.filter.svs.vcf.gz
        bgzip -@ 1 results/03.variants/SVs/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.cutesv.support.vcf
        rm -rf results/03.variants/SVs/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19_tmp
        

[Wed Apr 24 22:01:32 2024]
rule pafcall:
    input: fasta/Slycopersicum/SL5.fa.gz, fasta/Slycopersicum/TS421.fa.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.paf
    output: results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.all_variants.filter.svs.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.all_variants.filter.short.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.all_variants.raw.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.SiteDepth.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.depth_eq1.bed, results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.filter.50kb.paf, results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.filter.50kb.maf, results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.filter.50kb.rename.maf
    jobid: 34
    reason: Missing output files: results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.all_variants.filter.svs.vcf.gz; Input files updated by another job: results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.paf
    wildcards: species=Slycopersicum, sample=TS421, aligner=wfmash, params=p80s10000c20000k19
    resources: tmpdir=/tmp


        ~/software/wgatools/wgatools filter -f paf -a 50000 results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.paf|sort -k1,1V -k3,3n > results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.filter.50kb.paf

        ~/software/PanDepth/bin/pandepth -i results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.filter.50kb.paf -a -o results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19
        zcat results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.SiteDepth.gz|awk '$3<=1'|awk '{print $1"\t"$2"\t"$2+1}'|bedtools merge -d 1 > results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.depth_eq1.bed
        
        ~/software/wgatools/wgatools paf2maf -g fasta/Slycopersicum/SL5.fa.gz -q fasta/Slycopersicum/TS421.fa.gz results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.filter.50kb.paf > results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.filter.50kb.maf

        REF_NAME=$(basename fasta/Slycopersicum/SL5.fa.gz .fa.gz)
        ~/software/wgatools/wgatools rename --prefixs "${REF_NAME}#1#,TS421#1#" results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.filter.50kb.maf > results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.filter.50kb.rename.maf
        ~/software/wgatools/wgatools mi results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.filter.50kb.rename.maf
        ~/software/wgatools/wgatools call -s -l 1 -n TS421 results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.filter.50kb.rename.maf|sed "s/${REF_NAME}#1#//g"|bgzip -@ 1 -c > results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.all_variants.raw.vcf.gz
        bcftools view -T results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.depth_eq1.bed results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.all_variants.raw.vcf.gz|bcftools view -i 'SVLEN>=50' -O z -o results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.all_variants.filter.svs.vcf.gz 
        bcftools view -T results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.depth_eq1.bed results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.all_variants.raw.vcf.gz|bcftools view -e 'SVLEN>=50' -O z -o results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.all_variants.filter.short.vcf.gz 
        
        
Would remove temporary output results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.SiteDepth.gz
Would remove temporary output results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.filter.50kb.maf

[Wed Apr 24 22:01:32 2024]
rule sv_genotype:
    input: fasta/Slycopersicum/SL5.fa.gz, results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam, results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.all_variants.filter.svs.vcf.gz
    output: results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.cutesv.support.vcf.gz
    jobid: 35
    reason: Missing output files: results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.cutesv.support.vcf.gz; Input files updated by another job: results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.all_variants.filter.svs.vcf.gz, results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam
    wildcards: species=Slycopersicum, sample=TS421, aligner=wfmash, params=p80s10000c20000k19
    resources: tmpdir=/tmp


        mkdir -p results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19_tmp
        cuteSV results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam fasta/Slycopersicum/SL5.fa.gz results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.cutesv.support.vcf results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19_tmp -s 1 --genotype -mi 500 -md 500 --min_mapq 0 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --min_size 50 -Ivcf results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.all_variants.filter.svs.vcf.gz
        bgzip -@ 1 results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.cutesv.support.vcf
        rm -rf results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19_tmp
        

[Wed Apr 24 22:01:32 2024]
rule pafcall:
    input: fasta/Slycopersicum/SL5.fa.gz, fasta/Slycopersicum/TS421.fa.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.paf
    output: results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.all_variants.filter.svs.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.all_variants.filter.short.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.all_variants.raw.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.SiteDepth.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.depth_eq1.bed, results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.filter.50kb.paf, results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.filter.50kb.maf, results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.filter.50kb.rename.maf
    jobid: 29
    reason: Missing output files: results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.all_variants.filter.svs.vcf.gz; Input files updated by another job: results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.paf
    wildcards: species=Slycopersicum, sample=TS421, aligner=AnchorWave, params=default
    resources: tmpdir=/tmp


        ~/software/wgatools/wgatools filter -f paf -a 50000 results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.paf|sort -k1,1V -k3,3n > results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.filter.50kb.paf

        ~/software/PanDepth/bin/pandepth -i results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.filter.50kb.paf -a -o results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default
        zcat results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.SiteDepth.gz|awk '$3<=1'|awk '{print $1"\t"$2"\t"$2+1}'|bedtools merge -d 1 > results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.depth_eq1.bed
        
        ~/software/wgatools/wgatools paf2maf -g fasta/Slycopersicum/SL5.fa.gz -q fasta/Slycopersicum/TS421.fa.gz results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.filter.50kb.paf > results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.filter.50kb.maf

        REF_NAME=$(basename fasta/Slycopersicum/SL5.fa.gz .fa.gz)
        ~/software/wgatools/wgatools rename --prefixs "${REF_NAME}#1#,TS421#1#" results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.filter.50kb.maf > results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.filter.50kb.rename.maf
        ~/software/wgatools/wgatools mi results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.filter.50kb.rename.maf
        ~/software/wgatools/wgatools call -s -l 1 -n TS421 results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.filter.50kb.rename.maf|sed "s/${REF_NAME}#1#//g"|bgzip -@ 1 -c > results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.all_variants.raw.vcf.gz
        bcftools view -T results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.depth_eq1.bed results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.all_variants.raw.vcf.gz|bcftools view -i 'SVLEN>=50' -O z -o results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.all_variants.filter.svs.vcf.gz 
        bcftools view -T results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.depth_eq1.bed results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.all_variants.raw.vcf.gz|bcftools view -e 'SVLEN>=50' -O z -o results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.all_variants.filter.short.vcf.gz 
        
        

[Wed Apr 24 22:01:32 2024]
rule pafcall:
    input: fasta/Athaliana/Col-CC.fa.gz, fasta/Athaliana/Ler-0.fa.gz, results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.paf
    output: results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.all_variants.filter.svs.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.all_variants.filter.short.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.all_variants.raw.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.SiteDepth.gz, results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.depth_eq1.bed, results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.filter.50kb.paf, results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.filter.50kb.maf, results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.filter.50kb.rename.maf
    jobid: 10
    reason: Missing output files: results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.all_variants.filter.svs.vcf.gz; Input files updated by another job: results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.paf
    wildcards: species=Athaliana, sample=Ler-0, aligner=AnchorWave, params=default
    resources: tmpdir=/tmp


        ~/software/wgatools/wgatools filter -f paf -a 50000 results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.paf|sort -k1,1V -k3,3n > results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.filter.50kb.paf

        ~/software/PanDepth/bin/pandepth -i results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.filter.50kb.paf -a -o results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default
        zcat results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.SiteDepth.gz|awk '$3<=1'|awk '{print $1"\t"$2"\t"$2+1}'|bedtools merge -d 1 > results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.depth_eq1.bed
        
        ~/software/wgatools/wgatools paf2maf -g fasta/Athaliana/Col-CC.fa.gz -q fasta/Athaliana/Ler-0.fa.gz results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.filter.50kb.paf > results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.filter.50kb.maf

        REF_NAME=$(basename fasta/Athaliana/Col-CC.fa.gz .fa.gz)
        ~/software/wgatools/wgatools rename --prefixs "${REF_NAME}#1#,Ler-0#1#" results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.filter.50kb.maf > results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.filter.50kb.rename.maf
        ~/software/wgatools/wgatools mi results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.filter.50kb.rename.maf
        ~/software/wgatools/wgatools call -s -l 1 -n Ler-0 results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.filter.50kb.rename.maf|sed "s/${REF_NAME}#1#//g"|bgzip -@ 1 -c > results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.all_variants.raw.vcf.gz
        bcftools view -T results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.depth_eq1.bed results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.all_variants.raw.vcf.gz|bcftools view -i 'SVLEN>=50' -O z -o results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.all_variants.filter.svs.vcf.gz 
        bcftools view -T results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.depth_eq1.bed results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.all_variants.raw.vcf.gz|bcftools view -e 'SVLEN>=50' -O z -o results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.all_variants.filter.short.vcf.gz 
        
        
Would remove temporary output results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.SiteDepth.gz
Would remove temporary output results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.filter.50kb.maf
Would remove temporary output results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.SiteDepth.gz
Would remove temporary output results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.filter.50kb.maf

[Wed Apr 24 22:01:32 2024]
rule sv_genotype:
    input: fasta/Slycopersicum/SL5.fa.gz, results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam, results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.all_variants.filter.svs.vcf.gz
    output: results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.cutesv.support.vcf.gz
    jobid: 30
    reason: Missing output files: results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.cutesv.support.vcf.gz; Input files updated by another job: results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.all_variants.filter.svs.vcf.gz, results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam
    wildcards: species=Slycopersicum, sample=TS421, aligner=AnchorWave, params=default
    resources: tmpdir=/tmp


        mkdir -p results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_AnchorWave_default_tmp
        cuteSV results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam fasta/Slycopersicum/SL5.fa.gz results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.cutesv.support.vcf results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_AnchorWave_default_tmp -s 1 --genotype -mi 500 -md 500 --min_mapq 0 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --min_size 50 -Ivcf results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.all_variants.filter.svs.vcf.gz
        bgzip -@ 1 results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.cutesv.support.vcf
        rm -rf results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_AnchorWave_default_tmp
        

[Wed Apr 24 22:01:32 2024]
rule sv_genotype:
    input: fasta/Athaliana/Col-CC.fa.gz, results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam, results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.all_variants.filter.svs.vcf.gz
    output: results/03.variants/SVs/Athaliana/Athaliana_Ler-0_AnchorWave_default.cutesv.support.vcf.gz
    jobid: 11
    reason: Missing output files: results/03.variants/SVs/Athaliana/Athaliana_Ler-0_AnchorWave_default.cutesv.support.vcf.gz; Input files updated by another job: results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam, results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.all_variants.filter.svs.vcf.gz
    wildcards: species=Athaliana, sample=Ler-0, aligner=AnchorWave, params=default
    resources: tmpdir=/tmp


        mkdir -p results/03.variants/SVs/Athaliana/Athaliana_Ler-0_AnchorWave_default_tmp
        cuteSV results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam fasta/Athaliana/Col-CC.fa.gz results/03.variants/SVs/Athaliana/Athaliana_Ler-0_AnchorWave_default.cutesv.support.vcf results/03.variants/SVs/Athaliana/Athaliana_Ler-0_AnchorWave_default_tmp -s 1 --genotype -mi 500 -md 500 --min_mapq 0 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --min_size 50 -Ivcf results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.all_variants.filter.svs.vcf.gz
        bgzip -@ 1 results/03.variants/SVs/Athaliana/Athaliana_Ler-0_AnchorWave_default.cutesv.support.vcf
        rm -rf results/03.variants/SVs/Athaliana/Athaliana_Ler-0_AnchorWave_default_tmp
        

[Wed Apr 24 22:01:32 2024]
localrule all:
    input: results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam, results/03.variants/dv/Athaliana/Ler-0.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.paf, results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.all_variants.filter.svs.vcf.gz, results/03.variants/SVs/Athaliana/Athaliana_Ler-0_minimap2_asm20.cutesv.support.vcf.gz, results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam, results/03.variants/dv/Athaliana/Ler-0.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.paf, results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.all_variants.filter.svs.vcf.gz, results/03.variants/SVs/Athaliana/Athaliana_Ler-0_AnchorWave_default.cutesv.support.vcf.gz, results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam, results/03.variants/dv/Athaliana/Ler-0.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.paf, results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.all_variants.filter.svs.vcf.gz, results/03.variants/SVs/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.cutesv.support.vcf.gz, results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam, results/03.variants/dv/Athaliana/Ler-0.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.paf, results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.all_variants.filter.svs.vcf.gz, results/03.variants/SVs/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.cutesv.support.vcf.gz, results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam, results/03.variants/dv/Slycopersicum/TS421.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.paf, results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.all_variants.filter.svs.vcf.gz, results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.cutesv.support.vcf.gz, results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam, results/03.variants/dv/Slycopersicum/TS421.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.paf, results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.all_variants.filter.svs.vcf.gz, results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.cutesv.support.vcf.gz, results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam, results/03.variants/dv/Slycopersicum/TS421.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.paf, results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.all_variants.filter.svs.vcf.gz, results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.cutesv.support.vcf.gz, results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam, results/03.variants/dv/Slycopersicum/TS421.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.paf, results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.all_variants.filter.svs.vcf.gz, results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.cutesv.support.vcf.gz
    jobid: 0
    reason: Input files updated by another job: results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.paf, results/03.variants/dv/Slycopersicum/TS421.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.all_variants.filter.svs.vcf.gz, results/03.variants/SVs/Athaliana/Athaliana_Ler-0_minimap2_asm20.cutesv.support.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.paf, results/01.wga/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.all_variants.filter.svs.vcf.gz, results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.cutesv.support.vcf.gz, results/03.variants/SVs/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.cutesv.support.vcf.gz, results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam, results/03.variants/dv/Athaliana/Ler-0.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.all_variants.filter.svs.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.paf, results/03.variants/SVs/Athaliana/Athaliana_Ler-0_AnchorWave_default.cutesv.support.vcf.gz, results/02.reads_mapping/Athaliana/Ler-0.hifi.sorted.bam, results/01.wga/Slycopersicum/Slycopersicum_TS421_wfmash_p80s10000c20000k19.all_variants.filter.svs.vcf.gz, results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.cutesv.support.vcf.gz, results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_AnchorWave_default.cutesv.support.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.all_variants.filter.svs.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.paf, results/01.wga/Slycopersicum/Slycopersicum_TS421_nucmer_b500c500l100.all_variants.filter.svs.vcf.gz, results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.cutesv.support.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.all_variants.filter.svs.vcf.gz, results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.paf, results/01.wga/Athaliana/Athaliana_Ler-0_minimap2_asm20.paf, results/01.wga/Athaliana/Athaliana_Ler-0_wfmash_p80s10000c20000k19.all_variants.filter.svs.vcf.gz, results/01.wga/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.paf, results/01.wga/Athaliana/Athaliana_Ler-0_AnchorWave_default.paf, results/03.variants/SVs/Athaliana/Athaliana_Ler-0_nucmer_b500c500l100.cutesv.support.vcf.gz
    resources: tmpdir=/tmp

Job stats:
job                   count    min threads    max threads
------------------  -------  -------------  -------------
AnchorWave_align          2              1              1
AnchorWave_maf2paf        2              1              1
AnchorWave_pre            2              1              1
all                       1              1              1
dv                        2              1              1
dv_prep                   2              1              1
faindex                   2              1              1
hifi_mapping              2              1              1
minimap2                  2              1              1
nucmer                    2              1              1
pafcall                   8              1              1
sv_genotype               8              1              1
wfmash                    2              1              1
total                    37              1              1

Reasons:
    (check individual jobs above for details)
    input files updated by another job:
        AnchorWave_align, AnchorWave_maf2paf, all, dv, pafcall, sv_genotype, wfmash
    missing output files:
        AnchorWave_align, AnchorWave_maf2paf, AnchorWave_pre, dv, dv_prep, faindex, hifi_mapping, minimap2, nucmer, pafcall, sv_genotype, wfmash

This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
