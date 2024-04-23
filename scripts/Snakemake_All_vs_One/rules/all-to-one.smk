rule faindex:
    input:
        query=lambda wildcards: "fasta/{species}/{sample}.fa.gz".format(species=wildcards.species, sample=wildcards.sample)
    output:
        gzi="fasta/{species}/{sample}.fa.gz.gzi",
        fai="fasta/{species}/{sample}.fa.gz.fai"
    threads: 4
    shell:
        """
        samtools faidx {input.query} 
        """


rule minimap2:
    input:
        ref=lambda wildcards: "fasta/{species}/{ref}.fa.gz".format(species=wildcards.species, ref=config["references"][wildcards.species][0]),
        query=lambda wildcards: "fasta/{species}/{sample}.fa.gz".format(species=wildcards.species, sample=wildcards.sample)
    output:
        paf="results/01.wga/{species}/{species}_{sample}_minimap2_{params}.paf"
    params:
        par=lambda wildcards: untangle_parameters("minimap2", config["aligners"][wildcards.species]["minimap2"])
    threads: 12
    benchmark:
        "benchmark/{species}/{species}_{sample}_minimap2_{params}.txt"
    shell:
        """
        minimap2 {params.par} -c --eqx -t {threads} {input.ref} {input.query} > {output.paf}
        """

rule nucmer:
    input:
        ref=lambda wildcards: "fasta/{species}/{ref}.fa.gz".format(species=wildcards.species, ref=config["references"][wildcards.species][0]),
        query=lambda wildcards: "fasta/{species}/{sample}.fa.gz".format(species=wildcards.species, sample=wildcards.sample)
    output:
        paf="results/01.wga/{species}/{species}_{sample}_nucmer_{params}.paf"
    params:
        par=lambda wildcards: untangle_parameters("nucmer", config["aligners"][wildcards.species]["nucmer"]),
        prefix="results/{species}/{species}_{sample}_nucmer_{params}"
    threads: 12
    benchmark:
        "benchmark/{species}/{species}_{sample}_nucmer_{params}.txt"
    shell:
        """
        nucmer --maxmatch {params.par} -p {params.prefix} -t {threads} {input.ref} {input.query}
        """

rule wfmash:
    input:
        ref=lambda wildcards: "fasta/{species}/{ref}.fa.gz".format(species=wildcards.species, ref=config["references"][wildcards.species][0]),
        ref_fai=lambda wildcards: "fasta/{species}/{ref}.fa.gz.fai".format(species=wildcards.species, ref=config["references"][wildcards.species][0]),
        query=lambda wildcards: "fasta/{species}/{sample}.fa.gz".format(species=wildcards.species, sample=wildcards.sample),
        query_fai=lambda wildcards: "fasta/{species}/{sample}.fa.gz.fai".format(species=wildcards.species, sample=wildcards.sample)
    output:
        paf="results/01.wga/{species}/{species}_{sample}_wfmash_{params}.paf"
    params:
        par=lambda wildcards: untangle_parameters("wfmash", config["aligners"][wildcards.species]["wfmash"])
    threads: 12
    benchmark:
        "benchmark/{species}/{species}_{sample}_wfmash_{params}.txt"
    shell:
        """
        wfmash {params.par} -n 1 -k 19 -H 0.001 -t {threads} --hg-filter-ani-diff 30 {input.ref} {input.query} > {output.paf}
        """

rule AnchorWave_pre:
    input:
        ref=lambda wildcards: "fasta/{species}/{ref}.fa.gz".format(species=wildcards.species, ref=config["references"][wildcards.species][0]),
        ref_gff=lambda wildcards: "ref/{species}/{ref}.gene.gff3".format(species=wildcards.species, ref=config["references"][wildcards.species][0])
    output:
        ref_cds="ref/{species}/{ref}.filter.cds.fa",
        ref_sam="ref/{species}/{ref}.cds.sam"
    threads: 8
    benchmark:
        "benchmark/{species}/{species}_{ref}_AnchorWave_pre.txt"
    shell:
        """
        anchorwave gff2seq -i {input.ref_gff} -r {input.ref} -o {output.ref_cds}
        minimap2 -x splice -t {threads} -k 12 -a -p 0.4 -N 20 {input.ref} {output.ref_cds} > {output.ref_sam}
        """


rule AnchorWave_align:
    input:
        ref=lambda wildcards: "fasta/{species}/{ref}.fa.gz".format(species=wildcards.species, ref=config["references"][wildcards.species][0]),
        query=lambda wildcards: "fasta/{species}/{sample}.fa.gz".format(species=wildcards.species, sample=wildcards.sample),
        ref_cds=lambda wildcards: "ref/{species}/{ref}.filter.cds.fa".format(species=wildcards.species, ref=config["references"][wildcards.species][0]),
        ref_sam=lambda wildcards: "ref/{species}/{ref}.cds.sam".format(species=wildcards.species, ref=config["references"][wildcards.species][0]),
        ref_gff=lambda wildcards: "ref/{species}/{ref}.gene.gff3".format(species=wildcards.species, ref=config["references"][wildcards.species][0])
    output:
        sam=temp("results/01.wga/{species}/{species}_{sample}_AnchorWave_{params}.cds.sam"),
        maf="results/01.wga/{species}/{species}_{sample}_AnchorWave_{params}.maf",
        anchors=temp("results/01.wga/{species}/{species}_{sample}_AnchorWave_{params}.anchors")
    threads: 8
    params: 
        mmaf="results/01.wga/{species}/{species}_{sample}_AnchorWave_{params}.m.maf",
        par=lambda wildcards: untangle_parameters("AnchorWave", config["aligners"][wildcards.species]["AnchorWave"])
    benchmark:
        "benchmark/{species}/{species}_{sample}_AnchorWave_{params}.txt"
    shell:
        """
        minimap2 -x splice -t {threads} -k 12 -a -p 0.4 -N 20 {input.query} {input.ref_cds} > {output.sam}
        anchorwave {params.par} -t {threads} -i {input.ref_gff} -as {input.ref_cds} -r {input.ref} -a {output.sam} -ar {input.ref_sam} -s {input.query} -n {output.anchors} -o {output.maf} -f {params.mmaf}
        """

rule AnchorWave_maf2paf:
    input:
        maf=rules.AnchorWave_align.output.maf
    output: 
        paf="results/01.wga/{species}/{species}_{sample}_AnchorWave_{params}.paf"
    threads: 4
    shell:
        """
        wgatools maf2paf {input.maf} > {output.paf} 
        """

rule hifi_mapping:
    input:
        ref=lambda wildcards: "fasta/{species}/{ref}.fa.gz".format(species=wildcards.species, ref=config["references"][wildcards.species][0]),
        query=lambda wildcards: "fastq/{species}/{sample}.fastq.gz".format(species=wildcards.species, sample=wildcards.sample)
    output:
        bam="results/02.reads_mapping/{species}/{sample}.hifi.sorted.bam"
    threads: 24
    shell:
        """
        minimap2 -ax map-hifi -t {threads} {input.ref} {input.fastq}|samtools sort -@ {threads} -O bam -o {output.bam}
        """
rule dv:
    input:
        ref=lambda wildcards: "fasta/{species}/{ref}.fa.gz".format(species=wildcards.species, ref=config["references"][wildcards.species][0]),
        bam=ancient("results/02.reads_mapping/{species}/{sample}.hifi.sorted.bam")
    output:
        gvcf="results/03.variants/small/dv/{species}/{sample}.g.vcf.gz",
        vcf="results/03.variants/small/dv/{species}/{sample}.vcf.gz"
    threads: 48
    params:
        indir=wkdir + "/results/02.reads_mapping",
        outdir=wkdir + "/results/03.variants/small/dv/",
        sif=config["dv"]
    shell:
        """
        sample={wildcards.sample}
        singularity exec -B {params.indir}:/input -B {params.outdir}:/output {params.sif} /bin/bash -c "/opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref /input/{ref}.fa --reads /input/{wildcards.sample}.sorted.bam --regions {wildcards.chr} --output_vcf=/output/{wildcards.sample}.{wildcards.chr}.vcf.gz --output_gvcf=/output/{wildcards.sample}.{wildcards.chr}.g.vcf.gz --intermediate_results_dir=/output/{wildcards.sample}_{wildcards.chr} --num_shards={threads} --sample_name={wildcards.sample}"
        rm -rf {params.outdir}/{wildcards.sample}_{wildcards.chr}
        """

rule pafcall:
    input:
        ref=lambda wildcards: "fasta/{species}/{ref}.fa.gz".format(species=wildcards.species, ref=config["references"][wildcards.species][0]),
        query=lambda wildcards: "fastq/{species}/{sample}.fastq.gz".format(species=wildcards.species, sample=wildcards.sample),
        paf="results/01.wga/{species}/{species}_{sample}_{aligner}_{params}.paf"
    output:
        bam="results/02.reads_mapping/{species}/{sample}.hifi.sorted.bam"
    params:
        paf="results/01.wga/{species}/{species}_{sample}_{aligner}_{params}.filter.50kb.paf"
    threads: 24
    shell:
        """
        wgatools filter -p paf -i {input.paf} -r 50000 |sort -k1,1V -k3,3n > {params.paf}

        pandepth -i ./${sample}_${name}.filter.50kb.paf -a -o ${sample}_${name} 
        zcat ${sample}_${name}.SiteDepth.gz|awk '$3>1'|awk '{print $1"\t"$2"\t"$2+1}'|bedtools merge -d 1 > ${sample}_${name}_depth_lt1.bed
        zcat ${sample}_${name}.SiteDepth.gz|awk '$3==1'|awk '{print $1"\t"$2"\t"$2+1}'|bedtools merge -d 1 > ${sample}_${name}_depth_eq1.bed
        zcat ${sample}_${name}.SiteDepth.gz|awk '$3<1'|awk '{print $1"\t"$2"\t"$2+1}'|bedtools merge -d 1 > ${sample}_${name}_depth_ls1.bed
        
        wgatools paf2maf -g {input.ref} -q {input.query} ${sample}_${name}.filter.50kb.paf > ${sample}_${name}.filter.50kb.maf
        wgatools rename --prefixs "{wildcards.ref}#1#,{wildcardssample}#1#" ${sample}_${name}.filter.50kb.maf > ${sample}_${name}.filter.50kb.rename.maf
        wgatools mi ${sample}_${name}.filter.50kb.rename.maf
        wgatools call -s -l 1 -n ${sample}_${name} ${sample}_${name}.filter.50kb.rename.maf|sed "s/Col-CC#1#//g"|bgzip -@ 12 -c > ${sample}_${name}.vcf.gz
        bcftools view -T ^${sample}_${name}_depth_lt1.bed -O z -o ${sample}_${name}.filter.vcf.gz ${sample}_${name}.vcf.gz
        rm ${sample}_${name}.SiteDepth.gz
        rm ${sample}_${name}.filter.50kb.rename.maf ${sample}_${name}.filter.50kb.maf
        bgzip -@ 12 ${sample}_${name}.filter.50kb.paf
        """
    


rule sv_genotype:
    input:
        ref=lambda wildcards: "fasta/{species}/{ref}.fa.gz".format(species=wildcards.species, ref=config["references"][wildcards.species][0]),
        bam=lambda wildcards: "results/02.reads_mapping/{species}/{sample}.hifi.sorted.bam".format(species=wildcards.species, sample=wildcards.sample),
        vcf=rules.pafcall.outputs.vcf
    output:
        vcf="results/03.variants/SVs/{species}/{species}_{sample}_{aligner}_{params}.cutesv.support.vcf.gz"
    threads: 12
    params:    
        vcf="results/03.variants/SVs/{species}/{species}_{sample}_{aligner}_{params}.cutesv.support.vcf"
    shell:
        """
        cuteSV {input.bam} {input.ref} {params.vcf} {params.tmp} -s 1 --genotype -mi 500 -md 500 --min_mapq 0 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --min_size 50 -Ivcf {input.vcf}
        bgzip -@ {threads} {params.vcf}
        """

