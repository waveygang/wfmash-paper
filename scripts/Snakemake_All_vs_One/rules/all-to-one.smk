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
        delta-filter -m -i 90 -l 1000 {params.prefix}.delta > {params.prefix}_m_i90_l1k.delta
        paftools.js delta2paf {params.prefix}_m_i90_l1k.delta > {output.paf}
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
        wfmash {params.par} -n 1 -H 0.001 -t {threads} --hg-filter-ani-diff 30 {input.ref} {input.query} > {output.paf}
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
        fastq=lambda wildcards: "fastq/{species}/{sample}.hifi.fastq.gz".format(species=wildcards.species, sample=wildcards.sample)
    output:
        bam="results/02.reads_mapping/{species}/{sample}.hifi.sorted.bam"
    threads: 24
    shell:
        """
        minimap2 -ax map-hifi -t {threads} {input.ref} {input.fastq}|samtools sort -@ {threads} -O bam -o {output.bam}
        """
rule dv_prep:
    input:
        ref=lambda wildcards: "fasta/{species}/{ref}.fa.gz".format(species=wildcards.species, ref=config["references"][wildcards.species][0])
    output:
        fasta="results/02.reads_mapping/{species}/{ref}.fa",
        fai="results/02.reads_mapping/{species}/{ref}.fa.fai"
    threads: 2
    params:
    shell:
        """
        zcat {input.ref} > {output.fasta}
        samtools faidx {output.fasta}
        """

rule dv:
    input:
        ref=lambda wildcards: "results/02.reads_mapping/{species}/{ref}.fa".format(species=wildcards.species, ref=config["references"][wildcards.species][0]),
        bam=ancient("results/02.reads_mapping/{species}/{sample}.hifi.sorted.bam")
    output:
        gvcf="results/03.variants/dv/{species}/{sample}.g.vcf.gz",
        vcf="results/03.variants/dv/{species}/{sample}.vcf.gz"
    threads: 48
    params:
        indir=wkdir + "/results/02.reads_mapping/{species}",
        outdir=wkdir + "/results/03.variants/dv/{species}",
        sif=config["dv"]
    shell:
        """
        sample={wildcards.sample}
        REF_NAME=$(basename {input.ref})
        singularity exec -B {params.indir}:/input -B {params.outdir}:/output {params.sif} /bin/bash -c "/opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref /input/${{REF_NAME}} --reads /input/{wildcards.sample}.hifi.sorted.bam --output_vcf=/output/{wildcards.sample}.vcf.gz --output_gvcf=/output/{wildcards.sample}.g.vcf.gz --intermediate_results_dir=/output/{wildcards.sample}_tmp --num_shards={threads} --sample_name={wildcards.sample}"
        rm -rf {params.outdir}/{wildcards.sample}_tmp
        """

rule pafcall:
    input:
        ref=lambda wildcards: "fasta/{species}/{ref}.fa.gz".format(species=wildcards.species, ref=config["references"][wildcards.species][0]),
        query=lambda wildcards: "fasta/{species}/{sample}.fa.gz".format(species=wildcards.species, sample=wildcards.sample),
        paf="results/01.wga/{species}/{species}_{sample}_{aligner}_{params}.paf"
    output:
        sv="results/01.wga/{species}/{species}_{sample}_{aligner}_{params}.all_variants.filter.svs.vcf.gz",
        small="results/01.wga/{species}/{species}_{sample}_{aligner}_{params}.all_variants.filter.short.vcf.gz",
        raw="results/01.wga/{species}/{species}_{sample}_{aligner}_{params}.all_variants.raw.vcf.gz",
        site=temp("results/01.wga/{species}/{species}_{sample}_{aligner}_{params}.SiteDepth.gz"),
        bed="results/01.wga/{species}/{species}_{sample}_{aligner}_{params}.depth_eq1.bed",
        paf="results/01.wga/{species}/{species}_{sample}_{aligner}_{params}.filter.50kb.paf",
        maf=temp("results/01.wga/{species}/{species}_{sample}_{aligner}_{params}.filter.50kb.maf"),
        remaf="results/01.wga/{species}/{species}_{sample}_{aligner}_{params}.filter.50kb.rename.maf"
    params:
        prefix="results/01.wga/{species}/{species}_{sample}_{aligner}_{params}",
        pandepth=config["pandepth"],
        wgatools=config["wgatools"]
    threads: 4
    shell:
        """
        {params.wgatools} filter -f paf -a 50000 {input.paf}|sort -k1,1V -k3,3n > {output.paf}

        {params.pandepth} -i {output.paf} -a -o {params.prefix}
        zcat {output.site}|awk '$3<=1'|awk '{{print $1"\\t"$2"\\t"$2+1}}'|bedtools merge -d 1 > {output.bed}
        
        {params.wgatools} paf2maf -g {input.ref} -q {input.query} {output.paf} > {output.maf}

        REF_NAME=$(basename {input.ref} .fa.gz)
        {params.wgatools} rename --prefixs "${{REF_NAME}}#1#,{wildcards.sample}#1#" {output.maf} > {output.remaf}
        {params.wgatools} mi {output.remaf}
        {params.wgatools} call -s -l 1 -n {wildcards.sample} {output.remaf}|sed "s/${{REF_NAME}}#1#//g"|bgzip -@ {threads} -c > {output.raw}
        bcftools view -T {output.bed} {output.raw}|bcftools view -i 'SVLEN>=50' -O z -o {output.sv} 
        bcftools view -T {output.bed} {output.raw}|bcftools view -e 'SVLEN>=50' -O z -o {output.small} 
        
        """

rule sv_genotype:
    input:
        ref=lambda wildcards: "fasta/{species}/{ref}.fa.gz".format(species=wildcards.species, ref=config["references"][wildcards.species][0]),
        bam=lambda wildcards: "results/02.reads_mapping/{species}/{sample}.hifi.sorted.bam".format(species=wildcards.species, sample=wildcards.sample),
        vcf=rules.pafcall.output.sv
    output:
        vcf="results/03.variants/SVs/{species}/{species}_{sample}_{aligner}_{params}.cutesv.support.vcf.gz"
    threads: 12
    params:    
        vcf="results/03.variants/SVs/{species}/{species}_{sample}_{aligner}_{params}.cutesv.support.vcf",
        tmp="results/03.variants/SVs/{species}/{species}_{sample}_{aligner}_{params}_tmp"
    shell:
        """
        mkdir -p {params.tmp}
        cuteSV {input.bam} {input.ref} {params.vcf} {params.tmp} -s 1 --genotype -mi 500 -md 500 --min_mapq 0 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --min_size 50 -Ivcf {input.vcf}
        bgzip -@ {threads} {params.vcf}
        rm -rf {params.tmp}
        """

