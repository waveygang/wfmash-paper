# Evaluation

Variables:

```shell
DIR_BASE=/lizardfs/guarracino/wfmash-paper

PANDEPTH=/lizardfs/guarracino/git/PanDepth/bin/pandepth
BIN_VERSION="1.6.1" # DeepVariant version
THREADS=48

#WFMASH=/lizardfs/guarracino/git/wfmash/build/bin/wfmash-d7b696087f634f25e4b3de7dd521e1c4bfa3cf0e
#WFMASH=/lizardfs/guarracino/git/wfmash/build/bin/wfmash-c18520b786c2e8b1f1d19eb153d0908aa0a9b757
#WFMASH=/lizardfs/guarracino/git/wfmash/build/bin/wfmash-042386f01e6a5fdee5bb7d325529f035fbad0a29
WFMASH=/lizardfs/guarracino/git/wfmash/build/bin/wfmash-251f4e1d7770723d3495e9d7b443f941798d4174
WFMASH_PLOTS=/lizardfs/guarracino/git/wfmash/build/bin/wfmash-251f4e1d7770723d3495e9d7b443f941798d4174-with-plots

conda activate /lizardfs/guarracino/condatools/wfmash-paper
```

## Preparation

Tools:

```bash
conda create --prefix /lizardfs/guarracino/condatools/wfmash-paper/ -c conda-forge -c bioconda anchorwave=1.2.3 minimap2=2.28 wgatools=0.1.0 bcftools=1.20 samtools=1.20.0 syri=1.6.3 mummer4=4.0.0rc1 cutesv=2.1.1 -y

cd /lizardfs/guarracino/git

# PanDepth
git clone https://github.com/HuiyangYu/PanDepth.git
cd PanDepth
make


# DeepVariant (On lambda01)
wget https://github.com/apptainer/apptainer/releases/download/v1.2.5/apptainer_1.2.5_amd64.deb
sudo apt install -y ./apptainer_1.2.5_amd64.deb
rm ./apptainer_1.2.5_amd64.deb
# From https://github.com/google/deepvariant/blob/r1.6.1/docs/deepvariant-quick-start.md
apptainer pull docker://google/deepvariant:"${BIN_VERSION}"
```

Data:

```shell
# Assemblies
mkdir -p $DIR_BASE/assemblies/athaliana
cd $DIR_BASE/assemblies/athaliana
ATHALIANA82_FASTA=/lizardfs/guarracino/pggb-paper/assemblies/athaliana/athaliana82.fa.gz
SED_CMD=""
while IFS= read -r line; do
    search=$(echo "$line" | awk '{print $1}')
    replace=$(echo "$line" | awk '{print $2}')
    SED_CMD+="s/>$search/>$replace/g;"
done < $DIR_BASE/data/athaliana.contig2name.tsv
grep '^Athaliana' /lizardfs/guarracino/wfmash-paper/data/Atha_Sly_Hsp.data.tsv | cut -f 3 | while read ACC; do
  echo $ACC

  samtools faidx $ATHALIANA82_FASTA $(grep "^$ACC" $ATHALIANA82_FASTA.fai | cut -f 1) | sed "$SED_CMD" > $DIR_BASE/assemblies/athaliana/$ACC.fasta && samtools faidx $DIR_BASE/assemblies/athaliana/$ACC.fasta
done


mkdir -p $DIR_BASE/assemblies/tomato
cd $DIR_BASE/assemblies/tomato
grep '^Sly' /lizardfs/guarracino/wfmash-paper/data/Atha_Sly_Hsp.data.tsv | cut -f 3 | while read URL; do
  wget -c $URL
done
ls $DIR_BASE/assemblies/tomato/*fasta.gz | while read FASTA; do
  #NAME=$(basename $FASTA .fasta.gz | sed 's/\.0//g')
  #echo $FASTA $NAME
  #mv $FASTA x.gz
  #zcat x.gz | sed "s/>/>$NAME#1#chr/g" | bgzip -@ 48 -l 9 > $NAME.fa && samtools faidx $NAME.fa
  #rm x.gz
  gunzip $FASTA
done


# HiFi reads
mkdir -p $DIR_BASE/hifi/athaliana
cd $DIR_BASE/hifi/athaliana
grep '^Athaliana' /lizardfs/guarracino/wfmash-paper/data/Atha_Sly_Hsp.data.tsv | cut -f 3,4 | while read ACC URL; do
  echo $URL
  # Check if the fastq column contains a URL
  if [[ $URL == ftp* ]]; then
      wget $URL -P $DIR_BASE/hifi/athaliana
      NAME=$(basename $URL)
      mv $DIR_BASE/hifi/athaliana/$NAME $DIR_BASE/hifi/athaliana/$ACC.fq.gz
  fi
done

mkdir -p $DIR_BASE/hifi/tomato
cd $DIR_BASE/hifi/tomato
grep '^Sly' /lizardfs/guarracino/wfmash-paper/data/Atha_Sly_Hsp.data.tsv | cut -f 2,4 | while read ACC URL; do
  echo $URL
  # Check if the fastq column contains a URL
  if [[ $URL == ftp* ]]; then
      wget $URL -P $DIR_BASE/hifi/tomato
      NAME=$(basename $URL)
      mv $DIR_BASE/hifi/tomato/$NAME $DIR_BASE/hifi/tomato/$ACC.fq.gz
  fi
done


# Annotation
mkdir -p $DIR_BASE/annotation/athaliana
cd $DIR_BASE/annotation/athaliana
cat Col-CC.liffoff.gff3 | sed -e 's/^Chr1/GCA_028009825.2#1#CP116280.1/g' -e 's/^Chr2/GCA_028009825.2#1#CP116281.2/g' -e 's/^Chr3/GCA_028009825.2#1#CP116282.1/g' -e 's/^Chr4/GCA_028009825.2#1#CP116283.2/g' -e 's/^Chr5/GCA_028009825.2#1#CP116284.1/g' > GCA_028009825.2.gff3

mkdir -p $DIR_BASE/annotation/tomato
cd $DIR_BASE/annotation/tomato
wget -c https://solgenomics.net/ftp/genomes/TGG/gff/SL5.0.gff3.gz
gunzip SL5.0.gff3.gz
#zcat SL5.0.gff3.gz | sed 's/^/SL5#1#chr/g' | sort -k 1,1 -k 4,4n -k 5,5n > SL5.gff3
#rm SL5.0.gff3.gz
```

## Evaluation

### All-vs-all performance

```shell
#########################
PANGENOME=athaliana
PANGENOME=tomato
#########################

# List all fasta files in the specified directory
FASTA_FILES=($(ls $DIR_BASE/assemblies/$PANGENOME/*fasta))

# Get the number of files
NUM_FILES=${#FASTA_FILES[@]}

# Upper triangular matrix
# for ((i=1; i < NUM_FILES; i++)); do
#     for ((j=i + 1; j < NUM_FILES + 1; j++)); do
#         FASTA1="${FASTA_FILES[i]}"
#         FASTA2="${FASTA_FILES[j]}"

# Loop through the FASTA files to get only the lower triangular matrix of pairs
# Lower matrix to have GCA_028009825.2 and SL5.0 as target assemblies

for ((i = 1; i < NUM_FILES + 1; i++)); do
    for ((j = 1; j < i; j++)); do
        FASTA1="${FASTA_FILES[i]}"
        FASTA2="${FASTA_FILES[j]}"

        if [[ -f "$FASTA1" && -f "$FASTA2" ]]; then
          c=20k
          for k in 13 17; do
          #for k in 13 17 19; do
            for filter in "no1to1" "yes1to1"; do
              for p in 95 90 80 70; do
                for s in 10k 5k; do
                  for l in 0 5s; do
                    for hg in 0; do
                      echo $FASTA1 $FASTA2 $p $s $l $c $k $hg $filter
                    done
                  done
                done
              done
            done
          done
        fi
    done
done | tr ' ' '\t' > $DIR_BASE/$PANGENOME.combinations.tsv

mkdir -p $DIR_BASE/alignments/$PANGENOME
cd $DIR_BASE/alignments/$PANGENOME
sbatch -c 48 -p allnodes -x octopus03,octopus05 --array=1-$(wc -l < $DIR_BASE/$PANGENOME.combinations.tsv)%48 $DIR_BASE/scripts/job-array.wfmash.sh $WFMASH $DIR_BASE/athaliana.combinations.tsv $DIR_BASE/alignments/$PANGENOME $THREADS

for ((i = 1; i < NUM_FILES + 1; i++)); do
    for ((j = 1; j < i; j++)); do
        FASTA1="${FASTA_FILES[i]}"
        FASTA2="${FASTA_FILES[j]}"

        if [[ -f "$FASTA1" && -f "$FASTA2" ]]; then
          NAME1=$(basename $FASTA1 .fasta)
          NAME2=$(basename $FASTA2 .fasta)
          echo $NAME1 $NAME2

          # minimap2
          DIR_OUTPUT_MM2=$DIR_BASE/alignments/$PANGENOME/minimap2.xasm20.c.eqx.secondary-no/$NAME1-vs-$NAME2
          mkdir -p $DIR_OUTPUT_MM2
          cd $DIR_OUTPUT_MM2
          sbatch -c $THREADS -p allnodes --job-name "$NAME1-vs-$NAME2-minimap2" --wrap "hostname; cd /scratch; \time -v minimap2 -x asm20 -c --eqx --secondary=no -t $THREADS $FASTA2 $FASTA1 > $NAME1-vs-$NAME2.minimap2.asm20.paf; mv $NAME1-vs-$NAME2.minimap2.asm20.paf $DIR_OUTPUT_MM2"

          # # wfmash
          # c=20k
          # k=19
          # for p in 95 90 80 70; do
          #   for s in 10k 5k; do
          #     for hg in 30 0; do
          #       DIR_OUTPUT_WF=$DIR_BASE/alignments/$PANGENOME/wfmash.p${p}.s${s}.c${c}.k${k}.hg${hg}/$NAME1-vs-$NAME2
          #       mkdir -p $DIR_OUTPUT_WF
          #       cd $DIR_OUTPUT_WF
          #       sbatch -c $THREADS -p allnodes --job-name "$NAME1-vs-$NAME2-wfmash" --wrap "hostname; cd /scratch; \time -v $WFMASH -p $p -s $s -c $c -k $k --hg-filter-ani-diff $hg -t $THREADS $FASTA2 $FASTA1 > $NAME1-vs-$NAME2.wfmash.p${p}s${s}c${c}k${k}.hg${hg}.paf; mv $NAME1-vs-$NAME2.wfmash.p${p}s${s}c${c}k${k}.hg${hg}.paf $DIR_OUTPUT_WF"
          #     done
          #   done
          # done
        fi
    done
done




# Get runtime and memory
#UPDATE log2info TO MANAGE ALSO -l and [yes|no]1to1
(cat $DIR_BASE/alignments/$PANGENOME/slurm*; cat $DIR_BASE/alignments/$PANGENOME/*/*/slurm*.out) | python3 /lizardfs/guarracino/wfmash-paper/scripts/log2info.py > $DIR_BASE/alignments/$PANGENOME/$PANGENOME.runtime+memory.tsv

ls $DIR_BASE/alignments/$PANGENOME | grep '.tsv' -v | grep '.out' -v | while read TOOL; do
  ls $DIR_BASE/alignments/$PANGENOME/$TOOL/*/*paf | while read PAF; do
    NAME=$(basename $PAF .paf)

    #NUM_ALIGNMENT=$(cat $PAF | wc -l)
    AVG_ALIGNMENT_LEN=$(awk '{lenq+=($4-$3); lent+=($9-$8); n+=1} END {printf "%.0f\n", ((lenq + lent) / 2)/n}' $PAF)
    #TOTAL_ALIGNMENT_LEN=$(awk '{lenq+=($4-$3); lent+=($9-$8); n+=1} END {printf "%.0f\n", lenq + lent}' $PAF)

    QUERY_COVERED=$(join \
      <(cut -f 1,2 $PAF | sort -V | uniq) \
      <(awk -v OFS='\t' '{print($1,$3,$4)}' $PAF | bedtools sort | bedtools merge | awk 'BEGIN { FS="\t"; OFS="\t" } {
          seq_name = $1
          start = $2
          end = $3
          sum[seq_name] += (end - start)
      }
      END {
          for (seq_name in sum) {
              print seq_name, sum[seq_name]
          }
      }' | sort -V | uniq) | grep '^0' -v | awk '{ sum_len += $2; sum_covered += $3 } END { print(sum_covered / sum_len)}')

    TARGET_COVERED=$(join \
      <(cut -f 6,7 $PAF | sort -V | uniq) \
      <(awk -v OFS='\t' '{print($6,$8,$9)}' $PAF | bedtools sort | bedtools merge | awk 'BEGIN { FS="\t"; OFS="\t" } {
          seq_name = $1
          start = $2
          end = $3
          sum[seq_name] += (end - start)
      }
      END {
          for (seq_name in sum) {
              print seq_name, sum[seq_name]
          }
      }' | sort -V | uniq) | grep '^0' -v | awk '{ sum_len += $2; sum_covered += $3 } END { print(sum_covered / sum_len)}')

    echo $TOOL $NAME $AVG_ALIGNMENT_LEN $QUERY_COVERED $TARGET_COVERED
  done
done | tr ' ' '\t' > $DIR_BASE/alignments/$PANGENOME/$PANGENOME.contigs+coverage+totlen.tsv
```

Call variants from PAF files:

```shell
ls $DIR_BASE/alignments/$PANGENOME | grep '.tsv' -v | while read TOOL; do
  ls $DIR_BASE/alignments/$PANGENOME/$TOOL/*/*paf | while read PAF; do
    NAME=$(basename $PAF .paf)
    NAME12=${NAME%%.minimap2.*}
    NAME12=${NAME12%%.wfmash.*}
    IFS="-vs-" read -r NAME1 NAME2 <<< $NAME12
    NAME2="${NAME2#vs-}" # trim prefix
    echo $NAME $NAME12 $NAME1 $NAME2

    FASTA1=$DIR_BASE/assemblies/$PANGENOME/$NAME1.fasta
    FASTA2=$DIR_BASE/assemblies/$PANGENOME/$NAME2.fasta
    
    DIR_OUTPUT=$DIR_BASE/variants/$PANGENOME/$TOOL/$NAME1-vs-$NAME2
    mkdir -p $DIR_OUTPUT
    cd $DIR_OUTPUT

    sbatch -c 8 -p allnodes --job-name "paf2vcf-$NAME" --wrap "bash $DIR_BASE/scripts/paf2vcf.sh $PAF $NAME $FASTA1 $FASTA2 $PANDEPTH $DIR_OUTPUT"
  done
done
```

Call variants from HiFi reads:

```shell
ls $DIR_BASE/assemblies/$PANGENOME/*fasta | while read FASTA; do
  NAME=$(basename $FASTA .fasta)

  FASTQ=$DIR_BASE/hifi/$NAME.fq.gz
  if [[ -f "$FASTQ" ]]; then
    echo $NAME

    DIR_OUTPUT=$DIR_BASE/hifi-vs-assembly/$PANGENOME/
    mkdir -p $DIR_OUTPUT
    cd $DIR_OUTPUT
    sbatch -c 48 -p allnodes --job-name "hifi-vs-$NAME" --wrap "hostname; cd /scratch; minimap2 -ax map-hifi -t 42 $FASTA $FASTQ | samtools sort -m 8G -@ 6 -T /scratch/$NAME.HiFi.tmp -O bam -o $NAME.hifi.mm2.bam; samtools faidx $NAME.hifi.mm2.bam; mv $NAME.hifi.mm2.bam* $DIR_OUTPUT"
  fi
done

# On lambda01
cd ~
BIN_VERSION="1.6.1"
ls $DIR_BASE/assemblies/$PANGENOME/*fasta | while read FASTA; do
  NAME=$(basename $FASTA .fasta)

  FASTQ=$DIR_BASE/hifi/$NAME.fq.gz
  if [[ -f "$FASTQ" ]]; then
    echo $NAME

    DIR_OUTPUT=$DIR_BASE/deepvariant-vs-assembly/$PANGENOME/
    mkdir -p $DIR_OUTPUT
    cd $DIR_OUTPUT

    apptainer run \
        -B /usr/lib/locale/:/usr/lib/locale/ \
        -B "${DIR_BASE}":"/input" \
        -B "${DIR_OUTPUT}":"/output" \
        -B "/scratch":"/temporary" \
        deepvariant_"${BIN_VERSION}".sif \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=PACBIO \
        --ref=/input/assemblies/$PANGENOME/$NAME.fasta \
        --reads=/input/hifi-vs-assembly/$PANGENOME/$NAME.hifi.mm2.bam \
        --output_vcf=/output/$NAME.hifi.mm2.dv.vcf.gz \
        --output_gvcf=/output/$NAME.hifi.mm2.dv.g.vcf.gz \
        --intermediate_results_dir /temporary/$NAME2.DeepVariant.tmp \
        --num_shards=48

    rm /scratch/$NAME2.DeepVariant.tmp -rf
  fi
done
```


### All-vs-one

Alignments:

```shell
#########################
PANGENOME=athaliana
REF_FASTA=$DIR_BASE/assemblies/$PANGENOME/GCA_028009825.2.fasta

PANGENOME=tomato
REF_FASTA=$DIR_BASE/assemblies/$PANGENOME/SL5.0.fasta
#########################
REFNAME=$(basename $REF_FASTA .fasta)

p=95
s=10k
c=20k
k=19

ls $DIR_BASE/assemblies/$PANGENOME/*fasta | grep -v $REFNAME | while read FASTA; do
  NAME=$(basename $FASTA .fasta)

  echo $PANGENOME $NAME

  # # nucmer TO DO #########################################################################################################################################
  # DIR_OUTPUT_NC=$DIR_BASE/alignments/$PANGENOME/nucmer/$NAME
  # mkdir -p $DIR_OUTPUT_NC
  # cd $DIR_OUTPUT_NC
  # nucmer --maxmatch -c 500 -b 500 -l 100 -p $DIR_OUTPUT_NC/$NAME-vs-$REFNAME.nucmer.b500c500l100 -t $THREADS $DIR_BASE/assemblies/$PANGENOME/$REFNAME.fa $FASTA

  # minimap2
  DIR_OUTPUT_MM2=$DIR_BASE/alignments/$PANGENOME/minimap2/$NAME
  mkdir -p $DIR_OUTPUT_MM2
  cd $DIR_OUTPUT_MM2
  sbatch -c 48 -p allnodes --job-name "$NAME-minimap2" --wrap "hostname; cd /scratch; \time -v minimap2 -x asm20 -c --eqx --secondary=no -t $THREADS $REF_FASTA $FASTA > /scratch/$NAME.minimap2.asm20.paf; mv /scratch/$NAME.minimap2.asm20.paf $DIR_OUTPUT_MM2/"

  # # AnchorWave
  # DIR_OUTPUT_AW=$DIR_BASE/alignments/$PANGENOME/anchorwave/$NAME
  # mkdir -p $DIR_OUTPUT_AW
  # cd $DIR_OUTPUT_AW
  # sbatch -c 48 -p workers --job-name "$NAME-anchorwave" --wrap "hostname; mkdir -p /scratch/$NAME.anchorwave; cd /scratch/$NAME.anchorwave; \time -v anchorwave gff2seq -i $DIR_BASE/annotation/$PANGENOME/$REFNAME.gff3 -r $REF_FASTA -o $REFNAME.filter.cds.fa; \time -v minimap2 -x splice -t $THREADS -k 12 -a -p 0.4 -N 20 $REF_FASTA $REFNAME.filter.cds.fa > $REFNAME.cds.sam; \time -v minimap2 -x splice -t $THREADS -k 12 -a -p 0.4 -N 20 $FASTA $REFNAME.filter.cds.fa > $NAME.cds.sam; \time -v anchorwave genoAli -IV -t $THREADS -i $DIR_BASE/annotation/$PANGENOME/$REFNAME.gff3 -as $REFNAME.filter.cds.fa -r $REF_FASTA -a $NAME.cds.sam -ar $REFNAME.cds.sam -s $FASTA -n $NAME.AnchorWave.default.anchors -o $NAME.AnchorWave.default.maf -f $NAME.AnchorWave.default.m.maf; mv /scratch/$NAME.anchorwave/* $DIR_OUTPUT_AW/; rm /scratch/$NAME.anchorwave -rf"

  # wfmash
  DIR_OUTPUT_WF=$DIR_BASE/alignments/$PANGENOME/wfmash/$NAME
  mkdir -p $DIR_OUTPUT_WF
  cd $DIR_OUTPUT_WF
  sbatch -c 48 -p allnodes --job-name "$NAME-wfmash" --wrap "hostname; cd /scratch; \time -v $WFMASH -p $p -s $s -c $c -k $k -t $THREADS --hg-filter-ani-diff 30 $REF_FASTA $FASTA > $NAME.wfmash.p${p}s${s}c${c}k${k}.paf; mv /scratch/$NAME.wfmash.wfmash.p${p}s${s}c${c}k${k}.paf $DIR_OUTPUT_WF/"
done
```








# sv_genotype
mkdir -p results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_minimap2_asm20_tmp
cuteSV results/02.reads_mapping/Slycopersicum/TS421.hifi.sorted.bam fasta/Slycopersicum/SL5.fa.gz results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.cutesv.support.vcf results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_minimap2_asm20_tmp -s 1 --genotype -mi 500 -md 500 --min_mapq 0 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --min_size 50 -Ivcf results/01.wga/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.all_variants.filter.svs.vcf.gz
bgzip -@ 1 results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_minimap2_asm20.cutesv.support.vcf
rm -rf results/03.variants/SVs/Slycopersicum/Slycopersicum_TS421_minimap2_asm20_tmp




mkdir -p /lizardfs/guarracino/wfmash-paper/images.athaliana.p70-k1M
$WFMASH_PLOTS /lizardfs/guarracino/wfmash-paper/assemblies/athaliana/GCA_028009825.2.fasta /lizardfs/guarracino/wfmash-paper/assemblies/athaliana/GCA_903064275.1.fasta -t 40 -p 70 -l 1000k --wfplot-max-size 10000 -u /lizardfs/guarracino/wfmash-paper/images.athaliana.p70-k1M/athaliana. > /dev/null

mkdir -p /lizardfs/guarracino/wfmash-paper/images.acros
$WFMASH_PLOTS /lizardfs/guarracino/robertsonian_translocation/mappings/verkko/GM03786/GM03786.haplotype1-0000026.fa.gz /lizardfs/guarracino/robertsonian_translocation/NGenomeSyn/GM03786/chm13v2.0.13+14c.fa.gz -t 40 --wfplot-max-size 1000 -u /lizardfs/guarracino/wfmash-paper/images.acros/chm13-vs-rob. > /dev/null

mkdir -p /lizardfs/guarracino/wfmash-paper/images.scerevisiae.p70-N
$WFMASH_PLOTS /lizardfs/guarracino/wfmash-paper/assemblies/scerevisiae/scerevisiae7.fa.gz -t 40 -p 70 -Y '#' -N -u /lizardfs/guarracino/wfmash-paper/images.scerevisiae.p70-N/scerevisiae7. > /dev/null

mkdir -p /lizardfs/guarracino/wfmash-paper/images.primates-hsa6.p70-k5M
$WFMASH_PLOTS /lizardfs/guarracino/pggb-paper/assemblies/primates/primates16.hsa6.fa.gz -t 40 -p 70 -c 100k -l 5000k -u /lizardfs/guarracino/wfmash-paper/images.primates-hsa6.p70-k10M/primates-hsa6. > /dev/null