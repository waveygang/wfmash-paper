# Download the T2T assembly
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz

gunzip chm13.draft_v1.1.fasta.gz
samtools faidx chm13.draft_v1.1.fasta
samtools faidx chm13.draft_v1.1.fasta $(grep chr22 chm13.draft_v1.1.fasta.fai | cut -f 1) > chr22.fa && samtools faidx chr22.fa

run_splitfa=~/git/fork/splitfa/target/release/splitfa # https://github.com/AndreaGuarracino/splitfa.git, random_segment_length branch


threads=16
fasta=chr22.fa
fasta_name=$(basename $fasta .fa)


# Create simulated reads
for len in 200 500 1000 2000 5000 10000 20000 50000 100000 500000 1000000; do
  ~/git/PBSIM-PacBio-Simulator/src/pbsim --depth 2 --data-type CLR --accuracy-mean 0.99 --accuracy-min 0.99 --length-min $len --length-mean $len --length-max $len --model_qc ~/git/PBSIM-PacBio-Simulator/data/model_qc_clr $fasta
  ~/git/minimap2/misc/paftools.js pbsim2fq $fasta sd_0001.maf | sed 's/!>/!/g' | bgzip -@ $threads -c > pbsim.${fasta_name}.${len}bp.fa.gz && samtools faidx pbsim.${fasta_name}.${len}bp.fa.gz
  rm sd_0001.*
done

# Mapping and alignment
for len in 200 500 1000 2000 5000 10000 20000 50000 100000 500000 1000000; do
  s=$(echo "$len/10" | bc)
  s=$(( $s > 100 ? $s : 100 ))
  s=$(( $s < 50000 ? $s : 50000 ))

  \time -v minimap2 $fasta pbsim.${fasta_name}.${len}bp.fa.gz -t $threads -c -x asm5 > pbsim.${fasta_name}.${len}bp.mp2.paf
done

for len in 200 500 1000 2000 5000 10000 20000 50000 100000 500000 1000000; do
  s=$(echo "$len/10" | bc)
  s=$(( $s > 100 ? $s : 100 ))
  s=$(( $s < 50000 ? $s : 50000 ))

  \time -v wfmash $fasta pbsim.${fasta_name}.${len}bp.fa.gz -t $threads -s $s -l 0 -p 98 -n 5 -N -m > pbsim.${fasta_name}.${len}bp.N.wfm.approx.paf
done
for len in 200 500 1000 2000 5000 10000 20000 50000 100000 500000 1000000; do
  s=$(echo "$len/10" | bc)
  s=$(( $s > 100 ? $s : 100 ))
  s=$(( $s < 50000 ? $s : 50000 ))

  \time -v wfmash $fasta pbsim.${fasta_name}.${len}bp.fa.gz -t $threads -s $s -l 0 -p 98 -n 5 -N -i pbsim.${fasta_name}.${len}bp.N.wfm.approx.paf > pbsim.${fasta_name}.${len}bp.N.wfm.paf
done

# Evaluation (https://github.com/lh3/minimap2/blob/master/misc/README.md#mapeval)
echo read.len mapq.threshold num.reads num.wrong.reads ratio | tr ' ' '\t' > mapeval.tsv
for len in 200 500 1000 2000 5000 10000 20000 50000 100000 500000 1000000; do
  echo "read length: $len"
  #echo "  minimap2"
  #~/git/minimap2/misc/paftools.js mapeval -r 0.6 pbsim.${fasta_name}.${len}bp.mp2.paf
  #echo "  wfmash"
  #~/git/minimap2/misc/paftools.js mapeval -r 0.6 pbsim.${fasta_name}.${len}bp.N.wfm.paf
  ~/git/minimap2/misc/paftools.js mapeval -r 0.6 pbsim.${fasta_name}.${len}bp.N.wfm.paf | awk -v OFS='\t' -v len=$len '{print(len,$2,$3,$4,$4/$3)}' >> mapeval.tsv
  echo ""
  echo ""
done







$run_splitfa $ref -l $len-$len -s 1 2> $ref.l${len}.paf | bgzip -@ $threads -c > $ref.l${len}.fa.gz; samtools faidx $ref.l${len}.fa.gz


wfmash $ref $ref.l${len}.fa.gz -t $threads -s 500 -l 0 -p 99 -N -n 3 | sort -k 1 > $ref.l${len}.s500.l0.p99.n3.paf
