#!/bin/bash

PAF=$1
BED=$2
OUTDIR=$3

# Count total number of targets and genes
TOTAL_TARGETS=$(cut -f 6 "$PAF" | cut -f 1,2 -d '#' | sort | uniq | wc -l)
TOTAL_GENES=$(cut -f 4 "$BED" | sort | uniq | wc -l)

# Initialize counters
TARGET_COUNTER=0
GENE_COUNTER=0

cd /scratch
cut -f 6 "$PAF" | cut -f 1,2 -d '#' | sort | uniq | while read TARGET; do
    ((TARGET_COUNTER++))
    echo "Processing target $TARGET_COUNTER of $TOTAL_TARGETS: $TARGET" >&2

    cut -f 4 "$BED" | sort | uniq | while read GENE; do
        ((GENE_COUNTER++))
        
        # Calculate and print progress
        PROGRESS=$(( (TARGET_COUNTER - 1) * 100 / TOTAL_TARGETS + (GENE_COUNTER * 100 / TOTAL_GENES / TOTAL_TARGETS) ))
        echo "Progress: $PROGRESS% | Target: $TARGET_COUNTER/$TOTAL_TARGETS | Gene: $GENE_COUNTER/$TOTAL_GENES" >&2
        
        awk -v gene="$GENE" '$4 == gene' "$BED" | bedtools sort > temp.bed
        bedtools subtract \
            -a temp.bed \
            -b <(impg -p "$PAF" -b <(grep "^$TARGET" temp.bed | bedtools sort))
        rm temp.bed
    done
    
    # Reset gene counter for next target
    GENE_COUNTER=0
done > "$OUTDIR/$(basename "$PAF" .paf)-vs-$(basename "$BED" .bed).uncovered-regions.bed"

echo "Processing complete. Output file: $OUTDIR/$(basename "$PAF" .paf)-vs-$(basename "$BED" .bed).uncovered-regions.bed" >&2
