#!/bin/bash

# Input arguments
PATH_PAF=$1         # It assumes a WFMASH PAF. TODO: to make it general, to support also other sequence aligners
PATH_BED=$2         # It assumes that the BED file contains the feature names in the 4th column and strand information in the 6th column
MAX_INDEL_SIZE=$3   # Max indel size to allow in the alignment covering the feature intervals (-1 to disable such a filter)
OUTPUT_PREFIX=$4    # Prefix for the output files (/path/where/to/write/the/output/prefix)
DIR_TEMP_STUFF=$5   # Directory to store temporary files

# Check for file existence and non-emptiness
for file in "$PATH_PAF" "$PATH_BED"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: $file does not exist"
        exit 1
    elif [ ! -s "$file" ]; then
        echo "ERROR: $file is empty"
        exit 1
    fi
done

# Extract the directory path from the OUTPUT_PREFIX
DIR_OUTPUT=$(dirname "$OUTPUT_PREFIX")

# Check if the directories exists
for directory in "$DIR_TEMP_STUFF" "$DIR_OUTPUT"; do
    if [ ! -d "$directory" ]; then
        echo "ERROR: $directory does not exist"
        exit 1
    fi
done

# Get script's directory
DIR_SCRIPT=$( cd -- "$(dirname -- "$(readlink -f "${BASH_SOURCE[0]}" )" )" &> /dev/null && pwd )

# Extract unique query and target names
QUERY_TARGET_PAIRS=$(cut -f 1,6 "$PATH_PAF" | sort -u)

# Calculate total number of iterations
TOTAL=$(echo "$QUERY_TARGET_PAIRS" | wc -l)

# Initialize counter
COUNT=0

# Function to update progress bar
function update_progress {
    local PERCENT=$((100 * COUNT / TOTAL))
    printf "\r$COUNT / $TOTAL - $PERCENT%%"
}

# Feature-level report
PATH_FEATURE_TSV=$OUTPUT_PREFIX.report.features.tsv

# Genome-level report
PATH_GENOME_TSV=$OUTPUT_PREFIX.report.genomes.tsv

# Put the header
python3 $DIR_SCRIPT/feature_level_report.py > "$PATH_FEATURE_TSV"

# Check present pairs of query and target names
echo "$QUERY_TARGET_PAIRS" | while IFS=$'\t' read -r query target; do
    #echo "$query - $target"

    NAME=$(basename $PATH_PAF .paf).$query-$target

    # Extract the alignments for the current query-target pair
    grep -P "^$query\t" "$PATH_PAF" | grep -P "\t$target\t" > "$DIR_TEMP_STUFF/$NAME.paf"

    # Extract features to check for the current query-target pair
    grep -P "^$query\t" "$PATH_BED" > "$DIR_TEMP_STUFF/$NAME.query.bed"
    grep -P "^$target\t" "$PATH_BED" > "$DIR_TEMP_STUFF/$NAME.target.bed"

    # Take alignments covering the single-copy BUSCO genes
    bedtools intersect -a <(awk -v OFS='\t' '{print($1,$3,$4,$0)}' "$DIR_TEMP_STUFF/$NAME.paf") -b "$DIR_TEMP_STUFF/$NAME.query.bed"  -wa -wb | cut -f 4- > "$DIR_TEMP_STUFF/$NAME.query.paf"
    bedtools intersect -a <(awk -v OFS='\t' '{print($6,$8,$9,$0)}' "$DIR_TEMP_STUFF/$NAME.paf") -b "$DIR_TEMP_STUFF/$NAME.target.bed" -wa -wb | cut -f 4- > "$DIR_TEMP_STUFF/$NAME.target.paf"

    # TODO: what happens if the gene is covered by more mappings/alignments in one genome than in the other
    # Join the two files by the alignment columns + the gene name
    join -1 1 -2 1 \
        <(awk -v OFS='\t' '{concat=$1; for(i=2;i<=20;i++) concat=concat "___" $i; concat=concat "___" $24; print(concat,$21,$22,$23,$25,$26,$27)}' "$DIR_TEMP_STUFF/$NAME.query.paf"  | sort -T /scratch) \
        <(awk -v OFS='\t' '{concat=$1; for(i=2;i<=20;i++) concat=concat "___" $i; concat=concat "___" $24; print(concat,$21,$22,$23,$25,$26,$27)}' "$DIR_TEMP_STUFF/$NAME.target.paf" | sort -T /scratch) |
        # Deconcatenate the joined column back into original columns
        awk '{
            # Split the first field into array `a` using the "|" delimiter
            n=split($1, a, "___");
            # Reconstruct the original columns
            for(i=1; i<=20; i++) printf("%s\t", a[i]);
            printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t", $2, $3, $4, a[21], $5, $6, $7);
            printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $8, $9, $10, a[21], $11, $12, $13);
        }' | cut -f 1-12,20- | pigz -9 > "$DIR_TEMP_STUFF/$NAME.joined.paf.gz"

    # Generate gene-level report
    python3 $DIR_SCRIPT/feature_level_report.py -i "$DIR_TEMP_STUFF/$NAME.joined.paf.gz" -m $MAX_INDEL_SIZE | sed '1d' >> $PATH_FEATURE_TSV

    # Remove the temporary files
    rm "$DIR_TEMP_STUFF/$NAME.paf" "$DIR_TEMP_STUFF/$NAME.query.bed" "$DIR_TEMP_STUFF/$NAME.target.bed" "$DIR_TEMP_STUFF/$NAME.query.paf" "$DIR_TEMP_STUFF/$NAME.target.paf" "$DIR_TEMP_STUFF/$NAME.joined.paf.gz"

    ((COUNT++))
    update_progress
done

# generate genome-level report
echo "query num.features.in.query target  num.feature.in.target num.features.in.common aligned.bases not.aligned.in.query.bp not.aligned.in.target.bp" | tr ' ' '\t' > $PATH_GENOME_TSV

# Check present pairs of query and target names
echo "$QUERY_TARGET_PAIRS" | while IFS=$'\t' read -r query target; do
    NUM_FEATURES_QUERY=$(grep -P "^$query\t" "$PATH_BED" | wc -l)
    NUM_FEATURES_TARGET=$(grep -P "^$target\t" "$PATH_BED" | wc -l)

    awk -v q="$query" -v t="$target" '$2 == q && $6 == t' "$PATH_FEATURE_TSV" | awk -v OFS='\t' -v nq="$NUM_FEATURES_QUERY" -v nt="$NUM_FEATURES_TARGET" '{
            key = $2 FS $6; # Combine query and target as a unique key
            count[key]++;
            aligned[key] += $9;
            not_aligned_query[key] += $10;
            not_aligned_target[key] += $11;
        }
        END {
            for (k in count) {
                split(k, keys, FS); # Split the key back into query and target
                print keys[1], nq, keys[2], nt, count[k], aligned[k], not_aligned_query[k], not_aligned_target[k];
            }
        }' >> $PATH_GENOME_TSV
done
