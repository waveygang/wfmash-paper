#!/bin/bash

# Input arguments
PATH_PAF=$1         # It assumes a WFMASH PAF. TODO: to make it general, to support also other sequence aligners
PATH_BED=$2         # It assumes that the BED file contains the feature names in the 4th column and strand information in the 6th column
MAX_INDEL_SIZE=$3   # Max indel size to allow in the alignment covering the feature intervals (-1 to disable such a filter)
OUTPUT_PREFIX=$4    # Prefix for the output files
NUM_THREADS=$5      # Prefix for the output files
DIR_TEMP_STUFF=$6   # Directory to store temporary files

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

# Calculate total number of pairs
NUM_PAIRS=$(echo "$QUERY_TARGET_PAIRS" | wc -l)

# Feature-level report
PATH_FEATURE_TSV=$OUTPUT_PREFIX.report.features.tsv

# Genome-level report
PATH_GENOME_TSV=$OUTPUT_PREFIX.report.genomes.tsv

process_pair() {
    local query="$1"
    local target="$2"
    local DIR_TEMP_STUFF="$3"
    local PATH_PAF="$4"
    local PATH_BED="$5"
    local MAX_INDEL_SIZE="$6"
    local NUM_THREADS="$7"
    local OUTPUT_PREFIX="$8"
    local PATH_FEATURE_TSV="$9"

    local NAME=$(basename "$PATH_PAF" .paf)."$query"-"$target"

    # Extract the alignments for the current query-target pair
    grep -P "^$query\t" "$PATH_PAF" | grep -P "\t$target\t" > "$DIR_TEMP_STUFF/$NAME.paf"

    # Extract features to check for the current query-target pair
    grep -P "^$query\t" "$PATH_BED" > "$DIR_TEMP_STUFF/$NAME.query.bed"
    grep -P "^$target\t" "$PATH_BED" > "$DIR_TEMP_STUFF/$NAME.target.bed"

    # Take alignments covering the single-copy BUSCO genes and process further
    bedtools intersect -a <(awk -v OFS='\t' '{print($1,$3,$4,$0)}' "$DIR_TEMP_STUFF/$NAME.paf") -b "$DIR_TEMP_STUFF/$NAME.query.bed"  -wa -wb | cut -f 4- > "$DIR_TEMP_STUFF/$NAME.query.paf"
    bedtools intersect -a <(awk -v OFS='\t' '{print($6,$8,$9,$0)}' "$DIR_TEMP_STUFF/$NAME.paf") -b "$DIR_TEMP_STUFF/$NAME.target.bed" -wa -wb | cut -f 4- > "$DIR_TEMP_STUFF/$NAME.target.paf"

    # TODO: what happens if the gene is covered by more mappings/alignments in one genome than in the other
    # Join the two files by the alignment columns + the gene name
    join -1 1 -2 1 \
        <(awk -v OFS='\t' '{concat=$1; for(i=2;i<=20;i++) concat=concat "___" $i; concat=concat "___" $24; print(concat,$21,$22,$23,$25,$26,$27)}' "$DIR_TEMP_STUFF/$NAME.query.paf"  | sort -T $DIR_TEMP_STUFF) \
        <(awk -v OFS='\t' '{concat=$1; for(i=2;i<=20;i++) concat=concat "___" $i; concat=concat "___" $24; print(concat,$21,$22,$23,$25,$26,$27)}' "$DIR_TEMP_STUFF/$NAME.target.paf" | sort -T $DIR_TEMP_STUFF) |
        # Deconcatenate the joined column back into original columns
        awk '{
            # Split the first field into array `a` using the "|" delimiter
            n=split($1, a, "___");
            # Reconstruct the original columns
            for(i=1; i<=20; i++) printf("%s\t", a[i]);
            printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t", $2, $3, $4, a[21], $5, $6, $7);
            printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $8, $9, $10, a[21], $11, $12, $13);
        }' | cut -f 1-12,20- | pigz > "$DIR_TEMP_STUFF/$NAME.joined.paf.gz"

    # Append the feature-level report data to the final TSV
    feature_level_report -i "$DIR_TEMP_STUFF/$NAME.joined.paf.gz" -m "$MAX_INDEL_SIZE" | sed '1d'

    # Cleanup
    rm "$DIR_TEMP_STUFF/$NAME.paf" "$DIR_TEMP_STUFF/$NAME.query.bed" "$DIR_TEMP_STUFF/$NAME.target.bed" "$DIR_TEMP_STUFF/$NAME.query.paf" "$DIR_TEMP_STUFF/$NAME.target.paf" "$DIR_TEMP_STUFF/$NAME.joined.paf.gz"
}
export -f process_pair

# Put the header into the feature-level report file
feature_level_report > "$PATH_FEATURE_TSV"

# Execute processing in parallel
JOBLOG="$OUTPUT_PREFIX.parallel.log" # --progress creates issues with sbatch (sh: line 1: /dev/tty: No such device or address)
echo "Processing $NUM_PAIRS query-target pairs in parallel using $NUM_THREADS threads. The log is stored in $JOBLOG."

echo "$QUERY_TARGET_PAIRS" | parallel --no-notice --joblog "$JOBLOG" --colsep $'\t' -j "$NUM_THREADS" process_pair {1} {2} "$DIR_TEMP_STUFF" "$PATH_PAF" "$PATH_BED" "$MAX_INDEL_SIZE" "$NUM_THREADS" "$OUTPUT_PREFIX" "$PATH_FEATURE_TSV" >> "$PATH_FEATURE_TSV"

generate_genome_report() {
    local query="$1"
    local target="$2"
    local PATH_BED="$3"
    local PATH_FEATURE_TSV="$4"
    
    local NUM_FEATURES_QUERY=$(grep -P "^$query\t" "$PATH_BED" | wc -l)
    local NUM_FEATURES_TARGET=$(grep -P "^$target\t" "$PATH_BED" | wc -l)

    awk -v q="$query" -v t="$target" -v nq="$NUM_FEATURES_QUERY" -v nt="$NUM_FEATURES_TARGET" 'BEGIN {OFS="\t"} 
        $2 == q && $6 == t {
            key = $2 FS $6; # Combine query and target as a unique key
            count[key]++;
            aligned[key] += $9;
            not_aligned_in_query[key] += $10;
            not_aligned_in_target[key] += $11;
            indels_in_query[key] += $12;
            indels_in_target[key] += $13;
            ignored_bases_in_query[key] += $14;
            ignored_bases_in_target[key] += $15;
        }
        END {
            for (k in count) {
                split(k, keys, FS); # Split the key back into query and target
                print keys[1], nq, keys[2], nt, count[k], aligned[k], not_aligned_in_query[k], not_aligned_in_target[k], indels_in_query[k], indels_in_target[k], ignored_bases_in_query[k], ignored_bases_in_target[k];
            }
        }' "$PATH_FEATURE_TSV"
}
export -f generate_genome_report


# generate genome-level report
echo "query num.features.in.query target  num.feature.in.target num.features.in.common aligned.bases not.aligned.in.query.bp not.aligned.in.target.bp indels.in.query.bp indels.in.target.bp ignored.in.query.bp ignored.in.target.bp" | tr ' ' '\t' > "$PATH_GENOME_TSV"

echo "$QUERY_TARGET_PAIRS" | parallel --colsep $'\t' -j "$NUM_THREADS" generate_genome_report {1} {2} "$PATH_BED" "$PATH_FEATURE_TSV" >> "$PATH_GENOME_TSV"
