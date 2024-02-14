import sys
import re
import argparse
import gzip

def count_aligned_bases(query_start, query_end, query_strand, target_start, target_end, cigar, feature_in_query_start, feature_in_query_end, feature_in_target_start, feature_in_target_end, max_indel_size):
    aligned_bases = 0
    not_aligned_bases_query = feature_in_query_end - feature_in_query_start
    not_aligned_bases_target = feature_in_target_end - feature_in_target_start
    query_rev = query_strand == '-'

    # Initialize counters for the current position within the query and target sequences.
    query_pos = query_end if query_rev else query_start
    target_pos = target_start

    # Iterate over CIGAR operations.
    for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        length = int(length)

        # Handle match/mismatch, which affects both query and target.
        if op in 'M=X':
            if query_rev:
                overlap_query = max(0, min(query_pos, feature_in_query_end) - max(query_pos - length, feature_in_query_start))
                query_pos += -length
            else:
                overlap_query = max(0, min(query_pos + length, feature_in_query_end) - max(query_pos, feature_in_query_start))
                query_pos += length

            overlap_target = max(0, min(target_pos + length, feature_in_target_end) - max(target_pos, feature_in_target_start))
            aligned_bases += min(overlap_query, overlap_target)
            target_pos += length

        # Handle deletion in the target (insertion in the query).
        elif op == 'D':
            target_pos += length
            if target_pos > feature_in_target_start and target_pos - length < feature_in_target_end:
                if max_indel_size is None or length > max_indel_size:
                    not_aligned_bases_target += min(length, feature_in_target_end - target_pos, target_pos - length - feature_in_target_start)

        # Handle insertion in the query (gap in the target).
        elif op == 'I':
            query_pos += -length if query_rev else length
            if query_pos > feature_in_query_start and query_pos - length < feature_in_query_end:
                if max_indel_size is None or length > max_indel_size:
                    not_aligned_bases_query += min(length, feature_in_query_end - query_pos, query_pos - length - feature_in_query_start)

    return aligned_bases, not_aligned_bases_query, not_aligned_bases_target

def main(input_file, max_indel_size=None):
    print('\t'.join(['feature.name', 'query', 'query.feature.start', 'query.feature.end', 'query.strand', 'target', 'target.feature.start', 'target.feature.end', 'aligned.bp', 'not.aligned.in.query.bp', 'not.aligned.in.target.bp']))
    
    if input_file:
        open_func = gzip.open if input_file.endswith('.gz') else open
        with open_func(input_file, 'rt') as file:
            for line in file:
                query_name, query_len, query_start, query_end, query_strand, target_name, target_len, target_start, target_end, _, _, _, cigar, query_name_2, feature_in_query_start, feature_in_query_end, feature_in_query_name, _, feature_in_query_strand, feature_in_query_class, target_name_2, feature_in_target_start, feature_in_target_end, feature_in_target_name, _, feature_in_target_strand, feature_in_target_class = line.strip().split('\t')

                if query_name != query_name_2 or target_name != target_name_2 or feature_in_query_name != feature_in_target_name:
                    print('ERROR: query, target, and/or feature name do not match.')
                    sys.exit(1)
                if (feature_in_query_strand != feature_in_target_strand) and query_strand == '+':
                    # If the features are on different strands, the query should be reversed in order to align them
                    print('ERROR: <<<<<<<<TO CONFIRM>>>>>>>> the feature is on different strands in query and target, but query and target are in the same orientation!')
                    sys.exit(1)

                cigar = cigar.split('cg:Z:')[-1]

                aligned_bases, not_aligned_bases_query, not_aligned_bases_target = count_aligned_bases(
                    int(query_start), int(query_end), query_strand, int(target_start), int(target_end), cigar, int(feature_in_query_start), int(feature_in_query_end), int(feature_in_target_start), int(feature_in_target_end), max_indel_size
                )

                print('\t'.join([feature_in_query_name, query_name, feature_in_query_start, feature_in_query_end, query_strand, target_name, feature_in_target_start, feature_in_target_end, str(aligned_bases), str(not_aligned_bases_query - aligned_bases), str(not_aligned_bases_target - aligned_bases)]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process alignments for feature level report.")
    parser.add_argument("-i", "--input-alg-tsv", type=str, help="Input alignments, can be gzipped")
    parser.add_argument("-m", "--max-indel-size", type=int, help="Maximum size of indels to consider in feature intervals (optional)", default=None)
    args = parser.parse_args()

    main(args.input_alg_tsv, args.max_indel_size)
