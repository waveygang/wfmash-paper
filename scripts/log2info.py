# Usage:
#    cat *.log | python3 scripts/log2info.py

import sys
import os
import re
def get_filename_without_extension(file_path):
    # Get the base name of the file (i.e., the name with the path removed)
    base_name = os.path.basename(file_path)
    # Split the name by '.' and discard the extension
    file_name_without_extension = os.path.splitext(base_name)[0]
    return file_name_without_extension

elapsed_wall_clock_time = ''
max_resident_set_size = ''

reject_result = False

for line in sys.stdin:
    if 'Command terminated by signal' in line:
        reject_result = True
    elif 'Command being timed' in line:
        if 'minimap2' in line:
            tool = 'minimap2.xasm20.c.eqx.secondary-no'
        elif 'wfmash' in line:
            tool = 'wfmash'
            param_tuples = re.findall(r"-(p|s|c|k)\s+(\S+)", line)
            wfmash_params = ".".join([f"{p}{v}" for p, v in param_tuples])
            hg_filter_ani_diff = re.search(r"--hg-filter-ani-diff\s+(\d+)", line)
            if hg_filter_ani_diff:
                wfmash_params += f".hg{hg_filter_ani_diff.group(1)}"
            tool = f"wfmash.{wfmash_params}"
        else:
            reject_result = True

        # Find all arguments that end with '.fasta'
        fasta1 = get_filename_without_extension(line.split()[-1])
        fasta2 = get_filename_without_extension(line.split()[-2])

        elapsed_wall_clock_time = ''
        max_resident_set_size = ''
    elif 'Elapsed (wall clock) time' in line:
        elapsed_wall_clock_time = line.strip().split('): ')[-1]

        if len(elapsed_wall_clock_time.split(':')) == 3:
            # hh:mm:ss
            hh, mm, ss = elapsed_wall_clock_time.split(':')

            elapsed_wall_clock_time = str(int(hh) * 3600 + int(mm) * 60 + int(ss))
        elif len(elapsed_wall_clock_time.split(':')) == 2:
            # mm:ss.xx
            mm, ss = elapsed_wall_clock_time.split(':')
            ss = ss.split('.')[0]

            elapsed_wall_clock_time = str(int(mm) * 60 + int(ss))
        else:
            # Not expected
            elapsed_wall_clock_time = ''

        max_resident_set_size = ''
    elif 'Maximum resident set size' in line:
        max_resident_set_size = line.strip().split('): ')[-1]

        # Convert in Mbytes
        max_resident_set_size = '{:.4f}'.format(float(max_resident_set_size)/1024)

        # Check if all information are available
        if not reject_result and fasta1 and fasta2 and elapsed_wall_clock_time and max_resident_set_size:
            print('\t'.join([tool, fasta1, fasta2, elapsed_wall_clock_time, max_resident_set_size]))
        
        fasta1 = ''
        fasta2 = ''
        elapsed_wall_clock_time = ''
        max_resident_set_size = ''

        reject_result = False
