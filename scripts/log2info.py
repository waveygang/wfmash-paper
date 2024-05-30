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
            if 'asm20' in line:
                tool = 'minimap2.xasm20.c.eqx.secondary-no'
            elif 'asm10' in line:
                tool = 'minimap2.xasm10.c.eqx.secondary-no'
            elif 'asm5' in line:
                tool = 'minimap2.xasm5.c.eqx.secondary-no'
            else:
                reject_result = True
        elif 'wfmash' in line:
            tool = 'wfmash'
            param_tuples = re.findall(r"-(p|s|l|c|k)\s+(\S+)", line)
            wfmash_params = ".".join([f"{p}{v}" for p, v in param_tuples])
            hg_filter_ani_diff = re.search(r"--hg-filter-ani-diff\s+(\d+)", line)
            if hg_filter_ani_diff:
                wfmash_params += f".hg{hg_filter_ani_diff.group(1)}"
            one_to_one_filter = re.search(r"--one-to-one", line)
            if one_to_one_filter:
                wfmash_params += f".yes1to1"
            else:
                wfmash_params += f".no1to1"
            w_param = re.search(r"-w\s+(\d+)", line)
            if w_param:
                wfmash_params += f".w{w_param.group(1)}"
            else:
                wfmash_params += ".wauto"
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
