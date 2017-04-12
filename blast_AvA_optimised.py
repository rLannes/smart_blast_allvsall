#!usr/local/env python3


'''
  blast_AvA_optimised.py this script will make a blast all vs all comparison limiting the number of
   alignement done using the asumption
 that at posteriori the results are sort by mutual coverage. Usefull for homology detection

Copyright (C) 2017  Lannes Romain France Paris

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''



# this script will make a blast all vs all comparison limiting the number of alignement done using the asumption
# that at posteriori the results are sort by mutual coverage.
# using this see that some alignement are unecessary.


# algorithme:

# first compute nb_seq _ fasta_size
# first make a fasta by sequence size
# we couls also index the fasta by sequence size and get the size fasta each time wee ned it ( I got a script doing almost this)

# second compute and launch the blast (doing in the same time the required databases, erasing same after use)
# It will be very important to erase file and database that are no longer necessary
# we may use queu instead of pool for mp depending at which point we want to optimise it

# cat result at soon they arrive in a results file then delete them

# end

import os
import argparse
import subprocess
import utils.utils_fonction
import time


parser = argparse.ArgumentParser(description='')

parser.add_argument('-log', help='log file')

parser.add_argument('-input', help='input', required=True)

parser.add_argument('-wd', help='working_dir', required=True)

parser.add_argument('-cov', help='coverage threshold default 0.8', type=float, default=0.8)

parser.add_argument('-th', help='number_of_thread', type=int, required=True)

parser.add_argument('-mkdb', help='makeblastdb path', default='makeblastdb')

parser.add_argument('-blastp', help='blastp  path', default='blastp')

parser.add_argument('-output', help='output', required=True)

parser.add_argument('-blast_header', help='the blast output format default =\
 qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen', nargs='*',\
  default = "6 qseqid sseqid evalue pident bitscore qstart qend sstart send qlen slen")

parser.add_argument('-blast_other_arg', help='other arguments for blast default =\
 -seg yes -soft_masking true -max_hsps 1', nargs='*', default="-seg yes -soft_masking true -max_hsps 1")

args = parser.parse_args()

wd_dir = args.wd
os.mkdir(wd_dir)

split_fasta = wd_dir + '/spit_faa'
os.mkdir(split_fasta)


start_time = time.time()
dico_size_distribution = utils.utils_fonction.split_fasta_by_size_and_return_size_distribution(fasta_file=args.input,
                                                                                                out_dir=split_fasta)

print('making sub fasta files takes {}'.format(round(time.time() - start_time, 2)))

liste_size_sorted = sorted(list(dico_size_distribution.keys()))

dico_which_file_against = utils.utils_fonction.get_dico_which_file_for_which_id(dico_size_distribution=dico_size_distribution,
                                                                           out_dir=split_fasta, cov=args.cov)


liste_to_working_dir = utils.utils_fonction.launch_mp(dico_which_file_against, nb_thread=args.th, working_directory=args.wd,
                                                    makeblast_db=args.mkdb, blastp_=args.blastp,
                                                    output_format_arg=args.blast_header,
                                                    other_args=args.blast_other_arg, final_output_=args.output, cov=args.cov)


child = subprocess.Popen("rm -rf {}".format(' '.join(liste_to_working_dir)), shell=True)
child.wait()

print('done =)')