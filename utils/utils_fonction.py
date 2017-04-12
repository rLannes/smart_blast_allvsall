'''
<contains utils fonction for blast_AvA_optimised.py>
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






import subprocess
import os
from multiprocessing import Pool
import time



def time_reporteur_with_return(fonction):
	"""This decorator report the time a function amde to execute"""
	
	def wrapper(*args, **kwargs):
		start_time = time.time()
		execution = fonction(*args, **kwargs)
		elapsed_time = time.time() - start_time
		return execution
		print(fonction.__name__, ' :', elapsed_time, ' second')
	
	return wrapper


def time_reporteur_no_return(fonction):
	"""This decorator report the time a function amde to execute"""
	
	def wrapper(*args, **kwargs):
		start_time = time.time()
		execution = fonction(*args, **kwargs)
		elapsed_time = time.time() - start_time
		print(fonction.__name__, ' :', elapsed_time, ' second')
	
	return wrapper


def print_log(log_file, string):
	
	with open(log_file, 'a') as log:
		log.write(str(string) + '\n')


def get_seq_one_by_one(file_):
	"""Generator return prompt sequence for each sequence"""

	sequence = ''
	prompt = ''

	for line in file_:
		
		if line.startswith('>'):
			
			if sequence:
				yield [prompt, sequence]
				
			sequence = ''
			prompt = line.strip()
			
		else:
			sequence += line.strip()

	yield [prompt, sequence]


def chunck_sequence(sequence):
	
	inc = 80
	cpt = 0
	seq = ""
	while cpt <= len(sequence):
		seq += sequence[cpt:cpt + inc] + '\n'
		cpt += inc
	return seq.strip()

@time_reporteur_with_return
def split_fasta_by_size_and_return_size_distribution(fasta_file, out_dir):
	"""
	the fasta files will be not modify, this fonction will output in an already created directory
	fasta sequence by size -> mutliple file with a number corresponding to the size of the sequence inside
	"""
	
	dico_size_distribution = {}
	
	with open(fasta_file) as faa_file:
		
		for prompt, sequence in get_seq_one_by_one(file_=faa_file):
			
			len_sequence = len(sequence)
			dico_size_distribution[len_sequence] = dico_size_distribution.get(len_sequence, 0) + 1
			
			with open(out_dir + '/{}.faa'.format(len_sequence), 'a') as fasta_len_file:
				fasta_len_file.write('{}\n{}\n'.format(prompt, chunck_sequence(sequence)))
	
	return dico_size_distribution


def round_up(value):
	
	if int(value) != round(value):
		return int(value) + 1
	
	else:
		return int(value)


@time_reporteur_with_return
def get_dico_which_file_for_which_id(dico_size_distribution, cov,  out_dir=None):
	
	liste_size_sorted = sorted(list(dico_size_distribution.keys()))
	dico_which_files = {}
	
	for elem in liste_size_sorted:
		liste_match = []
		
		for sub_ in range(round_up(cov*elem), int(elem/cov), 1):
		
				if sub_ in liste_size_sorted:
					
					liste_match.append('{}/{}.faa'.format(out_dir, str(sub_)))
		
		dico_which_files['{}/{}.faa'.format(out_dir, str(elem))] = liste_match
	
	return dico_which_files


# I copied paste it
def blast(query, database, output_file, output_format_arg, other_args, blastp_):
	"""Just launch the blast"""
	
	cmd = ['{}'.format(blastp_), '-query', '{}'.format(query), '-db', '{}'.format(os.path.abspath(database)),
	       '-out', '{}'.format(output_file),
	       '-outfmt']
	cmd.append(output_format_arg + ' ')
	cmd.extend(other_args.split(' '))
	# if ban_list:
	# 	cmd.extend(['-negative_gilist', '{}'.format(ban_list)])
	# if db_size:
	# 	cmd.extend(['-dbsize', '{}'.format(db_size)])
	try:
		child = subprocess.Popen(cmd)
		child.wait()
	except:
		raise
	else:
		return 0
	
	
@time_reporteur_with_return
def split_dictionnary_for_each_thread(dictionary, number_of_thread):
	
	list_keys = list(dictionary.keys())
	indice = 0
	
	list_dict = []
	for i in range(number_of_thread):
		list_dict.append({})
	
	for key in list_keys:
		
		if indice >= number_of_thread:
			indice = 0
		
		actual_dico = list_dict[indice]
		actual_dico[key] = dictionary[key]
		del dictionary[key]
		
		indice += 1
	
	return list_dict

	

def process(args_tuple):
	directory, sub_dico, makeblast_db, blastp_, final_output, output_format_arg, other_args, cov = args_tuple
	# directory -> the directory where it will work.
	# sub_dico -> this argument is a sub dico made from the big dico.
	# it says which files it have to process against which other files.
	
	
	database_dir = directory + '/db'
	result_dir = directory + '/results'
	
	os.mkdir(directory)
	os.mkdir(database_dir)
	os.mkdir(result_dir)
	
	concatene_fasta_for_db = directory + '/concaten_fasta'
	tmp = directory + '/temporary'
	database = database_dir + "/db"
	big_blast = directory + '/{}'.format(os.path.basename(final_output))
	
	# Tha main loop
	
	for query_file in sub_dico:

		# regroup all the fasta for databas
		child = subprocess.Popen("cat {} > {}".format(' '.join(sub_dico[query_file]), concatene_fasta_for_db), shell=True)
		child.wait()
		
		# make the database
		child = subprocess.Popen("{} -dbtype prot -hash_index -in {} -out {}".format(makeblast_db,
		                                            concatene_fasta_for_db , database), shell=True)
		child.wait()

		
		# launch the blast
		blast(query=query_file, output_file=result_dir + '/{}'.format(os.path.basename(query_file)), database=database, blastp_=blastp_,
		      output_format_arg=output_format_arg, other_args=other_args)
	
		# may be cleaning here to limit size of the output
		# construct the awk command
		cmd = "awk "
		cmd += " \'($4>=20 && ($7-$6)/$10>={0} && ($9-$8)/$11>={0} && $1!=$2) ".format(cov)
		cmd += " {print $0}\' "
		
		cmd += "{0} > {1} && mv {1} {0}".format(result_dir + '/{}'.format(os.path.basename(query_file)), tmp)
		
		# launch the clean
		child = subprocess.Popen(cmd, shell=True)
		child.wait()
		
		# cat on the big blast
		child = subprocess.Popen("cat {} >> {}".format(result_dir + '/{}'.format(os.path.basename(query_file)), big_blast), shell=True)
		child.wait()
		
		# some clean
		child = subprocess.Popen("rm {} {}* {}".format(result_dir + '/{}'.format(os.path.basename(query_file)), database, concatene_fasta_for_db),
		                         shell=True)
		child.wait()
	
	return 0

@time_reporteur_with_return
def launch_mp(dico, nb_thread, working_directory, makeblast_db, blastp_, output_format_arg, other_args, final_output_, cov):
	
	sub_dico_list = split_dictionnary_for_each_thread(dictionary=dico, number_of_thread=nb_thread)
	
	# list with tuple for each process as argument # because of map
	liste_argument = []
	liste_final = []
	liste_wd = []
	#directory, sub_dico, makeblast_db, blastp_, final_output, output_format_arg, other_args
	
	indice = 1
	
	# we will launch as many thread than dictionary
	for sub_dico in sub_dico_list:
		directory_wd = working_directory +'/{}_working'.format(indice)
		
		liste_wd.append(directory_wd)
		
		liste_wd.append(directory_wd)
		
		final_output = working_directory + '/sub_blast{}.blastp'.format(indice)
		liste_final.append(working_directory + '/{}_working'.format(indice) + '/sub_blast{}.blastp'.format(indice))
		
		args_tuple = (directory_wd, sub_dico, makeblast_db, blastp_, final_output, output_format_arg, other_args, cov)
		
		liste_argument.append(args_tuple)
		indice +=1
	
	pool = Pool(processes=nb_thread)
	out_code = pool.map(process, liste_argument)
	
	child = subprocess.Popen("cat {} >> {}".format(' '.join(liste_final), final_output_), shell=True)
	child.wait()
	
	child = subprocess.Popen("rm {}".format(' '.join(liste_final)), shell=True)
	child.wait()
	
	return liste_wd