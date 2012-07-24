#!/usr/bin/env python

import subprocess
import os
import sqlite3  
from string import whitespace
from sys import stdout
from argparse import ArgumentParser
from operator import itemgetter
from collections import defaultdict, OrderedDict


class PairwiseAlignment:
	"""Structure to store the numbers associated with a Needleman Wunch pairwise alignment.
	Keeps track of two sequence ids (represented by Sequence object) and associated alignment score (float)"""
	def __init__(self, query, hit, score):	
		self.query = query
		self.hit = hit
		self.score = score

class Sequence: 
	"""Represents a FASTA sequence string with its header"""
	def __init__(self, header, data):
		self.header=header
		self.data=data
	
	"""Return an indetifier from the fasta sequence
	First non-whitespace string, excluding the first character (>)"""
	def name (self):
		return self.header.split()[0][1:]

	"""Return the complete FASTA sequence with both header and sequence data
	"""
	def fasta (self):
		return "\n".join([self.header,self.data])

class Identification:
	"""Represents a functional identification of a sequence toward an EC number
	Hypotheses is a possibly redundant list list of Hypothesis objects.
	The probability of a hypothesis being correct is calculated using the Bayes theorem.
	Please address Hung et al. (2010) for more details.
		Hung S, Wasmuth J, Sanford C & Parkinson J.
		DETECT - A Density Estimation Tool for Enzyme ClassificaTion and its application to Plasmodium falciparum.
		Bioinformatics. 2010 26:1690-1698
	This probability represents a singular alignment match event.
	Predictions is a non-redundant set of ec numbers associated with cumulative probabilities.
	The probability of a prediction is a cumulative probability of all hypotheses with the same EC number.
	"""
	def __init__(self, query_id):
		self.query_id = query_id
		self.hypotheses = list()
		self.predictions = defaultdict(self.__one)
		self.prediction_count = defaultdict(int)

	"""A callable function to initiate new values in a defaultdict to float 1
	"""
	def __one(self):
		return 1.0

class Hypothesis:
	"""Represents a single alignemnt result with an associated probability, as calculted using the Bayes theorem.
	EC is retrived from the swiss-to-EC mapping database.
	"""
	def __init__(self,swissprot_id,score):
		self.swissprot_id = swissprot_id
		self.score = score
		self.ec = "unknown"
		self.probability= 0.0
verbose=False
zero_density = 1e-10
"""Small number that is used as zero"""

def run_pair_alignment (seq, blast_db, num_threads, e_value_min, bitscore_cutoff):
	"""Core alignment routine.
	1) Takes a single sequence, acquires multiple BLASTp alignemnts to the swissprot enzyme database.
	2) Canonical sequences of the resutls from (1) are retrieved with blastdbcmd
	3) Original query is globally aligned versus sequences from (2)

	"""
	
	#First pass cutoff with BLAST alignments
	if verbose: print "[DETECT]: Running BLASTp for {} ...".format(seq.name())
	p = subprocess.Popen(("blastp", "-query", "-", 
					"-out", "-",
					"-db", blast_db,
					"-outfmt", "6 sseqid bitscore",
					"-max_target_seqs", "100000",
					"-num_threads",str(num_threads),
					"-evalue", str(e_value_min)),
				stdin=subprocess.PIPE,	
				stdout=subprocess.PIPE)
	stdout,stderr = p.communicate(seq.data)
	
	with open("blast_hits","w") as blast_hits:
		
		blast_hit_list = list()	
		for line in stdout.split("\n"):
			if not line in whitespace:
				swissprot_id,bitscore = line.split("\t")
				#sprot identifiers are sp|<ID>|<extra>
				seq_id = swissprot_id.split("|")[1]
				if float(bitscore) > bitscore_cutoff:
					blast_hit_list.append(seq_id)
		blast_hit_ids = "\n".join(blast_hit_list)
			
		if verbose: print "[DETECT]: Found {} hits for {} ...".format(len(blast_hit_ids),seq.name())
		
		#stop if nothing found
		if len(blast_hit_ids) == 0:
			return list()
		
		p = subprocess.Popen(("blastdbcmd", "-db", blast_db,
							"-entry_batch", "-"),
						stdout=subprocess.PIPE,
						stdin=subprocess.PIPE)
		
		stdout,stderr = p.communicate(blast_hit_ids)
		blast_hits.write(stdout)

	if verbose: print "[DETECT]: Running Needleman-Wunch alignments for {} ...".format(seq.name())

	#Run Needleman-Wunch alignment on the results of the BLAST search
	p = subprocess.Popen(("needle", "-filter",
					"-bsequence", "blast_hits",
					"-gapopen", "10",
					"-gapextend", "0.5",
					"-aformat_outfile", "score"),
				stdin=subprocess.PIPE,
				stdout=subprocess.PIPE)
		
	stdout,stderr = p.communicate(seq.fasta())
	os.remove("blast_hits")
	return parse_needle_results(stdout)
		

"""Split a fasta file into separate sequences, 
	return a list of Sequence objects.
	See class definition below"""

def split_fasta(input_file):
	#resultant array of peptide sequences
	sequences=list()

	#temporary array of lines for each sequence as it is being read
	seq_buffer=list()
	
	header = ""
	for line in open(input_file):
			
		#if we find a header line and there are already lines in sequence buffer
		if ">" in line and seq_buffer:
			
			#flush lines from the buffer to the sequences array
			sequences.append(Sequence(header,"".join(seq_buffer)))
			seq_buffer=list()
			header = line.strip()

		#skip the ID line
		elif ">" in line and not seq_buffer:
			header = line.strip()
		
		#add next line to sequence buffer
		else:
			seq_buffer.append(line.strip())

	#dont forget the last sequence
	sequences.append(Sequence(header,"".join(seq_buffer)))
	return sequences

"""Parse tab-delimeted BLAST results,
	return only hit IDs.
	Output generated with blastp -outfmt 6 (NCBI BLAST 2.2.26)
	BLAST docs http://www.ncbi.nlm.nih.gov/books/NBK1763 
		output format arguments: section 4.2.26"""
def get_blast_hit_identifiers (input_file):
	
	#resultant array of hit identifiers
	hit_IDs=list()
	
	#results are stored in-string as <query_id>\t<hit_id>\t...

	for line in open(input_file):
		hit_ID = line.strip().split("\t")[1]
		hit_IDs.append(hit_ID)
	return hit_IDs

"""Parse EMBOSS-needle output , 
	return list of structs.
	Output generated with needle -auto -aformat_outfile score (EMBOSS suite 6.4.0.0)
	Details of SCORE format: http://emboss.sf.net/docs/themes/AlignFormats.html
	Needle docs at http://emboss.sourceforge.net/apps/cvs/emboss/apps/needle.html"""
def parse_needle_results (needle_results):
	results=list()

	#results are stored in-string as <query> <hit> <alignment_length> (<score>)
	for line in needle_results.split("\n"):
		#ignore comment lines
		if not "#" in line and not line in whitespace:
			fields = line.strip().split()
			query = fields[0]
			hit = fields[1]
			#Score is printed in brackets - take them off, parse to float
			score = float(fields[3][1:-1])
			h = Hypothesis(hit,score)
			results.append(h)
	return results

def calculate_probability (hypothesis,db_connection):
	

	score = hypothesis.score
		
	cursor = db_connection.cursor()

	#Fetch the data from tables
	
	#Fetch the EC mapped to the swissprot ID
	cursor.execute("SELECT ec FROM swissprot_ec_map WHERE swissprot_id = '{}'".format(hypothesis.swissprot_id))
	#sqlite3.fetchone() returns a tupule. Since only one value <ec> was requested, this is a one-member tuple. Still, it is important to subset [0]
	mapping = cursor.fetchone()
	if mapping:
		ec = mapping[0]
		hypothesis.ec = ec

		#Get Prior probabilities for that EC
		cursor.execute("SELECT probability FROM prior_probabilities WHERE ec = '{}'".format(ec))
		prior = cursor.fetchone()[0]

		#Get positive density for the given score and EC
		cursor.execute("SELECT density FROM positive_densities WHERE ec = '{}' AND score < {} LIMIT 1".format(ec,score))
		previous_point = cursor.fetchone()

		cursor.execute("SELECT density FROM positive_densities WHERE ec = '{}' AND score > {} LIMIT 1".format(ec,score))
		next_point = cursor.fetchone()
		
		if previous_point and next_point:
			positive = (previous_point[0] + next_point[0])/2
		else:
			positive = 0
			

		#Get negative density for the given score and EC
		cursor.execute("SELECT density FROM negative_densities WHERE ec = '{}' AND score < {} LIMIT 1".format(ec,score))
		previous_point = cursor.fetchone()
	
		cursor.execute("SELECT density FROM negative_densities WHERE ec = '{}' AND score > {} LIMIT 1".format(ec,score))
		next_point = cursor.fetchone()
		
		if previous_point and next_point:
			negative = (previous_point[0] + next_point[0])/2
		else:
			negative = 0

		if positive == 0 and negative == 0:
			probability = zero_density
		else:
			positive_hit = prior * positive
			probability = positive_hit / (positive_hit + ((1.0-prior) * negative ))

		hypothesis.probability = probability
	else:
		probability=0


	return probability
	

if __name__=="__main__":
	parser = ArgumentParser(description="DETECT - Density Estimation Tool for Enzyme ClassificaTion. Version 2.0. June 2012")
	
	parser.add_argument("target_file",type=str,help="Path to the file containing the target FASTA sequence(s)")
	parser.add_argument("--output_file",type=str,help="Path of the file to contain the output of the predictions")
	parser.add_argument("--tabular_output",help="Print output in tab-separated form",action="store_true")
	parser.add_argument("--verbose",help="Print verbose output",action="store_true")
	parser.add_argument("--sort_output",help="Sort the prediction probabiliites",action="store_true")
	parser.add_argument("--num_threads",type=int,help="Number of threads used by BLASTp")
	parser.add_argument("--bit_score",type=float,help="The cutoff for BLASTp alignment bitscore")
	parser.add_argument("--e_value",type=float,help="The cutoff for BLASTp alignment E-value")
	
	args = parser.parse_args()
	script_path = os.path.dirname(os.path.realpath(__file__))

	verbose = args.verbose
	num_threads = args.num_threads if args.num_threads else 1
	bit_score = args.bit_score if args.bit_score else 50
	e_value = args.e_value if args.e_value else 1
	
	sequences = split_fasta(args.target_file)
	if verbose: print "Found {} sequences in file.".format(len(sequences))
	blast_db = script_path+"/data/uniprot_sprot.fsa"
	
	final_predictions = list()

	connection = sqlite3.connect(script_path+"/data/detect.db")
	for i,seq in enumerate(sequences):
		
		if verbose: 
			print "[DETECT]: Analyzing {} ({}/{}) ...".format(seq.name(),i+1,len(sequences))
		identification = Identification(seq.name())
		identification.hypotheses = run_pair_alignment(seq,blast_db,num_threads,e_value,bit_score)
		
		if not identification.hypotheses:
			if verbose: 
				print "[DETECT]: No BLASTp hits for {}".format(seq.name())
			continue

		if verbose: 
			print "[DETECT]: Running density estimations for {} ...".format(seq.name())
		for hypothesis in identification.hypotheses:

			probability = calculate_probability(hypothesis,connection)
			if not (hypothesis.ec == "unknown" or hypothesis.ec == "2.7.11.1" or hypothesis.ec == "2.7.7.6" or hypothesis.ec == "2.7.13.3"):
				identification.predictions[hypothesis.ec] *= (1.0-probability)	
				identification.prediction_count[hypothesis.ec] += 1

		for ec,probability in identification.predictions.items():
			cumulative = 1.0 - probability
			if (cumulative > zero_density):
				identification.predictions[ec] = cumulative
			else:
				del identification.predictions[ec]
		
		#sort dict by values
		if (args.sort_output):
			identification.predictions = OrderedDict(sorted(identification.predictions.iteritems(),key=itemgetter(1),reverse=True))
		
		final_predictions.append(identification)
		if verbose: 
			print "[DETECT]: Identified {} predictions for {}".format(len(identification.predictions.keys()),seq.name())
		
	if (args.output_file):
		output = open(args.output_file,"w")
	else:
		output = stdout
	if args.tabular_output:
		output.write("ID\tEC\tprobability\tpositive_hits\tnegative_hits\n")
		for identification in final_predictions:
			for ec in identification.predictions:
				output.write("{seq_id}\t{pred_ec}\t{prob:.3e}\t{pos_hits}\t{neg_hits}\n".format(
							seq_id=identification.query_id,
							pred_ec=ec,
							prob=identification.predictions[ec],
							pos_hits=identification.prediction_count[ec],
							neg_hits= len(identification.hypotheses)-identification.prediction_count[ec]))
	else:
		for identification in final_predictions:
			output.write("{:^48}\n".format("DETECT report for {}".format(identification.query_id)))
			output.write("Predicted EC number:  |  Cumulative probability:  |  Positive Hits:\n")
			for ec in identification.predictions:
				output.write("{pred_ec:^20}  |  {prob:^23.3e}  |  {pos_hits:^16}\n".format(
							pred_ec=ec,
							prob=identification.predictions[ec],
							pos_hits=identification.prediction_count[ec]))

	output.close()
	connection.close()
