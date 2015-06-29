#!/usr/bin/python

# See http://bib.oxfordjournals.org/content/13/5/569.full
# For definitions

import sys
import time
import networkx as nx
import math

# Helper functions
def announce(message):
	"""  Timestamped output to stdout """
	print time.strftime('%H:%M'),message
	sys.stdout.flush()

def open_or_abort(filename, option='r'):
	""" Output error message to stderr if file opening failed """
	try:
		newfile = open(filename, option)
	except IOError:
		sys.stderr.write("Could not open {} -- Aborting\n".format(filename))
		exit()
	return newfile


# NOTE(tfs): Accepted GO file format:
#
#	! comments
#
#	[Term]
#	id: GO_term
#	...
#	is_a: GO_term
#	is_a: GO_term
#
#	[Term]
#	...
#
# 	[Typedef]
#
# The [Typedef] tag signals end of GO terms.
# It is necessary in the current implementation
#
def parse_go_file(go_file_name):
	""" Parses and returns (does not natively store) GO data """
	
	go_file = open_or_abort(go_file_name)

	# Setup
	go_file.seek(0)
	go_graph = nx.DiGraph()

	alt_ids = {}

	go_term = ''
	parents = []
	is_obsolete = False

	# Don't start paying attention until we see '[Term]'
	valid_to_read = False

	# Main parsing loop
	for line in go_file:

		# Only if we're within a '[Term]' header
		if valid_to_read:
			if line.startswith('alt_id:'):
				alt_id = line.strip()[8:]
				alt_ids[alt_id] = go_term

			elif line.startswith('id:'):

				# Only log if the entry is valid
				if not is_obsolete:
					if go_term != '':

						# Only add connected node
						if len(parents) > 0:
							go_graph.add_node(go_term)

							for parent in parents:
								if parent != '':
									go_graph.add_edge(parent, go_term)

				# Reset regardless of logging status
				parents = []
				is_obsolete = False

				go_term = line.strip()[4:]
			elif line.startswith('is_a:'):

				# Store is_a as a parent
				parents.append(line.split('!')[0].strip()[6:])

			elif line.startswith('is_obsolete: true'):

				# Do not store the data under this '[Term]' header
				is_obsolete = True

			elif line.startswith('[Typedef]'):
				# Write if the previous entries were valid

				# Only log if the entry is valid
				if not is_obsolete:
					if go_term != '':

						# Only add connected node
						if len(parents) > 0:

							go_graph.add_node(go_term)

							for parent in parents:
								if parent != '':
									go_graph.add_edge(parent, go_term)

				# Reset regardless of logging status
				parents = []
				is_obsolete = False

				go_term = '' # No valid ID to reset with under a '[Typedef]' header

				# Stop paying attention
				valid_to_read = False
		else:
			if '[Term]' in line:

				# Start paying attention
				valid_to_read = True

	go_file.close()

	return (go_graph, alt_ids)


# NOTE(tfs): Accepted AC file format:
#
#	-
#	protein_name
#	GO_term
#	GO_term
#	GO_term
# 	-
#
def parse_annotation_corpus(ac_file_name):
	""" Parses annotation corpus. Returns a dictionary of { gene: [terms] } """

	ac_file = open_or_abort(ac_file_name)

	# Setup
	prot_to_gos = {}
	go_to_prots = {}

	ac_file.seek(0)

	curr_prot = ''
	curr_gos = []
	new_entry = True

	for line in ac_file:
		
		# Start information from new entry
		if line.startswith('-'):

			# Only update if we have enough information for the last entry
			if curr_prot != '' and len(curr_gos) > 0:

				# Update prot_to_gos
				if curr_prot in prot_to_gos:
					prot_to_gos[curr_prot] = prot_to_gos[curr_prot] + curr_gos
				else:
					prot_to_gos[curr_prot] = curr_gos

				# Update go_to_prots
				for go in curr_gos:
					if go in go_to_prots:
						go_to_prots[go].append(curr_prot)
					else:
						go_to_prots[go] = [curr_prot]

			# Reset, regardless of whether or not we updated
			curr_prot = ''
			curr_gos = []
			new_entry = True

		# If we've just started looking at a new entry, parse as protein name
		# DON'T do this if we're still on the delimiter line ('-')
		elif new_entry:
			curr_prot = line.strip().strip(';')
			new_entry = False

		# Otherwise, parse as GO term
		else:
			curr_gos.append(line.strip().strip(';'))

	ac_file.close()

	return (prot_to_gos, go_to_prots)





###############################
### SemSim_Calculator class ###
###############################


class SemSimCalculator():
	"""
		Stores GO and annotation corpus data internally.
		Calculates different semantic similarity metrics.
	"""

	def __init__(self, go_file_name, ac_file_name):
		""" Initialize using GO and annotation corpus files (pass in file name, not file object) """
		self._go_graph, self._alt_list = parse_go_file(go_file_name)
		self._prot_to_gos, self._go_to_prots = parse_annotation_corpus(ac_file_name)
		self._proteins = [x[0] for x in self._prot_to_gos.items()]
		self._num_proteins = len(self._proteins)
		self._ic_vals = {} # For memoizing IC values (they are unchanging given an ontology and annotation corpus)
		
		self._go_terms = self._go_graph.nodes()

	def get_go_graph(self):
		""" Return nx graph for GO """

		return nx.DiGraph(self._go_graph)

	def get_alt_list(self):
		""" Return alt_list """

		return dict(self._alt_list)

	def get_ptg(self):
		""" Return copy of prot_to_gos """

		return dict(self._prot_to_gos)

	def get_gtp(self):
		""" Return copy of go_to_prots """

		return dict(self._go_to_prots)

	def get_proteins(self):
		""" Return copy of proteins """

		return list(self._proteins)

	def get_num_proteins(self):
		""" Return number of proteins """

		return int(self._num_proteins)

	def get_ic_vals(self):
		"""
			Return all stored ic_vals.
			Not all values are guaranteed to exist.
			Consider running precompute_ic_vals first.
		"""

		return self._ic_vals

	def calc_term_prob(self, term):
		""" Probability of term or desc(term) to occur as a label within the annotation corpus """

		if term == None or (not term in self._go_graph):
			return None

		# Find all descendants of term, including term
		terms = nx.algorithms.dag.descendants(self._go_graph, term)
		terms.add(term)

		annotated_proteins = {}

		# Mark any protein labeled with term or a descendant of term
		for term in terms:
			if term in self._go_to_prots:
				for prot in self._go_to_prots[term]:
					annotated_proteins[prot] = True

		prob = float(len(annotated_proteins.items())) / float(self._num_proteins)

		return prob

	def IC(self, term):
		""" Information content: IC(c) = -log(p(c)) """

		# Check if IC has been computed for term already
		if not (term in self._ic_vals):
			# If not seen before, compute IC
			prob = self.calc_term_prob(term)

			if prob == 0 or prob == None:
				self._ic_vals[term] = None
				return None
			else:
				ic = (-1) * math.log(prob)
				self._ic_vals[term] = ic # Memoize IC value
				return ic
		else:
			# If seen before, return memoized value
			return self._ic_vals[term]

	def precompute_ic_vals(self):
		""" Compute and store IC values for all ontology terms """

		for term in self._go_graph.nodes():
			self.IC(term)

	def MICA(self, left, right):
		"""
			Maximum Informative Common Ancestor:
			MICA(t1, t2) = arg max, IC(tj)
							tj in ancestors(t1, t2)

			(returns a term, common ancestor of left and right)
		"""

		if not left in self._go_terms:
			if left in self._alt_list:
				left = self._alt_list[left]
			else:
				return None

		if not right in self._go_terms:
			if right in self._alt_list:
				right = self._alt_list[right]
			else:
				return None

		# Find common ancestors as intersection of two ancestor sets
		# NOTE(tfs): Python sets are very slow. List comprehensions are faster
		left_ancs = nx.algorithms.dag.ancestors(self._go_graph, left)
		right_ancs = nx.algorithms.dag.ancestors(self._go_graph, right)
		ancestors = [a for a in left_ancs if a in right_ancs]

		# Edge case where left and right are the same. Treat left and right as a common ancestor
		if left == right:
			ancestors.append(left)

		max_term = None
		max_IC = 0

		# Calculate IC for all ancestors; store maximum IC value and term
		for ancestor in ancestors:
			anc_IC = self.IC(ancestor)
			if anc_IC != None and anc_IC > max_IC:
				max_IC = anc_IC
				max_term = ancestor

		return max_term

	def simRes(self, left, right):
		"""
			simRes(t1, t2) = IC[MICA(t1, t2)]
			Returns a value (IC result)
		"""

		return self.IC(self.MICA(left, right))

	def simLin(self, left, right):
		"""
			simLin(t1, t2) = [IC[MICA(t1, t2)]] / [IC(t1) + IC(t2)]
			Returns a value
			Currently untested
		"""

		return self.IC(self.MICA(left, right)) / (IC(go_graph, go_to_prots, left) + IC(go_graph, go_to_prots, right))


	def simJC(self, left, right):
		"""
			simJC(t1, t2) = 1 - IC(t1) + IC(t2) - 2xIC[MICA(t1, t2)]
			Returns a value
			Currently untested
		"""

		return 1 - self.IC(left) + self.IC(right) - (2*self.IC(self.MICA(left, right)))

	def pairwise_average_term_comp(self, lefts, rights, metric):
		"""
			Compares each pair of terms in two sets or lists of terms.
			Returns the average of these comparison scores.
			Uses metric(left, right) to make each comparison.
			metric must take in two ontology terms (left and right) and return a numeric score.
		"""

		total_score = 0
		num_scores = 0
		for left in lefts:
			for right in rights:
				new_score = metric(left, right)

				# Count a new_score of None in the denominator, but treat it as a value of 0
				# This mimics a dummy root node if there are multiple roots in the GO DiGraph
				if new_score != None:
					total_score += new_score
				num_scores += 1

		if total_score == 0:
			return None
		else:
			return total_score / num_scores

	def pairwise_max_term_comp(self, lefts, rights, metric):
		"""
			Compares each pair of terms in two sets or lists of terms.
			Returns the maximum score found in these comparisons.
			Uses metric(left, right) to make each comparison.
			metric must take in two ontology terms (left and right) and return a numeric score.
		"""

		max_score = 0

		for left in lefts:
			for right in rights:
				temp_score = metric(left, right)
				if temp_score != None and temp_score > max_score:
					max_score = temp_score

		return max_score

	def average_protein_comp(self, left_prot, right_prot, metric):
		"""
			Looks up all go terms for left_prot and right_prot.
			Uses pairwise_average_term_comp to compare the above sets of terms.
			metric must take in two ontology terms (left and right) and return a numeric score.
		"""

		left_terms = self._prot_to_gos[left_prot]
		right_terms = self._prot_to_gos[right_prot]

		return self.pairwise_average_term_comp(left_terms, right_terms, metric)

	def max_protein_comp(self, left_prot, right_prot, metric):
		"""
			Looks up all terms for left_prot and right_prot.
			Uses pairwise_max_term_comp to compare the above sets of terms.
			metric must take in two ontology terms (left and right) and return a numeric score.
		""" 

		left_terms = self._prot_to_gos[left_prot]
		right_terms = self._prot_to_gos[right_prot]

		return self.pairwise_max_term_comp(left_terms, right_terms, metric)



####################################
###  End SemSim_Calculator class ###
####################################
