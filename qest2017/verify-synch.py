#!/usr/bin/env python3
################################################################################
# synchronisation of pulse-coupled oscillators
#
#
# author	: Paul Gainer
# email		: p.gainer@liverpool.ac.uk
################################################################################

################################################################################
# pre-python changelog
#
# 20/05/16 - initial test model
# 02/06/16 - added improved synchronisation (extra alphas)
# 01/12/16 - added refractory period
# 11/01/17 - removed epsilons as parameters
# 19/01/17 - added message loss
# 26/01/17 - added order parameter calculations
################################################################################
# python changelog
#
# 27/01/17 - migrated model generation to python
# 03/02/17 - message loss now applies to all messages
# 15/02/17 - removed order parameter constants/formulae & some reward structures
################################################################################
################################################################################
# python bugfixes
#
# 05/02/17 - fixed rounding error when generating parameter sequences
################################################################################

import errno
import itertools
import math
import numpy
import operator
import os
import re
import subprocess
import sys
import time


################################################################################
# configuration for:
#
# Firefly-Inspired Synchronization in Swarms of Mobile Agents
# Perez-Diaz, Zillmer, Gross
################################################################################
MODEL_NAME = 'swarm-perez-diaz'
#===============================================================================
# value = [value_min, value_max, value_step]
#
# Other parameter values can be used to define the values for any parameter,
# for instance, RP = [0, 't - 2', 1] will generate models with # refractory
# periods from 0 up to T - 2 etc. When referring to other parameter values put
# the expression in a string and use lowercase letters, so n, t, epsilon etc....
#===============================================================================
N = [6, 6, 1]
T = [11, 11, 1]
EPSILON = [0.1, 1.0, 0.1]
RP = [1, 1, 1]
ML = [0.1, 0.1, 0.1]
#===============================================================================
# evol(x)
#===============================================================================
def evolution(x):
	return x + 1
#===============================================================================
# pert(x)
#===============================================================================
def perturbation(t, x, epsilon, alpha):
	return x * epsilon * alpha
#===============================================================================
# properties to check in the models
#===============================================================================
def create_properties(n, t, rp, epsilon, ml):
	properties = ['P=? [F synchronised]', 'R{"time_to_synch"}=? [F synchronised]']
	return properties
################################################################################

'''
################################################################################
# configuration for:
#
# Bio-Inspired Tools for a Distributed Wireless Sensor Network Operating System
# Breza
################################################################################
MODEL_NAME = 'figo-breza'
#===============================================================================
# value = [value_min, value_max, value_step]
#
# Other parameter values can be used to define the values for any parameter,
# for instance, RP = [0, 't - 2', 1] will generate models with # refractory
# periods from 0 up to T - 2 etc. When referring to other parameter values put
# the expression in a string and use lowercase letters, so n, t, epsilon etc....
#===============================================================================
N = [5, 5, 1]
T = [10, 10, 2]
EPSILON = [1.0, 1.0, 0.05]
RP = [0, 't', 1]
ML = [0.0, 0.0, 0.02]
#===============================================================================
# evol(x)
#===============================================================================
def evolution(x):
	return x + 1
#===============================================================================
# pert(x)
#===============================================================================
def perturbation(t, x, epsilon, alpha):
	return ((x + (((2 ** alpha) - 1) * t)) / (2 ** alpha)) - x
#===============================================================================
# properties to check in the models
#===============================================================================
def create_properties(n, t, rp, epsilon, ml):
	properties = ['P=? [F synchronised]', 'R{"time_to_synch"}=? [F synchronised]']
	return properties
################################################################################	
'''

#-------------------------------------------------------------------------------
# prism configuration
#-------------------------------------------------------------------------------
PRISM = 'prism64'
PRISM_ENGINE = 'hybrid'
PRISM_PARAMS = ['-cuddmaxmem', '4g', '-javamaxmem', '12g', '-' + PRISM_ENGINE, '-maxiters', '100000']
TIME_COMMAND = '/usr/bin/time'
TIME_PARAMS = '--verbose'

#-------------------------------------------------------------------------------
# set to True to write PRISM output to terminal
#-------------------------------------------------------------------------------
VERBOSE = False

#-------------------------------------------------------------------------------
# set to True to ignore non-flashing states
#-------------------------------------------------------------------------------
SKIP_NON_FLASHING_STATES = True

#-------------------------------------------------------------------------------
# remove generated PRISM models
#-------------------------------------------------------------------------------
REMOVE_GENERATED_MODELS = True

#-------------------------------------------------------------------------------
# directories, sys date and time, results files
#-------------------------------------------------------------------------------
MODELS_DIR = './models/%s/' % MODEL_NAME
RESULTS_DIR ='./results/%s/' % MODEL_NAME
DATE = time.strftime('%d-%m-%Y')
TIME = time.strftime('%H-%M-%S')
RESULTS_PREFIX = '%sresults-%s_%s_%s' % (RESULTS_DIR, MODEL_NAME, DATE, TIME)
FULL_RESULTS_FILE_NAME = '%s.full' % RESULTS_PREFIX
CSV_RESULTS_FILE_NAME = '%s.csv' % RESULTS_PREFIX

#-------------------------------------------------------------------------------
# terms to search for in prism/time output
#-------------------------------------------------------------------------------
FULL_RESULTS_TERMS = ['Result', 'Reachability', 'Time for model construction', 'States', 'Transitions', 'Transition matrix', 'Time for model checking', 'Maximum resident set size']
CSV_RESULTS_TERMS = ['Result', 'Time for model construction', 'States', 'Transitions', 'Transition matrix', 'Time for model checking', 'Maximum resident set size']

#-------------------------------------------------------------------------------
# global variables
#-------------------------------------------------------------------------------
full_results_file = None
csv_results_file = None
prism_model_file = None
binomials = None
min_alphas = None
max_alphas = None
write_file = None

#------------------------------------------------------------------------------
# 'get_sequence' statements to calculate intervals of parameter values
#-------------------------------------------------------------------------------
N_VALUES = 'get_sequence(%s, int)' % ','.join(map(str, N)) 
T_VALUES = 'get_sequence(%s, int)' % ','.join(map(str, T)) 
EPSILON_VALUES = 'get_sequence(%s, float)' % ','.join(map(str, EPSILON)) 
RP_VALUES = 'get_sequence(%s, int)' % ','.join(map(str, RP)) 
ML_VALUES = 'get_sequence(%s, float)' % ','.join(map(str, ML)) 


################################################################################
#===============================================================================
# builds models for given parameter ranges and checks properties in the models
#===============================================================================
################################################################################
def main(argv):
	global full_results_file, csv_results_file, results_prefix
	# create subdirectories if non exist, and create the results file
	create_dir(MODELS_DIR)
	create_dir(RESULTS_DIR)	
	full_results_file = open(FULL_RESULTS_FILE_NAME, 'w')
	csv_results_file = open(CSV_RESULTS_FILE_NAME, 'w')
	# build a new model for each set of parameter values
	for n in eval(N_VALUES):
		for t in eval(T_VALUES):
			calculate_binomials(n, t);
			for rp in eval(RP_VALUES):
				for epsilon in eval(EPSILON_VALUES):
					build_alphas(n, t, epsilon)
					# generate the prism model...
					PRISM_FILE_NAME = '%smodel-%s-%d-%d-%g-%d.prism' % (MODELS_DIR, MODEL_NAME, n, t, epsilon, rp)						
					generate_prism_model_file(PRISM_FILE_NAME, n, t, epsilon, rp)
					for ml in eval(ML_VALUES):
						# ...and check the properties
						check_properties_in_prism(PRISM_FILE_NAME, n, t, epsilon, rp, ml)
					if REMOVE_GENERATED_MODELS:
						os.remove(PRISM_FILE_NAME)
	# close the results files
	full_results_file.close()
	csv_results_file.close()
	
#===============================================================================
# check all properties in the given prism model
#===============================================================================
def check_properties_in_prism(prism_file_name, n, t, epsilon, rp, ml):
	set_write_file(full_results_file)
	writeln('=====================================================================')
	writeln('N=%d T=%d epsilon=%g RP=%d ML=%g' % (n, t, epsilon, rp, ml))
	set_write_file(csv_results_file)
	write('%d,%d,%g,%d,%g' % (n, t, epsilon, rp, ml))
	print('for parameters: N=%d, T=%d, epsilon=%g, RP=%d, ML=%g' % (n, t, epsilon, rp, ml))
	for property in create_properties(n, t, rp, epsilon, ml):		
		print('\tverifying property: %s' % property)
		# time the prism process, we want some memory consumption info from time
		results = run_process(TIME_COMMAND, [TIME_PARAMS, PRISM, prism_file_name, '-pf', property, *PRISM_PARAMS, '-const', 'ML=%g' % ml]).split('\n')
		set_write_file(full_results_file)
		writeln('_____________________________________________________________________')
		writeln('results for checking property:')
		writeln(property)
		writeln('_____________________________________________________________________')							
		# extract all the relevant infprmation from the output
		for term in FULL_RESULTS_TERMS:
			for line in results:
				if line.startswith(term):
					writeln(line)
					break
		set_write_file(csv_results_file)
		# for the csv file we just want the numerical output
		for term in CSV_RESULTS_TERMS:
			found = False
			for line in results:
				line = line.strip()
				if line.startswith(term):
					write(',' + re.search('[0-9]+[0-9.E-]*|Infinity', line).group())
					found = True
					break
			if not found:
				write(',')
	set_write_file(csv_results_file)
	writeln('')										

#===============================================================================
# generates a prism model given a file name and parameters
#===============================================================================
def generate_prism_model_file(prism_file_name, n, t, epsilon, rp):
	if os.path.exists(prism_file_name):
		print('using preexisting model file: %s' % prism_file_name)
		return
	prism_model_file = open(prism_file_name, 'w')
	set_write_file(prism_model_file)
	print('building model file: %s' % prism_file_name);
	# add some info about built date/time and parameters used
	build_comment_lines(
		[
			'file generated on %s at %s' % (DATE, TIME),
			'',
			'model name: %s' % MODEL_NAME
		], "=")
	writeln('dtmc')
	empty_lines(1)
	build_comment_lines(
		[
			'parameters',
			'',
			'N       = %d' % n,
			'T       = %d' % t,
			'EPSILON = %g' % epsilon,
			'RP      = %d' % rp,										
		], "-")
	empty_lines(1)
	build_comment_lines(['rewards'], '=')
	build_rewards(n, t, epsilon, rp)
#	build_comment_lines(['order parameter'], '=')
#	build_order_parameter(n, t, epsilon, rp)
	build_comment_lines(['constants'], '=')
	build_constants(n, t, epsilon, rp)
	build_comment_lines(['formulae'], '=')
	build_formulae(n, t, epsilon, rp)
	build_comment_lines(['modules'], '=')
	build_initial_state_indicator_module(n, t, epsilon, rp)
	build_oscillator_population_module(n, t, epsilon, rp)
	prism_model_file.close()

#===============================================================================
# build the reward structures
#===============================================================================	
def build_rewards(n, t, epsilon, rp):
	build_comment_lines(['records the duration'], '-')
	#---------------------------------------------------------------------------
	# time to synch reward
	#---------------------------------------------------------------------------
	writeln('rewards "time_to_synch"')
	if not SKIP_NON_FLASHING_STATES:
		writeln('assigned & (!synchronised) : 1 / %d;' % t, tabs = 1)
	else:
		writeln('assigned & (n_%d > 0) & (!synchronised) : 1 / %d;' % (t, t), tabs = 1)
		for i in range(1, t):
			write('assigned & (n_%d = 0 & n_%d > 0' % (t, i), tabs = 1)
			for j in range(i + 1, t):
				write(' & n_%d = 0' % j)
			writeln(') & (!synchronised) : %d / %d;' % (t - i, t))
	writeln('endrewards')
	empty_lines(1)
	#---------------------------------------------------------------------------
	# cycles reward	
	#---------------------------------------------------------------------------
#	writeln('rewards "cycles"')
#	if not SKIP_NON_FLASHING_STATES:
#		writeln('assigned : 1 / %d;' % t, tabs = 1)
#	else:
#		writeln('assigned & (n_%d > 0) : 1 / %d;' % (t, t), tabs = 1)
#		for i in range(1, t):
#			write('assigned & (n_%d = 0 & n_%d > 0' % (t, i), tabs = 1)
#			for j in range(i + 1, t):
#				write(' & n_%d = 0' % j)
#			writeln(') : %d / %d;' % (t - i, t))
#	writeln('endrewards')
#	empty_lines(1)
	#---------------------------------------------------------------------------
	# mean flashes reward	
	#---------------------------------------------------------------------------
#	writeln('rewards "mean_flashes"')
#	writeln('is_non_initial_state : n_1 / %d;' % n, tabs = 1)
#	writeln('endrewards')
#	empty_lines(1)
	
	
#===============================================================================
# build the constants and forumlae for the order parameter
#===============================================================================
def build_order_parameter(n, t, epsilon, rp):
	# calculate x and y components of unit vectors
	for i in range(0, t):
		write('const double unit_vector_x_%d = %.8f; ' % (i + 1, math.cos(2 * math.pi * float((1 / t) * i))));
		writeln('const double unit_vector_y_%d = %.8f;' % (i + 1, math.sin(2 * math.pi * float((1 / t) * i))));
	empty_lines(1)
	write('formula unit_vector_x_avg_squared = pow(((');
	for i in range(1, t + 1):
		write('(unit_vector_x_%d * n_%d)%s' % (i, i, '' if i == t else ' + ' ));
	writeln(') / %d), 2);' % n);
	write('formula unit_vector_y_avg_squared = pow(((');
	for i in range(1, t + 1):
		write('(unit_vector_y_%d * n_%d)%s' % (i, i, '' if i == t else ' + '));
	writeln(') / %d), 2);' % n);
	writeln('formula order_parameter = pow(unit_vector_x_avg_squared + unit_vector_y_avg_squared, 0.5);');
	empty_lines(1)
	
#===============================================================================
# build the formulae
#===============================================================================	
def build_formulae(n, t, epsilon, rp):
	#---------------------------------------------------------------------------
	# num oscillators formula
	#---------------------------------------------------------------------------
	build_comment_lines(['the total number of oscillators at any given time'], '-')
	write('formula num_oscillators = ')
	for i in range(1, t):
		write('n_%d + ' % i)
	writeln('n_%d;' % t)
	empty_lines(1)
	#---------------------------------------------------------------------------
	# assigned formula
	#---------------------------------------------------------------------------
	build_comment_lines(['true after all oscillators have been assigned phases'], '-')
	write('formula assigned = num_oscillators = %d;' % n)
	empty_lines(1)
	#---------------------------------------------------------------------------
	# synchronised formula
	#---------------------------------------------------------------------------
	build_comment_lines(['true if all oscillators have the same phase'], '-')
	write('formula synchronised = ')
	for i in range(1, t):
		write('n_%d = %d | ' % (i, n))
	write('n_%d = %d;' % (t, n))
	empty_lines(1)
	#---------------------------------------------------------------------------
	# message loss probability formulae
	#---------------------------------------------------------------------------
	build_comment_lines(['message loss probabilities'], '-')
	for i in range(pow(2, n)):
		sigma = '_'.join(tuple('{:0{}b}'.format(i, n)))
		num_lost = sigma.count('1')
		writeln('formula pr_%s = (pow(ML, %d) * pow((1 - ML), %d));' % (sigma, num_lost, n - num_lost))
	empty_lines(1)
					
#===============================================================================
# build the constants
#===============================================================================
def build_constants(n, t, epsilon, rp):
	#---------------------------------------------------------------------------
	# infinity constant
	#---------------------------------------------------------------------------
	build_comment_lines(['infinity'], '-')
	writeln('const int INFINITY = %d;' % (n + 1))
	empty_lines(1)
	#---------------------------------------------------------------------------
	# message loss constant
	#---------------------------------------------------------------------------
	build_comment_lines(['message loss chance'], '-')
	writeln('const double ML;')
	empty_lines(1)

#===============================================================================
# build the min_alphas and max_alphas
#===============================================================================
def build_alphas(n, t, epsilon):
	global min_alphas, max_alphas
	min_alphas = numpy.empty([t + 1, t + 1], dtype = int)
	max_alphas = numpy.empty([t + 1, t + 1], dtype = int)
	l = []
	s_prime = -1
	for i in range(1, t):
		alpha = 0
		minimum = alpha
		s = evolution(i) + int(round_half_up(perturbation(t, i, epsilon, alpha)))
		while(s <= t and alpha < n):
			alpha = alpha + 1
			s_prime = evolution(i) + int(round_half_up(perturbation(t, i, epsilon, alpha)))
			maximum = alpha
			if s != s_prime:
				min_alphas[i][s] = minimum
				max_alphas[i][s] = maximum - 1								
				l = l + [(i, s)]
				minimum = maximum
				s = s_prime
		if s > t:
			min_alphas[i][1] = minimum;
			max_alphas[i][1] = n;
			l = l + [(i, 1)]
		else:
			min_alphas[i][s] = minimum;
			max_alphas[i][s] = maximum;		
			l = l+ [(i, s)]
	for i in range (1, t + 1):
		for j in range(1, t + 1):
			if not ((i, j) in l) and (j > i or j == 1):
				if i == t:
					min_alphas[i][j] = 0;
				else:
					min_alphas[i][j] = n + 1;	
				max_alphas[i][j] = n + 1;
					
#===============================================================================
# build the oscillator population module
#===============================================================================
def build_oscillator_population_module(n, t, epsilon, rp):
	build_comment_lines(['oscillator population module'], '-')
	writeln('module oscillator_population')
	#---------------------------------------------------------------------------
	# counters
	#---------------------------------------------------------------------------
	build_comment_lines(['counters for oscillators in each state'], '-', tabs = 1)	
	for i in range(1, t + 1):
		writeln('n_%d : [0..%d] init 0;' % (i, n), tabs = 1)
	empty_lines(1)
	#---------------------------------------------------------------------------
	# transitions from starting state
	#---------------------------------------------------------------------------
	build_comment_lines(['transitions from starting state'], '-', tabs = 1)
	build_initial_transitions(n, t, epsilon, rp)
	build_comment_lines(['update oscillator phases, no flashes'], '-', tabs = 1)
	build_update_oscillator_phases_no_flashes(n, t, epsilon, rp)
	build_comment_lines(['update oscillator phases, flashes'], '-', tabs = 1)
	build_update_oscillator_phases_flashes(n, t, epsilon, rp)
	writeln('endmodule')

#===============================================================================
# build the transitions for when alpha_T > 0
#===============================================================================
def build_update_oscillator_phases_flashes(n, t, epsilon, rp):
	pow_2_n = pow(2, n)
	# create the strings for the probability loss formulae
	# and the message loss masks
	masks = numpy.empty((pow_2_n,), dtype = list)
	pr_strings = numpy.empty((pow_2_n,), dtype = object)
	for r in range(pow_2_n):
		bin_r = '{:0{}b}'.format(r, n)
		pr_strings[r] = 'pr_%s' % ('_'.join(bin_r))
		masks[r] = list(map(int, '{:0{}b}'.format(r, n)))
	for sigma in combinations_with_replacement_counts(t, n):
		# for every configuration of oscillators with at least one oscillator
		# in the firing state
		if(sigma[t - 1] != 0):
			write('[step] (assigned', tabs = 1)
			for i in range(t):
				write(' & n_%d = %d' % (i + 1, sigma[i]))
			writeln(') ->')
			next_sigmas = []
			pr_string_list = []
			for r in range(pow_2_n):
				mask = masks[r]
				# to record oscillators in the next state
				next_sigma = numpy.zeros((t,), dtype = int)
				# to record how many flashed in each state				
				flashed = numpy.zeros((t,), dtype = int)
				for i in range(t, 0, -1):
					# calculate how many actually fashed after message loss					
					flashed[i - 1] = sigma[i - 1] - sum(mask[:sigma[i - 1]])
					mask = mask[sigma[i - 1]:]
				alpha = 0
				# calculate where all the oscillators will move to
				for i in range(t - 1, -1, -1):
					s = i + 1
					if s > rp:
						for k in range(t, i, -1):
							s_primed = (k % t) + 1
							if min_alphas[s][s_primed] <= alpha <= max_alphas[s][s_primed]:
								next_sigma[s_primed - 1] += sigma[s - 1]
								if s_primed == 1:
									alpha += flashed[s - 1]
									break
					else:
						# in the refractory period so move to the next state
						next = (i + 1) % t
						next_sigma[next] += sigma[s - 1]
				index = [i for i, j in enumerate(next_sigmas) if numpy.array_equal(next_sigma, j)]
				if not index:
					next_sigmas.append(next_sigma)
					pr_string_list.append(pr_strings[r])
				else:
					pr_string_list[index[0]] += (' + %s' % pr_strings[r])
			l = len(next_sigmas)
			for i in range(l):
				next_sigma = next_sigmas[i]
				if l > 1:
					writeln('(%s):' % pr_string_list[i], tabs = 2)
				write('(n_1\' = %d)' % next_sigma[0], tabs = 3)
				for j in range(2, t + 1):
					write(' & (n_%d\' = %d)' % (j, next_sigma[j - 1]))
				writeln(' +' if i < l - 1 else ';')
			
#===============================================================================
# build the transitions for when alpha_T = 0
#===============================================================================	
def	build_update_oscillator_phases_no_flashes(n, t, epsilon, rp):
	if not SKIP_NON_FLASHING_STATES:
		write('[step] (assigned & n_%d = 0) -> ' % t, tabs = 1)
		for i in range(2, t + 1):
			write('(n_%d\' = n_%d) & ' % (i, i - 1))
		writeln('(n_1\' = n_%d);' % t)
	else:
		for i in range(1, t):
			write('[step] (assigned & n_%d = 0 & n_%d > 0' % (t, i), tabs = 1)
			for j in range(i + 1, t):
				write(' & n_%d = 0' % j)
			write(') -> ')
			for k in range(1, t - (i - 1)):
				write('(n_%d\' = 0) & ' % k);				
			for k in range((t - i) + 1, t):
				write('(n_%d\' = n_%d) & ' % (k, k - (t - i)))
			writeln('(n_%d\' = n_%d);' % (t, t - (t - i)))			
	empty_lines(1)

#===============================================================================
# build the initial state indicator module
#===============================================================================		
def build_initial_state_indicator_module(n, t, epsilon, rp):
	build_comment_lines(['initial state indicator module'], '-')
	writeln('module initial_state_indicator')
	writeln('is_initial_state : bool init false;', tabs = 1)
	writeln('is_non_initial_state : bool init false;', tabs = 1)
	empty_lines(1)
	writeln('[step] (!assigned) ->', tabs = 1)
	writeln('(is_initial_state\' = true) &', tabs = 2)
	writeln('(is_non_initial_state\' = false);', tabs = 2)
	writeln('[step] (assigned) ->', tabs = 1)
	writeln('(is_initial_state\' = false) &', tabs = 2)
	writeln('(is_non_initial_state\' = is_initial_state | is_non_initial_state);', tabs = 2)
	writeln('endmodule')
	empty_lines(1)

#===============================================================================
# builds the initial transitions
#===============================================================================	       
def build_initial_transitions(n, t, epsilon, rp):
	writeln('[step] (!assigned) ->', tabs = 1)
	num_tuples = 0;
	for tuple in combinations_with_replacement_counts(t, n):
		num_tuples = num_tuples + 1;
	it_count = 0;	
	for tuple in combinations_with_replacement_counts(t, n):
		it_count = it_count + 1;
		sigma_n_i = 0;		
		binomial_product = 1;
		for i in range(t):
			binomial_product *= binomials[n - sigma_n_i, tuple[i]];
			sigma_n_i = sigma_n_i + tuple[i];
		write('(%d / pow(%d, %d)): ' % (binomial_product, t, n), tabs = 2);
		for i in range(t - 1):
			write('(n_%d\' = %d) & ' % (i + 1, tuple[i]))
		if it_count == num_tuples:
			writeln('(n_%d\' = %d);' % (t, tuple[t - 1]))
		else:
			writeln('(n_%d\' = %d) +' % (t, tuple[t - 1]))
	empty_lines(1)
			
#===============================================================================
# create a directory if it does not already exist
#===============================================================================
def create_dir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

#===============================================================================
# generate a sequence from v_min to v_max with v_step of type v_type
#===============================================================================
def get_sequence(v_min, v_max, v_step, v_type):
	return numpy.linspace(v_min, v_max, num = round((v_max - v_min) / (1 if v_step == 0 else v_step)) + 1, dtype = v_type)

#===============================================================================
# runs the process with the params and returns the output
#===============================================================================
def run_process(process_name, params):
	process = subprocess.Popen([process_name, *params], stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
	stdout = ''
	if VERBOSE:
		print('%s' % ('_' * 80))
		print('%sPRISM output' % (' ' * 34))
		print('%s' % ('_' * 80))
		while process.poll() is None:
			line = process.stdout.readline().decode('utf-8')
			if line:
				sys.stdout.write(line)
				stdout += line
		print('%s' % ('_' * 80))
	else:
		stdout, stderr = process.communicate()
		stdout = stdout.decode('utf-8')
	return stdout

#===============================================================================
# calculate all binomials up to n + t choose n + t
#===============================================================================
def calculate_binomials(n, t):
	global binomials
	factorials = [1]
	for i in range(1, n + t + 1):
		factorials.append(i * factorials[i - 1])
	binomials = numpy.empty([n + t + 1, n + t + 1], dtype = int)
	binomials[0][0] = 1
	for i in range(1, n + t + 1):
		for j in range(0, i + 1):
			binomials[i][j] = factorials[i] / (factorials[j] * factorials[i - j])
	return binomials

#===============================================================================
# return all permutations of n unlaballed balls into r boxes
#===============================================================================	
def combinations_with_replacement_counts(n, r):
   size = n + r - 1
   for indices in itertools.combinations(range(size), n-1):
       starts = [0] + [index+1 for index in indices]
       stops = indices + (size,)
       yield tuple(map(operator.sub, stops, starts))

#===============================================================================
# write n empty lines to the active file
#===============================================================================
def empty_lines(n):
	for i in range(0, n):
		writeln('');

#===============================================================================
# set the file to write to
#===============================================================================
def set_write_file(f):
	global write_file
	write_file = f

#===============================================================================
# write a line to the active file, then new line
#===============================================================================	
def writeln(string, tabs = 0):
	write(string + '\n', tabs = tabs)

#===============================================================================
# write a line to the active file, no new line
#===============================================================================	
def write(string, tabs = 0):
	write_file.write((tabs * '\t') + string)
	
#===============================================================================
# add some comments surrounded by dividers
#===============================================================================	
def build_comment_lines(lines, divider_char, tabs = 0):
	build_divider(divider_char, tabs = tabs)
	for line in lines:
		writeln('// %s' % line, tabs = tabs)
	build_divider(divider_char, tabs = tabs)
	
#===============================================================================
# build a divider with the given character
#===============================================================================	
def build_divider(divider_char, tabs = 0):
	writeln('//' + (divider_char * 78), tabs = tabs)

#===============================================================================
# override default python3 rounding and always round 0.5 up
#===============================================================================
def round_half_up(x):
	return math.ceil(x) if x - math.floor(x) >= 0.5 else math.floor(x)
	
if __name__ == '__main__':
#	build_alphas(5, 10, 0.1)
#	print(min_alphas[9][1])
	#print(max_alphas[9][1])
	#sys.exit(0)
	main(sys.argv[1:])

