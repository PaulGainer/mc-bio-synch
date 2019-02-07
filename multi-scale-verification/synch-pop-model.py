#!/usr/bin/env python3
################################################################################
# synchronisation of pulse-coupled oscillators
#     - PRISM model generation script
#
#
# authors	: Paul Gainer, Sven Linker
# email		: p.gainer@liverpool.ac.uk, s.linker@liverpool.ac.uk
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
# 05/02/17 - fixed rounding error when generating parameter sequences
# 15/02/17 - removed order parameter constants/formulae & some reward structures
# 24/04/17 - implemented more efficient failure vector calculations
# 25/04/17 - added U parameter
# 28/04/17 - improved recursive calcuations for update/fire/alpha
# 08/05/17 - added clock drift parameter
# 09/05/17 - added command line switch list to property tuples
# 10/05/17 - convert rationals to floats before writing to csv
# 11/05/17 - added power consumption reward structure
# 31/05/17 - added flag to use non-determinism for initial state selection
# 24/08/17 - fixed generation of initial states for u > 0
# 30/08/17 - added prime factorisation of large binomial coefficients
################################################################################

import errno
import itertools
import math
import numpy
from scipy.special import binom
from decimal import Decimal
import operator
import os
import re
import subprocess
import sys
import time
import signal
from sympy import factorint

#-------------------------------------------------------------------------------
# reward structure constants
#-------------------------------------------------------------------------------
R_TIME_TO_SYNCH = 0
R_FIRE_COUNT = 1
R_POWER_CONSUMPTION = 2


################################################################################
# MICAz Wireless Measurement System
################################################################################
VOLTAGE = 3.0		# 2.7V - 3.3V
DRAW_IDLE_AMPS = 0.000020
DRAW_TRANSMIT_AMPS = 0.011
DRAW_RECEIVE_AMPS = 0.0197
CYCLE_LENGTH_SECONDS = 10.0
################################################################################

'''
################################################################################
# Analog Devices ADF7023 High Performance Low Power Transceiver IC
################################################################################
VOLTAGE = 2.9		# 2.2V - 3.6V
DRAW_IDLE_AMPS = 0.00000128
DRAW_TRANSMIT_AMPS = 0.0241
DRAW_RECEIVE_AMPS = 0.0128
CYCLE_LENGTH_SECONDS = 10.0
################################################################################
'''

################################################################################
# Mirollo & Strogatz synchronisation model
################################################################################
MODEL_NAME = 'mirollo-strogatz'
#===============================================================================
# value = [value_min, value_max, value_step]
#
# Other parameter values can be used to define the values for any parameter,
# for instance, RP = [0, 't - 2', 1] will generate models with # refractory
# periods from 0 up to T - 2 etc. When referring to other parameter values put
# the expression in a string and use lowercase letters, so n, t, epsilon etc....
#===============================================================================
N = [8, 8, 1]
T = [10, 10, 1]
EPSILON = [0.1, 0.1, 0.1]
RP = [5, 5, 1]
ML = [0.1, 0.1, 0.1]
U = [0.1, 0.1, 0.1]
CD = [0, 0, 0]
REWARDS = [R_TIME_TO_SYNCH] #R_POWER_CONSUMPTION
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
def create_properties(n, t, rp, epsilon, ml, u, cd):
	properties = [('P=? [F synchronised]', ['-hybrid'])]
	properties.append(('R{"time_to_synch"}=? [F synchronised]', ['-hybrid']))
	
	#properties = [('filter(avg, R{"power_consumption"}=? [F order_parameter >= 1], is_initial_state)', ['-hybrid', '-maxiters', '100000'])]
	'''
	op = 0.1
	for i in range(10):
		properties.append(('filter(avg, R{"power_consumption"}=? [F order_parameter >= %.3g], is_initial_state)' % op, ['-hybrid', '-maxiters', '100000']))
		op += 0.1
	op = 0.1
	for i in range(10):
		properties.append(('filter(max, R{"power_consumption"}=? [F order_parameter >= %.3g], is_initial_state)' % op, ['-hybrid', '-maxiters', '100000']))
		op += 0.1
	'''		
	return properties
################################################################################


'''
################################################################################
# Mean Phase synchronisation model
################################################################################
MODEL_NAME = 'mean-phase'
#===============================================================================
# value = [value_min, value_max, value_step]
#
# Other parameter values can be used to define the values for any parameter,
# for instance, RP = [0, 't - 2', 1] will generate models with # refractory
# periods from 0 up to T - 2 etc. When referring to other parameter values put
# the expression in a string and use lowercase letters, so n, t, epsilon etc....
#===============================================================================
N = [2, 2, 1]
T = [50, 50, 2]
EPSILON = [0.1, 0.5, 0.05]
RP = [10, 24, 14]
ML = [0.0, 0.0, 0.05]
U = [0, 0, 1]
CD = [0, 0, 0]
REWARDS = [R_TIME_TO_SYNCH, R_FIRE_COUNT, R_POWER_CONSUMPTION]
#===============================================================================
# evol(x)
#===============================================================================
def evolution(x):
	return x + 1
#===============================================================================
# pert(x)
#===============================================================================
def perturbation(t, x, epsilon, alpha):
	alpha = min(alpha, 1)
	return ((x + (((2 ** alpha) - 1) * t)) / (2 ** alpha)) - x
#===============================================================================
# properties to check in the models
#===============================================================================
def create_properties(n, t, rp, epsilon, ml, u, cd):
	properties = [('P=? [F synchronised]', ['-hybrid'])]
	return properties
################################################################################	
'''

#-------------------------------------------------------------------------------
# prism configuration
#-------------------------------------------------------------------------------
PRISM = 'prism'
PRISM_MEM_PARAMS = ['-cuddmaxmem', '4g', '-javamaxmem', '12g']
TIME_COMMAND = '/usr/bin/time'
TIME_PARAMS = '--verbose'

#-------------------------------------------------------------------------------
# set to True to use non-determinism to select initial states
#-------------------------------------------------------------------------------
INITIAL_NON_DETERMINISM = False

#-------------------------------------------------------------------------------
# set to True to use prime factorisation for large binomial coefficients
#
#	- used to overcome a limitation of PRISM version 4.3.1 where 64-bit integer
#	  literals are not recognised by the parser
#-------------------------------------------------------------------------------
USE_PRIME_FACTORISATION = True
PRIME_FACTORISATION_THRESHOLD = 2147483647 # (2^31) - 1

#-------------------------------------------------------------------------------
# set to True to write debug info to terminal
#-------------------------------------------------------------------------------
DEBUG = False

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
REMOVE_GENERATED_MODELS = False

#-------------------------------------------------------------------------------
# value used to represent \star in failure vector
#-------------------------------------------------------------------------------
F_VEC_NON_FIRED = -1

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
CSV_RESULTS_TERMS = ['Result']
#CSV_RESULTS_TERMS = ['Result', 'Time for model construction', 'States', 'Transitions', 'Transition matrix', 'Time for model checking', 'Maximum resident set size']

#-------------------------------------------------------------------------------
# global variables
#-------------------------------------------------------------------------------
full_results_file = None
csv_results_file = None
prism_model_file = None
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
U_VALUES = 'get_sequence(%s, int)' % ','.join(map(str, U))
CD_VALUES = 'get_sequence(%s, float)' % ','.join(map(str, CD))

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
		for u in eval(U_VALUES):
			for t in eval(T_VALUES):
				for rp in eval(RP_VALUES):
					for epsilon in eval(EPSILON_VALUES):
						for cd in eval(CD_VALUES):
							build_alphas(n, t, epsilon)
							# generate the prism model...
							PRISM_FILE_NAME = '%smodel-%s-%d-%d-%g-%d-%d-%r.prism' % (MODELS_DIR, MODEL_NAME, n, t, epsilon, rp, u, (cd != 0))
							generate_prism_model_file(PRISM_FILE_NAME, n, t, epsilon, rp, u, cd)
							for ml in eval(ML_VALUES):
								# ...and check the properties
								check_properties_in_prism(PRISM_FILE_NAME, n, t, epsilon, rp, ml, u, cd)
							if REMOVE_GENERATED_MODELS:
								os.remove(PRISM_FILE_NAME)
	# close the results files
	full_results_file.close()
	csv_results_file.close()
	
#===============================================================================
# check all properties in the given prism model
#===============================================================================
def check_properties_in_prism(prism_file_name, n, t, epsilon, rp, ml, u, cd):
	set_write_file(full_results_file)
	writeln('=====================================================================')
	writeln('N=%d T=%d epsilon=%g RP=%d ML=%g U=%d CD=%g' % (n, t, epsilon, rp, ml, u, cd))
	set_write_file(csv_results_file)
	write('%d,%d,%g,%d,%g,%d,%g' % (n, t, epsilon, rp, ml, u, cd))
	print('for parameters: N=%d, T=%d, epsilon=%g, RP=%d, ML=%g U=%d CD=%g' % (n, t, epsilon, rp, ml, u, cd))
	for property, params in create_properties(n, t, rp, epsilon, ml, u, cd):
		print('\tverifying property: %s' % property)
		# time the prism process, we want some memory consumption info from time
		results = run_process(TIME_COMMAND, [TIME_PARAMS, PRISM, prism_file_name, '-pf', property, *PRISM_MEM_PARAMS, *params, '-const', 'ML=%g,CD=%g' % (ml, cd)]).split('\n')
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
					num_string = re.search('[0-9]+([0-9.E-]|\/)*|Infinity', line).group() 
					if '/' in num_string:
						numerator, denominator = num_string.split('/')
						num_string = str(Decimal(numerator) / Decimal(denominator))
					write(',' + num_string)
					found = True
					break
			if not found:
				write(',')
	set_write_file(csv_results_file)
	writeln('')										

#===============================================================================
# generates a prism model given a file name and parameters
#===============================================================================
def generate_prism_model_file(prism_file_name, n, t, epsilon, rp, u, cd):
	if os.path.exists(prism_file_name):
		print('using preexisting model file: %s' % prism_file_name)
		return
	prism_model_file = open(prism_file_name, 'w')
	set_write_file(prism_model_file)
	print('building model file: %s' % prism_file_name)
	# add some info about built date/time and parameters used
	build_comment_lines(
		[
			'file generated on %s at %s' % (DATE, TIME),
			'',
			'model name: %s' % MODEL_NAME
		], "=")
	if INITIAL_NON_DETERMINISM:
		writeln('mdp')
	else:
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
			'U       = %d' % u,
			'CD      = %r' % (cd > 0),
		], "-")
	empty_lines(1)
	build_comment_lines(['rewards'], '=')
	build_rewards(n, t, epsilon, rp, cd)
	build_comment_lines(['order parameter'], '=')
	build_order_parameter(n, t, epsilon, rp)
	build_comment_lines(['constants'], '=')
	build_constants(n, t, epsilon, rp, cd)
	build_comment_lines(['formulae'], '=')
	build_formulae(n, t, epsilon, rp)
	build_comment_lines(['modules'], '=')
	build_initial_state_indicator_module(n, t, epsilon, rp)
	if cd > 0:
		build_clock_drift_module()
	build_oscillator_population_module(n, t, epsilon, rp, u, cd)
	prism_model_file.close()

#===============================================================================
# build the reward structures
#===============================================================================	
def build_rewards(n, t, epsilon, rp, cd):
	if R_TIME_TO_SYNCH in REWARDS:
		#---------------------------------------------------------------------------
		# time to synch reward
		#---------------------------------------------------------------------------
		build_comment_lines(['records cycles until synchrony'], '-')
		writeln('rewards "time_to_synch"')
		if (not SKIP_NON_FLASHING_STATES) or (cd > 0):
			writeln('%sassigned & (!synchronised) : 1 / %d;' % ('!drift & ' if cd > 0 else '', t), tabs = 1)
		else:
			writeln('assigned & (n_%d > 0) & (!synchronised) : 1 / %d;' % (t, t), tabs = 1)
			for i in range(1, t):
				write('assigned & (n_%d = 0 & n_%d > 0' % (t, i), tabs = 1)
				for j in range(i + 1, t):
					write(' & n_%d = 0' % j)
				writeln(') & (!synchronised) : %d / %d;' % (t - i, t))
		writeln('endrewards')
		empty_lines(1)
	if R_FIRE_COUNT in REWARDS:
		#---------------------------------------------------------------------------
		# fire count reward	
		#---------------------------------------------------------------------------
		build_comment_lines(['records number of oscillator firings'], '-')
		writeln('rewards "fire_count"')
		writeln('%sis_non_initial_state : n_1;' % ('!drift & ' if cd > 0 else ''), tabs = 1)
		writeln('endrewards')
		empty_lines(1)
	if R_POWER_CONSUMPTION in REWARDS:
		#---------------------------------------------------------------------------
		# power consumption reward	
		#---------------------------------------------------------------------------
		build_comment_lines(['records power consumption in watt hours'], '-')
		writeln('rewards "power_consumption"')
		time_in_hours = CYCLE_LENGTH_SECONDS / t / 3600.0
		transmit_time_in_hours = 0.001 / 3600.0
		draw_idle_watt_hours = DRAW_IDLE_AMPS * VOLTAGE * time_in_hours
		draw_transmit_watt_hours = DRAW_TRANSMIT_AMPS * VOLTAGE * time_in_hours
		draw_receive_watt_hours = DRAW_RECEIVE_AMPS * VOLTAGE * time_in_hours
		if (not SKIP_NON_FLASHING_STATES) or (cd > 0):
			writeln('%sassigned :' % ('(!drift) & ' if cd > 0 else ''), tabs = 1)
			write('(n_%d * %.10f)' % (1, transmit_time_in_hours), tabs = 2)
			if rp > 0:
				writeln('')
				write('+ ((n_1', tabs = 2)
				for i in range(rp - 1):
					write(' + n_%d' % (i + 2))
				write(') * %.10f)' % (draw_idle_watt_hours))
			if rp < t - 1:
				writeln('')
				write('+ ((n_%d' % (t - 1), tabs = 2)
				for i in range(t - rp - 2):
					write(' + n_%d' % ((t - 2) - i))
				write(') * %.10f)' % (draw_receive_watt_hours))
			writeln(';')
		else:
			writeln('%sassigned & (n_%d > 0) :' % ('(!drift) & ' if cd > 0 else '', t), tabs = 1)
			write('(n_%d * %.10f)' % (1, transmit_time_in_hours), tabs = 2)
			if rp > 0:
				writeln('')
				write('+ ((n_1', tabs = 2)
				for i in range(rp - 1):
					write(' + n_%d' % (i + 2))
				write(') * %.10f)' % (draw_idle_watt_hours))
			if rp < t - 1:
				writeln('')
				write('+ ((n_%d' % (t - 1), tabs = 2)
				for i in range(t - rp - 2):
					write(' + n_%d' % ((t - 2) - i))
				write(') * %.10f)' % (draw_receive_watt_hours))
			writeln(';')
			for i in range(1, t):
				shift = t - i
				write('%sassigned & (n_%d = 0 & n_%d > 0' % ('(!drift) & ' if cd > 0 else '', t, i), tabs = 1)				
				for j in range(i + 1, t):
					write(' & n_%d = 0' % j)
				write(') : ')
				writeln('')
				write('(n_%d * %.10f)' % (1, transmit_time_in_hours), tabs = 2)
				if rp > 0:
					writeln('')
					write('+ ((', tabs = 2)
					for j in range(t - shift):
						phi = j + 1
						idle_len = min(max((rp + 1) - phi, 0), shift)
						write('(%d * n_%d)' % (idle_len, phi))
						if j < t - shift - 1:
							write(' + ')
					write(') * %.10f)' % draw_idle_watt_hours)
				if rp < t - 1:
					writeln('')
					write('+ ((', tabs = 2)
					for j in range(t - shift):
						phi = j + 1
						idle_len = min(max((rp + 1) - phi, 0), shift)
						receive_len = shift - idle_len
						write('(%d * n_%d)' % (receive_len, phi))
						if j < t - shift - 1:
							write(' + ')
					write(') * %.10f)' % draw_receive_watt_hours)
				writeln(';')					
		writeln('endrewards')
		empty_lines(1)				
	
#===============================================================================
# build the constants and forumlae for the order parameter
#===============================================================================
def build_order_parameter(n, t, epsilon, rp):
	# calculate x and y components of unit vectors
	for i in range(0, t):
		write('const double unit_vector_x_%d = %.12f; ' % (i + 1, math.cos(2 * math.pi * float((1 / t) * i))))
		writeln('const double unit_vector_y_%d = %.12f;' % (i + 1, math.sin(2 * math.pi * float((1 / t) * i))))
	empty_lines(1)
	write('formula unit_vector_x_avg_squared = pow(((')
	for i in range(1, t + 1):
		write('(unit_vector_x_%d * n_%d)%s' % (i, i, '' if i == t else ' + ' ))
	writeln(') / %d), 2);' % n)
	write('formula unit_vector_y_avg_squared = pow(((')
	for i in range(1, t + 1):
		write('(unit_vector_y_%d * n_%d)%s' % (i, i, '' if i == t else ' + '))
	writeln(') / %d), 2);' % n)
	writeln('formula order_parameter = pow(unit_vector_x_avg_squared + unit_vector_y_avg_squared, 0.5);')
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
	writeln('formula assigned = num_oscillators = %d;' % n)
	empty_lines(1)
	#---------------------------------------------------------------------------
	# synchronised formula
	#---------------------------------------------------------------------------
	build_comment_lines(['true if all oscillators have the same phase'], '-')
	write('formula synchronised = ')
	for i in range(1, t):
		write('n_%d = %d | ' % (i, n))
	writeln('n_%d = %d;' % (t, n))
	empty_lines(1)
	#---------------------------------------------------------------------------
	# message loss probability formulae
	#---------------------------------------------------------------------------
	build_comment_lines(['message loss probabilities'], '-')
	writeln('formula pr_fail_0_0 = 1;')
	for i in range(n):
		for j in range(i + 2):
			bc = binom(i + 1, j)
			if USE_PRIME_FACTORISATION and bc > PRIME_FACTORISATION_THRESHOLD:
				prime_string = ''
				prime_dict = factorint(int(bc))
				for prime in prime_dict:
					if prime > PRIME_FACTORISATION_THRESHOLD:
						error_string = 'Prime factorisation failed: a factor ({}) of {} was greater than {}'.format(bc, PRIME_FACTORISATION_THRESHOLD, prime)
						raise Exception(error_string)
				prime_power_pairs = list(prime_dict.items())
				num_pairs = len(prime_power_pairs)
				for p in range(num_pairs):
					k, v = prime_power_pairs[p]
					prime_string += 'pow({}, {})'.format(k, v)
					if p < num_pairs - 1:
						prime_string += ' * '						
				writeln('formula pr_fail_{}_{} = (pow(ML, {}) * pow(1 - ML, {}) * {});'.format(i + 1, j, j, (i + 1) - j, prime_string))
			else:
				writeln('formula pr_fail_{}_{} = (pow(ML, {}) * pow(1 - ML, {}) * {});'.format(i + 1, j, j, (i + 1) - j, binom(i + 1, j)))			
	empty_lines(1)
						
#===============================================================================
# build the constants
#===============================================================================
def build_constants(n, t, epsilon, rp, cd):
	#---------------------------------------------------------------------------
	# infinity constant
	#---------------------------------------------------------------------------
	build_comment_lines(['infinity'], '-')
	writeln('const int INFINITY = %d;' % (n + 1))
	empty_lines(1)
	#---------------------------------------------------------------------------
	# message loss constant
	#---------------------------------------------------------------------------
	build_comment_lines(['message loss probability'], '-')
	writeln('const double ML;')
	empty_lines(1)
	#---------------------------------------------------------------------------
	# clock drift constant
	#---------------------------------------------------------------------------
	if cd > 0:
		build_comment_lines(['clock drift probability'], '-')
		writeln('const double CD;')
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
			min_alphas[i][1] = minimum
			max_alphas[i][1] = n
			l = l + [(i, 1)]
		else:
			min_alphas[i][s] = minimum
			max_alphas[i][s] = maximum
			l = l+ [(i, s)]
	for i in range (1, t + 1):
		for j in range(1, t + 1):
			if not ((i, j) in l) and (j > i or j == 1):
				if i == t:
					min_alphas[i][j] = 0
				else:
					min_alphas[i][j] = n + 1
				max_alphas[i][j] = n + 1
					
#===============================================================================
# build the oscillator population module
#===============================================================================
def build_oscillator_population_module(n, t, epsilon, rp, u, cd):
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
	# clock drift transitions
	#---------------------------------------------------------------------------
	if cd > 0:
		build_clock_drift_transitions()
	#---------------------------------------------------------------------------
	# transitions from starting state
	#---------------------------------------------------------------------------
	build_comment_lines(['transitions from starting state'], '-', tabs = 1)
	build_initial_transitions(n, t, epsilon, rp, u)
	build_comment_lines(['update oscillator phases, no flashes'], '-', tabs = 1)
	build_update_oscillator_phases_no_flashes(n, t, epsilon, rp, u, cd)
	build_comment_lines(['update oscillator phases, flashes'], '-', tabs = 1)
	build_update_oscillator_phases_flashes(n, t, epsilon, rp, u, cd)
	writeln('endmodule')

def build_clock_drift_transitions():
	build_comment_lines(['clock drift transitions'], '-', tabs = 1)	
	writeln('[step] (assigned & drift & (n_1 > 0)) ->', tabs = 1)
	writeln('(CD): (n_1\' = n_1 - 1) & (n_2\' = n_2 + 1) + (1 - CD): true;', tabs = 2)
	writeln('[step] (assigned & drift & (n_1 = 0)) -> true;', tabs = 1)
	empty_lines(1)

#===============================================================================
# update function
#===============================================================================
def update(n, t, epsilon, rp, phi, g_state, f_vec, last_alpha, last_fire):
	alpha_phi = alpha(n, t, epsilon, rp, phi, g_state, f_vec, last_alpha, last_fire)
	if phi <= rp:
		next_phi = evolution(phi)
	else:
		next_phi = evolution(phi) + int(round_half_up(perturbation(t, phi, epsilon, alpha_phi)))
	return (next_phi, alpha_phi)

#===============================================================================
# fire predicate
#===============================================================================
def fire(n, t, epsilon, rp, phi, g_state, f_vec, last_alpha, last_fire):
	next_phi, alpha_phi = update(n, t, epsilon, rp, phi, g_state, f_vec, last_alpha, last_fire)
	return (next_phi > t, alpha_phi)

#===============================================================================
# alpha value calculation
#===============================================================================
def alpha(n, t, epsilon, rp, phi, g_state, f_vec, last_alpha, last_fire):
	if phi == t:
		return 0
	elif f_vec[phi] != F_VEC_NON_FIRED and last_fire: # fire(n, t, epsilon, rp, phi + 1, g_state, f_vec):
		return last_alpha + g_state[phi] - f_vec[phi]
	else:
		return last_alpha

#===============================================================================
# the transition function
#===============================================================================
def tau(n, t, epsilon, rp, phi, g_state, f_vec, last_alpha, last_fire):
	fire_phi, alpha_phi = fire(n, t, epsilon, rp, phi, g_state, f_vec, last_alpha, last_fire)
	if fire_phi:
		return (fire_phi, alpha_phi, 1)
	else:
		next_phi, alpha_phi = update(n, t, epsilon, rp, phi, g_state, f_vec, last_alpha, last_fire)
		return (fire_phi, alpha_phi, next_phi)

#===============================================================================
# the successor function
#===============================================================================
def succ(n, t, epsilon, rp, g_state, f_vec, last_alpha, last_fire):
	next_state = numpy.zeros((t,), dtype = int)
	for phi in reversed(range(t)):
		last_fire, last_alpha, updated_phi = tau(n, t, epsilon, rp, phi + 1, g_state, f_vec, last_alpha, last_fire)
		next_state[updated_phi - 1] += g_state[phi]
	return next_state

f_vec_count = 1
num_global_firing_states = 0
global_firing_state = 1
#===============================================================================
# builds the set of possible failure vectors for a global state
#===============================================================================
def build_failure_vectors(n, t, epsilon, rp, phi, g_state, partial_f_vec, last_alpha, last_fire):
	global f_vec_count
	if not partial_f_vec:
		f_vec_count = 1
	f_vec_list = []
	f_vec_phi_g_1_f = [F_VEC_NON_FIRED] * phi
	f_vec_phi_g_1_f.extend(partial_f_vec)
	fire_phi, alpha_phi = fire(n, t, epsilon, rp, phi, g_state, f_vec_phi_g_1_f, last_alpha, last_fire)
	if phi > 1 and fire_phi:
		for k in range(g_state[phi - 1] + 1):
			f_vec_prefix_k = [k]
			f_vec_prefix_k.extend(partial_f_vec)
			f_vec_list.extend(build_failure_vectors(n, t, epsilon, rp, phi - 1, g_state, f_vec_prefix_k, alpha_phi, fire_phi))
	elif phi == 1:
		f_vec_prefix_non_fired = [F_VEC_NON_FIRED]
		f_vec_prefix_non_fired.extend(partial_f_vec)
		fire_phi, alpha_phi = fire(n, t, epsilon, rp, 1, g_state, f_vec_prefix_non_fired, last_alpha, last_fire)
		if fire_phi:
			for k in range(g_state[0] + 1):
				f_vec_prefix_k = [k]
				f_vec_prefix_k.extend(partial_f_vec)
				f_vec_list.extend([tuple(f_vec_prefix_k)])
				if DEBUG:
					print('built failure vector(%d) for global state (%d/%d)' % (f_vec_count, global_firing_state, num_global_firing_states))
				f_vec_count += 1				
	if not fire_phi:
		f_vec = [F_VEC_NON_FIRED] * phi
		f_vec.extend(partial_f_vec)
		f_vec_list.extend([tuple(f_vec)])
		if DEBUG:
			print('built failure vector(%d) for global state (%d/%d)' % (f_vec_count, global_firing_state, num_global_firing_states))
		f_vec_count += 1				
	return f_vec_list

#===============================================================================
# returns the probability of b broadcast failures given n messages, as a string
#===============================================================================
def p_fail(n, b):
	return ('pr_fail_%d_%d' % (n, b))
	
#===============================================================================
# returns the probability of a transition being made w.r.t. a failure vector
#===============================================================================
def p_tau(n, t, epsilon, rp, g_state, f_vec):
	product = '('
	is_one = False
	for i in range(t):
		if f_vec[i] != F_VEC_NON_FIRED:
			product += p_fail(g_state[i], f_vec[i])
			is_one = True
		if i < t - 1 and is_one:
			product += ' * '
	if not is_one:
		product += '1'
	product += ')'
	return product

#===============================================================================
# build the transitions for when alpha_T > 0
#===============================================================================
def build_update_oscillator_phases_flashes(n, t, epsilon, rp, u, cd):
	global num_global_firing_states, global_firing_state
	global_firing_state = 1
	global_firing_states = generate_global_firing_states(n, t, u)
	for state in global_firing_states:
		write('[step] (%sassigned' % ('!drift & ' if cd > 0 else ''), tabs = 1)
		for i in range(t):
			write(' & n_%d = %d' % (i + 1, state[i]))
		writeln(') ->')
		next_states = []
		pr_string_list = []
		for f_vec in build_failure_vectors(n, t, epsilon, rp, t, state, [], 0, True):
			next_state = succ(n, t, epsilon, rp, state, f_vec, 0, True)
			index = [i for i, j in enumerate(next_states) if numpy.array_equal(next_state, j)]
			if not index:
				next_states.append(next_state)
				pr_string_list.append(p_tau(n, t, epsilon, rp, state, f_vec))
			else:
				pr_string_list[index[0]] += (' + %s' % p_tau(n, t, epsilon, rp, state, f_vec))
		l = len(next_states)
		for i in range(l):
			next_state = next_states[i]
			if l > 1:
				writeln('(%s):' % pr_string_list[i], tabs = 2)
			write('(n_1\' = %d)' % next_state[0], tabs = 3)
			for j in range(2, t + 1):
				write(' & (n_%d\' = %d)' % (j, next_state[j - 1]))
			writeln(' +' if i < l - 1 else ';')
		global_firing_state += 1
		
#===============================================================================
# build the transitions for when alpha_T = 0
#===============================================================================	
def	build_update_oscillator_phases_no_flashes(n, t, epsilon, rp, u, cd):
	if (not SKIP_NON_FLASHING_STATES) or (cd > 0):
		write('[step] (%sassigned & n_%d = 0) -> ' % ('!drift & ' if cd > 0 else '', t), tabs = 1)
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
				write('(n_%d\' = 0) & ' % k)
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
# build the clock drift module
#===============================================================================		
def build_clock_drift_module():
	build_comment_lines(['clock_drift module'], '-')
	writeln('module clock_drift')
	writeln('drift : bool init false;', tabs = 1)
	empty_lines(1)
	writeln('[step] (!assigned) -> true;', tabs = 1)
	writeln('[step] (assigned) -> (drift\' = !drift);', tabs = 1)
	writeln('endmodule')
	empty_lines(1)
	
#===============================================================================
# generates all possible initial states
#===============================================================================	       
def generate_initial_states(n, t, u):
	initial_states = []
	if u == 0:
		for state in combinations_with_replacement_counts(n, t):
			initial_states.append(state)
	else:
		for state in combinations_with_replacement_counts(u, t):
			for index in range(t):
				state_temp = list(state)
				state_temp[index] += n - u
				state_temp_tuple = tuple(state_temp)
				if state_temp_tuple not in initial_states:
					initial_states.append(state_temp_tuple)
	return initial_states
	
#===============================================================================
# generates global states
#===============================================================================	       
def generate_global_firing_states(n, t, u):
	global num_global_firing_states
	global_firing_states = []
	if u == 0:
		for state in combinations_with_replacement_counts(n, t):
			global_firing_states.append(state)
	else:
		for state in combinations_with_replacement_counts(u, t):
			for index in range(t):
				state_temp = list(state)
				state_temp[index] += n - u
				state_temp_tuple = tuple(state_temp)
				if state_temp_tuple not in global_firing_states:
					global_firing_states.append(state_temp_tuple)
	global_firing_states = [state for state in global_firing_states if state[t - 1] != 0]
	num_global_firing_states = len(global_firing_states)
	return global_firing_states
		
#===============================================================================
# builds the initial transitions
#===============================================================================	       
def build_initial_transitions(n, t, epsilon, rp, u):
	if INITIAL_NON_DETERMINISM:
		initial_states = generate_initial_states(n, t, u)
		num_tuples = len(initial_states)
		for state in initial_states:
			writeln('[step] (!assigned) ->', tabs = 1)
			for i in range(t - 1):
				write('(n_%d\' = %d) & ' % (i + 1, state[i]))
			writeln('(n_%d\' = %d);' % (t, state[t - 1]))
	else:
		writeln('[step] (!assigned) ->', tabs = 1)
		num_tuples = 0
		initial_states = generate_initial_states(n, t, u)
		num_tuples = len(initial_states)
		binomial_products = []
		for state in initial_states:
			binomial_product = 1
			sigma_n_i = 0
			for i in range(t):
				binomial_product *= binom(n - sigma_n_i, state[i])
				sigma_n_i = sigma_n_i + state[i]				
			binomial_products.append(binomial_product)
		binomial_product_sum = sum(binomial_products)
		it_count = 0
		for state in initial_states:
			it_count = it_count + 1
			sigma_n_i = 0
			if u > 0:
				write('(%d / %d): ' % (binomial_products[it_count - 1], binomial_product_sum), tabs = 2)
			else:
				write('(%d / pow(%d, %d)): ' % (binomial_products[it_count - 1], t, n), tabs = 2)
			for i in range(t - 1):
				write('(n_%d\' = %d) & ' % (i + 1, state[i]))
			if it_count == num_tuples:
				writeln('(n_%d\' = %d);' % (t, state[t - 1]))
			else:
				writeln('(n_%d\' = %d) +' % (t, state[t - 1]))
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
# return all permutations of k unlaballed balls into n boxes
#===============================================================================	
def combinations_with_replacement_counts(k, n):
   size = n + k - 1
   for indices in itertools.combinations(range(size), n - 1):
       starts = [0] + [index + 1 for index in indices]
       stops = indices + (size,)
       yield tuple(map(operator.sub, stops, starts))

#===============================================================================
# write n empty lines to the active file
#===============================================================================
def empty_lines(n):
	for i in range(0, n):
		writeln('')

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
	
#-------------------------------------------------------------------------------
# builds a fully expanded expression for 'synchronised'
#-------------------------------------------------------------------------------
def build_synchronised_string(n, t):
	synchronised_expanded = ''
	for i in range(t - 1):
		synchronised_expanded += 'n_%d = %d | ' % (i + 1, n)
	synchronised_expanded += 'n_%d = %d' % (t, n)
	return synchronised_expanded

def signal_handler(signal, frame):
	sys.exit(0)

if __name__ == '__main__':
	signal.signal(signal.SIGINT, signal_handler)
	main(sys.argv[1:])



