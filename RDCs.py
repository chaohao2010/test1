#This script is used to calculate the RDCs

import os, re, sys
import argparse
import numpy as np
import pandas as pd
import mdtraj as md

def InputOptions():
	options = argparse.ArgumentParser(description = 'Following is detail description')
	options.add_argument('-v', '--ver', action = 'store_true', help = 'version')
	options.add_argument('-d', '--directory', required = True, type = str, help = 'working derectory')
	options.add_argument('-s', '--system', required = True, type = str, help = 'simulation system')
	options.add_argument('-f', '--forcefield', required = True, type = str , help = 'force field')
	options.add_argument('-t', '--trajectory', type = str, default = '1-5', help = 'trajectories for analysis' )
	options.add_argument('-b', '--begin', required = True, type = int, help = 'time for beginning of analysis')
	options.add_argument('-e', '--end', required = True, type = int, help = 'time for ending of analysis')
	options.add_argument('-w', '--windows', type = int, default = 15, help = 'windows size')
	options.add_argument('-overwrite', action = 'store_true', help = 'Overwrite the result file')
	options.add_argument('-bestFit', type = str, default = 'Yes', help = 'using bestFit to calculate RDCs')
	options = options.parse_args()
	return options

def TrajParse(string):
	trajlist = []
	for s in string.split(','):
		if ('-') in s:
			trajlist.extend(range(int(s.split('-')[0]), int(s.split('-')[1])+1))
		else:
			trajlist.append(int(m))
	return trajlist


'''Using the number of pdb files to estimate the simulation time
   Assume that gap is 10 when generating the pdb files from trajectory
   and there are only pdbs files in the 'pdbs' directory'''

def TotalTime(ffdir, trajlist):
	pdbs = os.popen('ls %s/%dC/pdbs/ |wc -l' % (ffdir, trajlist[0])).readlines()
	total_pdb = int(pdbs[0])
	if total_pdb % 10 == 1:
		total_pdb -= 1
		flag = 'gromacs'
	else:
		flag = 'amber'
	total_time = total_pdb / 20
	return total_time, flag

def TotalResidues(ffdir, system, trajlist, software):
	if software == 'gromacs':
		pdb = md.load_pdb('%s/%dC/pdbs/%s_1.pdb' % (ffdir, trajlist[0], system))
	elif software == 'amber':
		pdb = pd.load_pdb('%s/%dC/pdbs/%s.pdb.1' % (ffdir, trajlist[0], system))
	return pdb.n_residues
	

def RDCsCalculation(ffdir, system, trajlist, ExpRDCs, software, total_time, begin, end, total_res, windows, bestFit):
	os.chdir(ffdir)
	total_pdb = total_time * 20
	begin_pdb = begin * 20 + 1
	end_pdb = end * 20 + 1
	windows_time = total_res - windows + 1 #how many time we move the windows
	Qvalues = []
	if not os.path.exists('%s/RDCs' % ffdir):
		os.makedirs('RDCs')
		print('  Created directory \033[34mRDCs\033[0m')
	else:
		print('  Directory \033[34mRDCs\033[0m exists') 
	
	print('  starting calculating RDCs..')
	for traj in trajlist:
		for pdb in range(begin_pdb, end_pdb): 
			if software == 'gromacs':
				pdbfile = '%s/%dC/pdbs/%s_%d.pdb' % (ffdir, traj, system, pdb) #modify here for gromacs
			elif software == 'amber':
				pdbfile = '%s/%dC/pdbs/%s.pdb.%d' % (ffdir, traj, system, pdb) #modify here for amber
			if not os.path.isfile(pdbfile):
				print('No such pdb file: %s' % (pdbfile))
				exit()
			else:
				for w in range(0, windows_time):
					res_begin = w+1
					res_end = w+windows
					if(bestFit == 'YES'):
						pales_out = os.popen('pales -bestFit -inD %s -pdb %s -s1 %d -sN %d 2>/dev/null'
										% (ExpRDCs, pdbfile, res_begin, res_end)).readlines()
					elif(bestFit == 'NO'):
						pales_out = os.popen('pales -inD %s -pdb %s -s1 %d -sN %d 2>/dev/null'
										% (ExpRDCs, pdbfile, res_begin, res_end)).readlines()
					for line in pales_out:
						pales_match = re.match(r'DATA Q SAUPE\s+(\S+)', line)
						if pales_match:
							Qvalues.append(float(pales_match.group(1)))
		print('    traj %d finished' % traj)
	#print(len(Qvalues))
	Qvalues_df = pd.DataFrame({'Qvalues':Qvalues})
	Qavg = Qvalues_df['Qvalues'].mean() #average values
	Qstd = Qvalues_df['Qvalues'].std()	#standard deviation
	Qsem = Qstd / np.sqrt(Qvalues_df['Qvalues'].count()) #standard error of mean
	return Qavg, Qstd, Qsem
	#print('Average Q value: %.4f' % Qavg)
	#print('Standard error of mean: %.4f' % Qsem)

def main():
	options = InputOptions()

	# check the options
	systemdir = options.directory + '/' + options.system
	forcefielddir = systemdir + '/' + options.forcefield
	if not os.path.exists(options.directory):
		print('Error: your working directory \033[31m%s\033[0m is not exists' % options.directory)
		return False
	if not os.path.exists(systemdir):
		print('Error: can not find the simulation system \033[31m%s\033[0m' % options.system)
		return False
	if not os.path.exists(forcefielddir):
		print('Error: can not find the force field \033[31m%s\033[0m of corresponding system \033[31m%s\033[0m' % (options.forcefield, options.system))
		return False
	if re.search(r'[^\d,-]', options.trajectory):
		print('Error: unknown trajectory format \033[31m%s\033[0m' % options.trajectory)
		return False
	pattonYes = re.compile(r'\by(es)?\b', re.I)
	pattonNo  = re.compile(r'\bn(o)?\b', re.I)
	if not re.match(pattonYes, options.bestFit) and not re.match(pattonNo, options.bestFit):
		print('Error: unknown bestFit options \033[31m%s\033[0m' % options.bestFit)
		return False

	# Required Arguments
	workdir = options.directory
	system = options.system
	forcefield = options.forcefield
	begin = options.begin
	end = options.end

	# Default Arguments
	windows = options.windows
	trajlist = TrajParse(options.trajectory)
	if re.match(pattonYes, options.bestFit):
		bestFit = 'YES'
	elif re.match(pattonNo, options.bestFit):
		bestFit = 'NO'
	
	#calulated Arguments
	time, software = TotalTime(forcefielddir, trajlist)
	residues = TotalResidues(forcefielddir, system, trajlist, software) 

	print(options.overwrite)
	# Output Arguments
	print('-'*60)
	print('  Input Arguments are following')
	print('  working directory: \033[34m%s\033[0m' % workdir)
	print('  system: \033[34m%s\033[0m' % system)
	print('  force field: \033[34m%s\033[0m' % forcefield)
	print('  trajectories: \033[34m%s\033[0m' % trajlist)
	print('  begin time(ns): \033[34m%s\033[0m' % begin)
	print('  end time(ns): \033[34m%s\033[0m' % end)
	print('  aligment windows: \033[34m%s\033[0m' % windows)
	print('  simulation time: \033[34m%d ns\033[0m' % time)
	print('  simulation software: \033[34m%s\033[0m' % software)
	print('  number of residues: \033[34m%d\033[0m' % residues)
	print('  bestFit option: \033[34m%s\033[0m' % bestFit)
	print('-'*60)
	
	# RDC calculation start
	print('Processing...')
	print('  Checking the Exp. data...')
	ExpRDCs = '/lustre/home/clschf/FFDB/%s/RDC_HN' % system #TODO modify the name to "~/FFDB/...."
	#print(ExpRDCs)
	if not os.path.isfile(ExpRDCs):
		print('Error: can not find the Exp. RDCs data')
		return False
	else:
		if not os.path.isfile('%s/RDCs/%s_%s_RDCs_%d-%d.dat' % (forcefielddir, system, forcefield, begin, end)):
			RDCout = open('%s/RDCs/%s_%s_RDCs_%d-%d.dat' % (forcefielddir, system, forcefield, begin, end), 'w')
			Qavg, Qstd, Qsem = RDCsCalculation(forcefielddir, system, trajlist, ExpRDCs, software, time , begin, end, residues, windows, bestFit)
			RDCout.write('Average Q value: %.4f' % Qavg)
			RDCout.write('Standard deviation: %.4f' % Qstd)
			RDCout.write('Standard error of mean: %.4f' % Qsem)
			RDCout.close()
		else:
			if not options.overwrite:
				print('  Find file \033[34m%s_%s_RDCs_%d-%d.dat\033[0m, finish calculation' % (system, forcefield, begin, end))
			else:
				RDCout = open('%s/RDCs/%s_%s_RDCs_%d-%d.dat' % (forcefielddir, system, forcefield, begin, end), 'w')
				Qavg, Qstd, Qsem = RDCsCalculation(forcefielddir, system, trajlist, ExpRDCs, software, time , begin, end, residues, windows, bestFit)
				RDCout.write('Average Q value: %.4f' % Qavg)
				RDCout.write('Standard deviation: %.4f' % Qstd)
				RDCout.write('Standard error of mean: %.4f' % Qsem)
				RDCout.close()
		
			

if __name__ == '__main__':
	main()
