import shutil
import random
import argparse
import itertools
import os
from stiff import generate_stiff

def generate_initial_molecule():
	# probably not needed? there's a cpp script that does this, could just wrap over it lol
	pass

def generate_master_lams_file():
	# to vary bond coefficients, for future proofing
	pass

def main(do_local, do_slurm, jobname):
	# generate many scripts
	# run all in a master lammps file using the `include` command
	# batches of 20 replicas

	# folder structure
	# SMC_single_polymer/N1000/lp10/repX

	# or....... just generate 
	
	# source files:
	# original molecule file
	# parameter file? for cleanliness i guess, but otherwise not super important

	# consider: non integer persistence lengths?
	# database of runs? 

	# where is the variable info stored?
	# separate parameter file starts to make more sense since you can quite literally parse it as a fixed width file 

	# then generate bash file that will execute all the scripts

	###########################3
	#
	#
	#

	if do_local:
		path_str = "local"
	else:
		path_str = "slurm"


	# master_folder = '/home/zy/Documents/tap/smc-single-polymer/' 
	# jobname = 'stiff-14-fix-edgecase'

	start_delta = 20 # hard coded in the fix... 

	lp_stiff = [50]
	n_stiff = [1, 2, 3, 7, 10, 14]

	lp_list = [5] # list of persistence lengths to run through

	parameter_set = itertools.product(lp_list, lp_stiff, n_stiff)

	lpol = 1000 # length of polymer

	ratesmc = 100 # smc movement attempt 
	run_duration = 100000 # total number of steps

	nrep = 96 # number of replicas to do
	npara = 16 # no. of parallel jobs
	start_index = 0

	#
	#
	#
	default_paths = {
		"local":{
			"script":'/home/zy/Documents/tap/masterfile.lam',
			"molecule_file" : '/home/zy/Documents/tap/smc-lammps/data/In_conf.Nb1000.RW4.fixed.dat',
			"script_stiff":"/home/zy/Documents/tap/masterfile_stiff.lam",
			"folder":"/home/zy/Documents/tap/smc-single-polymer/",
		},
		"slurm":{
			"script":'/home/v1zchoon/masterfile.lam',
			"molecule_file" : '/home/v1zchoon/smc-lammps/data/In_conf.Nb1000.RW4.fixed.dat',
			"script_stiff":"/home/v1zchoon/masterfile_stiff.lam",
			# "folder":"/home/v1zchoon/smc-single-polymer/",
			"folder":"/scratch/smc-single-polymer",
			"slurm_output": "/home/v1zchoon/slurm-wd",
		}
	}
	###############################
	master_folder = default_paths[path_str]["folder"]


	bash_target = os.path.join(master_folder, jobname, 'run.sh')

	if not os.path.exists(master_folder):
		os.makedirs(master_folder)
	if do_slurm:
		slurm_folder = os.path.join(default_paths["slurm"]["slurm_output"], )
		if not os.path.exists(slurm_folder):
			os.makedirs(slurm_folder)
	

	lp_bash_fps = [] # collect all the different bash files to execute all at once
	run_str = ["" for _ in range(npara)]

	for c, p in enumerate(parameter_set):
		print(c)
		lp = p[0]

		if p[1] == None or p[2] == None:
			master_script = default_paths[path_str]["script"] 
			master_molecule_file = default_paths[path_str]["molecule_file"] 

			lp_folder = os.path.join(master_folder, jobname, f"N{lpol}", f"lp{lp:02}")
			_start_position = 0

		else:
			master_script = default_paths[path_str]["script_stiff"] 

			_lp_stiff = p[1]
			_n_stiff = p[2]

			master_molecule_file = generate_stiff(_n_stiff, os.path.join(master_folder, "stiff_molecules"), default_paths[path_str]["molecule_file"])

			lp_folder = os.path.join(master_folder, jobname, f"N{lpol}", f"lp{lp:02}", f"lp-stiff{_lp_stiff:02}", f"n-stiff{_n_stiff:d}")

			_start_position = int((lpol - _n_stiff)/2 - 20) # ??? e.g. for 1000 - 100, start at 430, move until 450

		parameter_fp = os.path.join(lp_folder, 'parameters.dat')

		lp_bash_fp = os.path.join(lp_folder, 'run.sh')
		lp_bash_fps.append(lp_bash_fp)

		if not os.path.exists(lp_folder):
			os.makedirs(lp_folder)

		with open(parameter_fp, 'w') as f:
			f.write("variable lpol equal {}\n".format(lpol))
			f.write("variable ratesmc equal {}\n".format(ratesmc))
			f.write("variable persistence_length equal {:.2f}\n".format(lp))
			f.write("variable run_duration equal {:d}\n".format(run_duration))

			f.write("variable stiff_persistence_length equal {:.2f}\n".format(_lp_stiff))
			f.write("variable start_position equal {:d}\n".format(_start_position))
			f.write("variable n_stiff equal {:d}\n".format(_n_stiff))
			f.write("variable max_jump equal {:d}\n".format(int(lpol/2)))

			f.write("variable noiseseed equal {:d}\n".format((random.randint(1, 32768))))
			f.write("variable smcseed equal {:d}\n".format((random.randint(1, 32768))))

			# f.write(f"angle_coeff 1 {lp:d}") # assumes that persistence length lp is an integer

		if do_local:
			with open(lp_bash_fp, 'w') as f:
				f.write(f"nparajobs={npara}\n")
				# f.write("export OMP_NUM_THREADS=8\n")
				f.write("cd {}\n".format(lp_folder))
				f.write("for i in {{{}..{}}}; do\n".format(start_index, start_index + nrep - 1))
				f.write("mkdir -p rep$i\n")
				f.write("cd rep$i\n")
				f.write("vi=$(( ${i}%${nparajobs} ))\n")

				f.write("cp -R -p -u ../parameters.dat .\n".format(parameter_fp)) # parameter file
				f.write("cp -R -p -u {} masterfile.lam\n".format(master_script)) # lams script
				f.write("cp -R -p -u {} molecule.dat\n".format(master_molecule_file))
				# f.write("r1=$(shuf -i 1-32768 -n 1)\n")
				# f.write("r2=$(shuf -i 1-32768 -n 1)\n")

				# f.write("echo -e \"variable noiseseed equal $r1\\n\" >> parameters.dat\n")
				# f.write("echo -e \"variable smcseed equal $r2\\n\" >> parameters.dat\n")

				# f.write("echo -e \"variable replica_id equal $i\\n\" >> parameters.dat\n")
				# f.write("echo -e \"run {}\\n\" >> masterfile.lam\n".format(run_duration))

				f.write("mpirun --cpu-set $vi -display-map -n 1 -bind-to none ~/lmp -in masterfile.lam < /dev/null > out &\n")
				# f.write("wait\n")
				f.write("cd ..\n")
				f.write("if(( ${{vi}} == {} )); then\nwait\nfi\n".format(npara - 1))
				f.write("done\n")
				f.write("wait")

			os.chmod(lp_bash_fp, 0o755)

		if do_slurm:
			# generate sbatch script
			# use array size of 16 (flexible) to call different scripts running serially

			# create the replica folders on local scratch (/scratch) which is arbitrary and probably unretrievable 

			# master folder created by default


			for rep in range(nrep):
				rep_folder = os.path.join(lp_folder, f"rep{rep}")
				if not os.path.exists(rep_folder):
					os.makedirs(rep_folder)

				shutil.copy(parameter_fp, rep_folder)
				shutil.copy(master_molecule_file, rep_folder)
				shutil.copy(master_script, rep_folder)

			target_idx = c % npara

			run_str[target_idx] += "cd {}\n".format(rep_folder)
			run_str[target_idx] += "~/lmp -in masterfile.lam < /dev/null > out \n"
			print(run_str)



	if do_local:
		with open(bash_target, "w") as f:
			for _ in lp_bash_fps:
				f.write("echo " + _ + "\n")
				f.write(_ + "\n")
		os.chmod(bash_target, 0o755)

	elif do_slurm:
		slurm_job = os.path.join(slurm_folder, jobname)
		if not os.path.exists(slurm_job):
			os.makedirs(slurm_job)

		for i in range(npara):
			array_script_path = os.path.join(slurm_folder, jobname, f"run_{i}.sh")

			with open(array_script_path, "w") as f:
				print(run_str[i])
				f.write(run_str[i])

		slurm_path = os.path.join(slurm_folder, jobname + ".slurm")

		with open(slurm_path, "w") as f:
			f.write("#!/bin/bash\n")
			f.write("#SBATCH --ntasks=1\n")
			f.write("#SBATCH --cpus-per-task=1\n")
			f.write(f"#SBATCH --array=1-{npara}\n")

			f.write(f"{jobname}/run_{{SLURM_ARRAY_TASK_ID}}.sh")


		# parameters:
		# seed for langevin
		# seed for smcrun
		# lpol -> length of polymer
		# angle style, angle coeff (persistence length)
	
	# generate bash file to call all the bash files

	# wahoo

if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("-l", "--local", action = "store_true")
	ap.add_argument("-s", "--slurm", action = "store_true") # generate the SBATCH script
	ap.add_argument("job_name") # generate the SBATCH script
	args = ap.parse_args()

	main(args.local, args.slurm, args.job_name)