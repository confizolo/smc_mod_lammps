from curses import mousemask
from doctest import master
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


def main(do_local, do_slurm, jobname, nrep = 96, npara = 16, run_duration = 100000, n_stiff = [4, 8, 12, 28, 40], lp_list = [20], lp_stiff = [200], add_force = False, custom_masterfile = "", regenerate_stiff = False, use_long = False, equil = "", start_shift = 20, start_index = 0, do_relax = False, force_list = [], extend_boundary = 0, tangent_cutoff_list = [0], patterns = ["1.0"], grab_cutoff = 12,ratesmc = 1000, langfric = 1.0, fake_lpst = -1, overwrite_accessible = False, **kwargs):

	def check_equil_in_folder(equil_folder, lp, lp_stiff, n_stiff, force, pattern):
		target_file = "eq_lp{:d}-{:d}_nstiff-{:d}_f-{:.2f}.dat".format(lp, lp_stiff, n_stiff, force)
		target_path = os.path.join(equil_folder, "pat{}".format(pattern), target_file)

		if os.path.exists(target_path):
			return target_path
		else:
			raise Exception("Equilibrated molecule does not exist in specified folder")


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


	parameter_set = itertools.product(lp_list, lp_stiff, n_stiff, force_list, tangent_cutoff_list, patterns)

	lpol = 1000 # length of polymer

	default_paths = {
		"local":{
			"script":'/home/zy/Documents/tap/masterfile.lam',
			"molecule_file" : '/home/zy/Documents/tap/smc-lammps/data/In_conf.Nb1000.RW4.fixed.dat',
			"script_stiff":"/home/zy/Documents/tap/masterfile_stiff.lam",
			"folder":"/home/zy/Documents/tap/smc-single-polymer/",
		},
		"slurm":{
			"script":'/home/v1zchoon/masterfile.lam',
			"molecule_file" : '/home/v1zchoon/smc-lammps/initfiles/single.dat',
			"script_stiff":"/home/v1zchoon/masterfile_stiff.lam",
			# "folder":"/home/v1zchoon/smc-single-polymer/",
			"folder":"/storage/cmstore02/groups/TAPLab/zy-smc-single-polymer/scratch",
			"slurm_output": "/storage/cmstore02/groups/TAPLab/zy-smc-single-polymer/slurm-wd",
		}
	}
	###############################
	master_folder = default_paths[path_str]["folder"]
	global_counter = 0


	bash_target = os.path.join(master_folder, jobname, 'run.sh')

	if not os.path.exists(master_folder):
		os.makedirs(master_folder)
	if do_slurm:
		slurm_folder = os.path.join(default_paths["slurm"]["slurm_output"], )
		if not os.path.exists(slurm_folder):
			os.makedirs(slurm_folder)


	lp_bash_fps = [] # collect all the different bash files to execute all at once
	run_str = ["" for _ in range(npara)]

	for p in (parameter_set):
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
			_force = p[3]
			_tangent_cutoff = p[4]
			_pattern = p[5]

			if len(equil):
				if do_relax:
					master_molecule_file = check_equil_in_folder(equil, lp, _lp_stiff, _n_stiff, 0, _pattern)
				else:
					if fake_lpst < 0:
						master_molecule_file = check_equil_in_folder(equil, lp, _lp_stiff, _n_stiff, _force, _pattern)
					else:
						master_molecule_file = check_equil_in_folder(equil, lp, fake_lpst, _n_stiff, _force, _pattern)

					# figure out the correct template, then over write it
					if overwrite_accessible:
						new_file_name = os.path.join(os.path.dirname(master_molecule_file), "inaccess_"+os.path.basename(master_molecule_file) )
						master_molecule_file = generate_stiff(_n_stiff, equil, master_molecule_file, new_file_name, force = regenerate_stiff, pattern = _pattern, overwrite_accessible=True)


			else:
				# TODO i need to fix the logic for generating molecules...
				master_molecule_file = generate_stiff(_n_stiff, os.path.join(master_folder, "stiff_molecules"), default_paths[path_str]["molecule_file"], force = regenerate_stiff, extend_boundary=extend_boundary, pattern = _pattern)

			lp_folder = os.path.join(master_folder, jobname, f"N{lpol}", f"lp{lp:02}", f"lp-stiff{_lp_stiff:02}", f"n-stiff{_n_stiff:d}", f"force_{_force:.2f}", f"tan_{_tangent_cutoff:.2f}", f"pat_{_pattern}")

			_start_position = int((lpol - _n_stiff)/2 - start_shift) # ??? e.g. for 1000 - 100, start at 430, move until 450

		# override:
		if len(custom_masterfile):
			master_script = custom_masterfile

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

			f.write("variable lang_fric equal {:.1f}\n".format(langfric))

			f.write("variable stiff_persistence_length equal {:.2f}\n".format(_lp_stiff))
			f.write("variable start_position equal {:d}\n".format(_start_position))
			f.write("variable n_stiff equal {:d}\n".format(_n_stiff))
			f.write("variable max_jump equal {:d}\n".format(int(lpol/2)))
			f.write("variable init_force equal {:.2f}\n".format(_force))
			f.write("variable tangent_cutoff equal {:.2f}\n".format(_tangent_cutoff))
			f.write("variable pattern equal {}\n".format(_pattern))
			f.write("variable grab_cutoff equal {:.2f}\n".format(grab_cutoff))

			if add_force:
				f.write("group end1 id 1\n")
				f.write("group end2 id 1000\n")
				f.write("fix tension1 end1 addforce -1.0 0.0 0.0\n")
				f.write("fix tension2 end2 addforce 1.0 0.0 0.0\n")

		if do_local:
			with open(lp_bash_fp, 'w') as f:
				f.write(f"nparajobs={npara}\n")
				f.write("cd {}\n".format(lp_folder))
				f.write("for i in {{{}..{}}}; do\n".format(start_index, start_index + nrep - 1))
				f.write("mkdir -p rep$i\n")
				f.write("cd rep$i\n")
				f.write("vi=$(( ${i}%${nparajobs} ))\n")

				f.write("cp -R -p -u ../parameters.dat .\n".format(parameter_fp)) # parameter file
				f.write("cp -R -p -u {} masterfile.lam\n".format(master_script)) # lams script
				f.write("cp -R -p -u {} molecule.dat\n".format(master_molecule_file))
				f.write("r1=$(shuf -i 1-32768 -n 1)\n")
				f.write("r2=$(shuf -i 1-32768 -n 1)\n")

				f.write("echo -e \"variable noiseseed equal $r1\\n\" >> parameters.dat\n")
				f.write("echo -e \"variable smcseed equal $r2\\n\" >> parameters.dat\n")

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
				rep_folder = os.path.join(lp_folder, "rep{}".format(start_index + rep))
				if not os.path.exists(rep_folder):
					os.makedirs(rep_folder)


				shutil.copy(parameter_fp, rep_folder)

				with open(os.path.join(rep_folder, "parameters.dat"), "a") as f:
					f.write("variable noiseseed equal {:d}\n".format((random.randint(1, 32768))))
					f.write("variable smcseed equal {:d}\n".format((random.randint(1, 32768))))

				with open(os.path.join(rep_folder, "start_pos.dat"), "w") as f:
					f.write("{} {}".format(_start_position-2, _start_position))

				shutil.copy(master_molecule_file, os.path.join(rep_folder, "molecule.dat"))
				shutil.copy(master_script, os.path.join(rep_folder, "masterfile.lam"))

				target_idx = global_counter % npara
				global_counter += 1

				run_str[target_idx] += "cd {}\n".format(rep_folder)
				run_str[target_idx] += "~/lmp -in masterfile.lam < /dev/null > out \n"

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
				f.write(run_str[i])
			os.chmod(array_script_path, 0o755)

		slurm_path = os.path.join(slurm_folder, jobname + ".slurm")

		with open(slurm_path, "w") as f:
			f.write("#!/bin/bash\n")
			f.write("#SBATCH --ntasks=1\n")
			f.write("#SBATCH --cpus-per-task=1\n")
			f.write(f"#SBATCH --array=0-{npara-1}\n")
			f.write(f"#SBATCH --job-name={jobname}\n")
			if use_long:
				f.write(f"#SBATCH --partition=long\n")
			else:
				f.write(f"#SBATCH --partition=short\n")

			f.write("{}/{}/run_${{SLURM_ARRAY_TASK_ID}}.sh".format(
				slurm_folder,
				jobname,
			))

		print(os.path.abspath(os.path.join(master_folder, jobname)))

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
	ap.add_argument("-r", "--nrep", type = int, default = 96)
	ap.add_argument("-p", "--npara", type = int, default = 8)
	ap.add_argument("-t", "--run_time", type = int, default = 100000)
	ap.add_argument("-lpst", "--lp_stiff", type = int, nargs = "+", default = [200], help = "lp of stiff section bead")
	ap.add_argument("-lp", "--lp", type = int, nargs = "+", default = [20], help = "lp of default bead")
	ap.add_argument("-nst", "--n_stiff", type = int, nargs = "+", default = [4, 8, 12, 28, 40])

	ap.add_argument("-ff", "--forces", type = float, nargs="+", default = [0])

	ap.add_argument("-tan", "--tangent_cutoff", type = float, nargs="+", default = [1.0])

	ap.add_argument('-f', '--add_force', action = "store_true")
	ap.add_argument('-m', '--masterfile', default = "")
	ap.add_argument("-rgs", "--regenerate_stiff", action = "store_true")
	ap.add_argument("-pl", "--long", action = "store_true", help = "flag to use partition `long`")

	# for now, let this look from a folder of already equilibrated molecules

	ap.add_argument("-eq", "--equilibrate", type = str, help = "path to folder", default = "")
	ap.add_argument("-ss", "--start_shift", type = int, default = 20)
	ap.add_argument("-si", "--start_index", type = int, default = 0)

	ap.add_argument("-rx", "--do_relax", action = "store_true")
	ap.add_argument("-ex", "--extend_boundary", type = float, default = 0.0)

	ap.add_argument("-pat", "--patterns", nargs = "+", type = str, default = ["1.0"])

	ap.add_argument("-cut", "--cutoff", default = 12., type = float)
	ap.add_argument("-dt","--ratesmc", type = int, default = 1000)
	ap.add_argument("-lf", "--langfric", default = 1.0, type = float)
	ap.add_argument("-flpst", "--fake_lpst", default = -1, type = int)
	# ap.add_argument("-tem", "--custom_molecule_template", default = "")
	ap.add_argument("-oa", "--overwrite_accessible", action = "store_true")

	ap.add_argument("job_name") # generate the SBATCH script

	args = ap.parse_args()

	main(args.local, args.slurm, args.job_name,
		nrep = args.nrep,
		npara = args.npara,
		run_duration = args.run_time,
		n_stiff = args.n_stiff,
		lp_stiff = args.lp_stiff,
		lp_list= args.lp,
		add_force = args.add_force,
		custom_masterfile = args.masterfile,
		regenerate_stiff = args.regenerate_stiff,
		use_long = args.long,
		equil = args.equilibrate,
		start_shift = args.start_shift,
		start_index = args.start_index,
		do_relax = args.do_relax,
		force_list = args.forces,
		extend_boundary=args.extend_boundary,
		tangent_cutoff_list=args.tangent_cutoff,
		patterns = args.patterns,
		grab_cutoff = args.cutoff,
		ratesmc = args.ratesmc,
		langfric = args.langfric,
		fake_lpst = args.fake_lpst,
		overwrite_accesible = args.overwrite_accessible,
	)
