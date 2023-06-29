import itertools
import os
from stiff import generate_stiff

def generate_initial_molecule():
	# probably not needed? there's a cpp script that does this, could just wrap over it lol
	pass

def generate_master_lams_file():
	# to vary bond coefficients, for future proofing
	pass

def main():
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
	master_folder = '/home/zy/Documents/tap/smc-single-polymer/' 
	jobname = 'stiff-14-fix-edgecase'

	start_delta = 20 # hard coded in the fix... 

	lp_stiff = [50]
	n_stiff = [1, 2, 3, 7, 10, 14]

	lp_list = [5] # list of persistence lengths to run through

	parameter_set = itertools.product(lp_list, lp_stiff, n_stiff)

	lpol = 1000 # length of polymer

	ratesmc = 100 # smc movement attempt 
	run_duration = 100000 # total number of steps

	nrep = 96 # number of replicas to do
	npara = 8 # no. of parallel jobs
	start_index = 0
	#
	#
	#
	###############################

	bash_target = os.path.join(master_folder, jobname, 'run.sh')

	if not os.path.exists(master_folder):
		os.makedirs(master_folder)
	

	lp_bash_fps = [] # collect all the different bash files to execute all at once

	for p in parameter_set:
		lp = p[0]

		if p[1] == None or p[2] == None:
			master_script = '/home/zy/Documents/tap/masterfile.lam'
			master_molecule_file = '/home/zy/Documents/tap/smc-lammps/data/In_conf.Nb1000.RW4.fixed.dat'

			lp_folder = os.path.join(master_folder, jobname, f"N{lpol}", f"lp{lp:02}")
			_start_position = 0

		else:
			master_script = '/home/zy/Documents/tap/masterfile_stiff.lam'

			_lp_stiff = p[1]
			_n_stiff = p[2]
			master_molecule_file = generate_stiff(_n_stiff)

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
			# f.write(f"angle_coeff 1 {lp:d}") # assumes that persistence length lp is an integer

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

	with open(bash_target, "w") as f:
		for _ in lp_bash_fps:
			f.write("echo " + _ + "\n")
			f.write(_ + "\n")
	os.chmod(bash_target, 0o755)


		# parameters:
		# seed for langevin
		# seed for smcrun
		# lpol -> length of polymer
		# angle style, angle coeff (persistence length)
	
	# generate bash file to call all the bash files

	# wahoo

if __name__ == "__main__":
	main()