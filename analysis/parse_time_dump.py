import argparse
import os
from pathlib import Path
from stiff import parse_pattern
import numpy as np 
import pandas as pd

"""
parse the movement file
"""

def parse_lammps(input_file):
	t = 0
	n_atoms = 1000

	all_data = {}

	with open(input_file, "r") as f:
		for line in f:
			if "ITEM: TIMESTEP" in line.strip():
				line = next(f) # weird hack but ok
				t = int(line.strip())
			if "ITEM: NUMBER OF ATOMS" in line.strip():
				line = next(f)
				n_atoms = int(line.strip())
			if "ITEM: ATOMS" in line.strip():
				# figure out the format
				atom_data = {}
				fmt = line.strip().split(" ")[2:]
				for i in range(n_atoms):
					line = next(f).strip().split(" ")
					_data = {k:v for k, v, in zip(fmt[1:], line[1:], )}
					# the data is stored as strings, so they'll need to be manually parsed
					for k in _data:
						if k == "type":
							_data[k] = int(_data[k])
						elif k in ["x", "y", "z"]:
							_data[k] = float(_data[k])

					atom_data[int(line[0])] = _data
				
				all_data[t] = atom_data

	return all_data


def parse_positions(input_file, timestep = 1000):
	# sample the locations every % timestep

	all_data = {}

	with open(input_file, "r") as f:
		for line in f:
			x = [int(_) for _ in line.strip().split(" ")]
			if x[0] % timestep == 0:
				all_data[x[0]] = x[-1]

	return all_data
	

def worker(pos_file, trj_file, params):
	# location_data = parse_lammps("../../workbench/timedump/pattern4.3_rate200_cutoff12_force_0.lammpstrj", "")
	# smc_data = parse_positions("../../workbench/timedump/smc_pos.txt")

	location_data = parse_lammps(trj_file)
	smc_data = parse_positions(pos_file)
	# output_file = "../../workbench/test_output.csv"


	CUTOFF = params["grab_cutoff"] # can infer from parameter set
	pattern = "{:.1f}".format(params["pattern"])
	n_stiff = params["n_stiff"]

	stiff_start, stiff_end, build_str = parse_pattern(n_stiff, pattern)

	# print(len(location_data.keys()))
	# print(len(smc_data.keys()))

	# we are tracking: the total no. of accessible beads in FRONT (radius wise)
	# and out of those, how many beyond the array, before the array, and within the array

	# location will have t = 0, but not position (lol)

	df = pd.DataFrame(columns = ["time", "n_before", "n_middle", "n_after", 'n_total'])

	for k in smc_data:
		smc_pos = smc_data[k]
		if smc_pos > stiff_end:
			# no longer interested
			break
		# build a numpy array
		# compute the dist squared

		smc_loc = np.array([location_data[k][smc_pos][_] for _ in ["x", "y", "z"]])

		_data = [k, 0,0,0,0]
		for atom_no, atom in location_data[k].items():
			if atom_no <= smc_pos: continue # only beads in front
			if atom["type"] != 1: continue # only accessible beads

			atom_loc = np.array([location_data[k][atom_no][_] for _ in ["x", "y", "z"]])

			dist_sq = np.sum((atom_loc - smc_loc)**2)
			# print(dist_sq)
			if dist_sq < CUTOFF * CUTOFF:
				if atom_no < stiff_start: 
					_data[1] += 1
				elif atom_no > stiff_end:
					_data[3] += 1
				else:
					_data[2] += 1

				_data[-1] += 1
		df.loc[len(df.index)] = _data

	df = df.replace(0, np.nan,)

	df["middle_pct"] = (df["n_middle"]) / (df["n_middle"] + df["n_after"])

	return df

def control(target_folders, output_file):
	df_list = []
	rep_list = []
	for target_folder in target_folders:
		rep_folder_list = [str(f.parent.absolute()) for f in Path(os.path.join(target_folder)).rglob("smc_pos.txt")]

		for _, rep in enumerate(rep_folder_list):
			print(_, rep)
			pos_file = os.path.join(rep, "smc_pos.txt")

			params = {
				"ratesmc": 200,
				"run_duration": 1000000,
				"lpol": 1000,
				"persistence_length": 0,
				"stiff_persistence_length": 0,
				"n_stiff":0,
				"start_position":0,
				"init_force": 0,
				"tangent_cutoff":0,
				"pattern":"1.0",
				"grab_cutoff":0,
				"lang_fric":1,
			}

			params["replica_id"] = rep 
			rep_list.append(rep)

			_pdf = pd.read_csv(os.path.join(rep, "parameters.dat"), sep = " ", header = None, names = [ "var_name" ,"value"], usecols = [1,3])
			_pdf = _pdf.set_index("var_name")
			for p in params:
				if p in _pdf.index:
					params[p] = _pdf.loc[p]["value"]

			trj_file = [str(_) for _ in Path(rep).glob("*.lammpstrj")][0]

			output_df = worker(pos_file, trj_file, params)

			df_list.append(output_df)

	df = pd.concat(df_list, axis = 1, keys = rep_list)
	# save the df lol
	# this is going to be SLOW
	# unless i add multiprocessing
	df.to_csv(output_file, index = False)
			
if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("-i", "--target_folders", nargs = "+", default = [])
	ap.add_argument("-o", "--output_file", default = "")

	args = ap.parse_args()

	control(args.target_folders, args.output_file)


