import math
import subprocess
import os
from pathlib import Path
import argparse
import pandas as pd

def parse(target_folders, output_file, run_time = 100000):

	DEFAULT_RATE_SMC = 100
	DEFAULT_RUN_DURATION = run_time
	DEFAULT_LPOL = 1000

	df = pd.DataFrame()
	for target_folder in target_folders:
		rep_folder_list = [str(f.parent.absolute()) for f in Path(os.path.join(target_folder)).rglob("smc_pos.txt")]

		for _, rep in enumerate(rep_folder_list):
			pos_file = os.path.join(rep, "smc_pos.txt")

			params = {
				"ratesmc": DEFAULT_RATE_SMC,
				"run_duration": DEFAULT_RUN_DURATION,
				"lpol": DEFAULT_LPOL,
				"persistence_length": 0,
				"stiff_persistence_length": 0,
				"n_stiff":0,
				"start_position":0,
				"force": 0,
				"tangent_cutoff":0,
			}

			params["replica_id"] = rep 

			_pdf = pd.read_csv(os.path.join(rep, "parameters.dat"), sep = " ", header = None, names = [ "var_name" ,"value"], usecols = [1,3])
			_pdf = _pdf.set_index("var_name")
			for p in params:
				if p in _pdf.index:
					params[p] = _pdf.loc[p]["value"]

			n_stiff = params["n_stiff"]
			lpol = params["lpol"]

			stiff_start = math.floor((lpol - n_stiff)/2)
			stiff_end = stiff_start + n_stiff

			t0_candidate = 0
			flag_blocked = 1
			flag_one_before = 0
			dt = -1

			with open(pos_file) as f:
				last_t = 0
				last_x = 0
				for line in f:
					_line = line.strip().split(" ")
					if len(_line) < 4:
						break 

					try:
						int(_line[0])
						int(_line[3])
					except:
						print("Error for ", rep)
						break

					_t = int(_line[0])
					_x = int(_line[3])

					if (_x == stiff_start) and (not flag_one_before):
						flag_one_before = 1 
						# stop considering the time
						t0_candidate = _t

					if (_x > stiff_end):
						if flag_blocked:
							flag_blocked = 0
						if (not flag_one_before): # jumped over entirely
							t0_candidate = last_t
							flag_one_before = 1

						dt = _t - t0_candidate


					if _x != last_x:
						last_x = _x
						last_t = _t
			if _t < run_time:
				continue
			params["t0"] = t0_candidate
			params["dt"] = dt

			_df = pd.DataFrame([params])
			df = pd.concat([df, _df], ignore_index = True)

	df.to_csv(output_file, index = False)


if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("-i","--input_folder", nargs = "+")
	ap.add_argument("-o","--output_file")
	ap.add_argument("-t", "--time", type = int, default = 100000)

	args = ap.parse_args()
	
	parse(args.input_folder, args.output_file, args.time)