import math
import subprocess
import os
from pathlib import Path
import argparse
import pandas as pd

def parse_pattern(sl, x):
	bunch, gap = x.split(".")
	bunch, gap = int(bunch), int(gap)

	_sl = sl + ((sl/bunch)-1) * gap
	_start = int((1000 - _sl)/2 + 1) # if sl = 0, this is 501
	_end = int(_start + _sl) # if sl = 0, this is 501

	_build_str = ""

	while _sl > 0:
		_build_str += "o" * bunch
		_sl -= bunch
		_build_str += "x" * gap

	return _start, _end, _build_str

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
				"init_force": 0,
				"tangent_cutoff":0,
				"pattern":"1.0",
				"grab_cutoff":0,
			}

			params["replica_id"] = rep 

			_pdf = pd.read_csv(os.path.join(rep, "parameters.dat"), sep = " ", header = None, names = [ "var_name" ,"value"], usecols = [1,3])
			_pdf = _pdf.set_index("var_name")
			for p in params:
				if p in _pdf.index:
					params[p] = _pdf.loc[p]["value"]

			stiff_start, stiff_end, _build_str = parse_pattern(params["n_stiff"], "{:.1f}".format(params["pattern"]))
			print(stiff_start, stiff_end)
			print(pos_file)

			t0_candidate = 0
			dt = -1

			reached_start = False
			reached_end = False

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

					if (_x >= stiff_start - 1) and not reached_start:
						if _x < stiff_end:
							reached_start = True
							t0_candidate = _t
						else:
							t0_candidate = _t - params["ratesmc"]
							reached_end = True
							
					if _x > stiff_end and not reached_end:
						reached_end = True
						dt = _t - t0_candidate

			if _t < run_time:
				continue
			params["t0"] = t0_candidate
			params["dt"] = dt
			print(t0_candidate, dt)
			print()

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