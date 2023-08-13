import shutil
import argparse
import pandas as pd
import os
from pathlib import Path


def main(input_folder, output_folder, source_file_name = "equil.dat", ):
	# only need one input folder since we assume that we equilibrate all at once

	# search for equil.dat outputs
	# assumes that only ONE replica per parameter set

	# infer n_stiff and applied force
	# infer lp and lpstiff
	# /storage/scratch/v1zchoon/smc-single-polymer/equil_smallnstiff/N1000/lp20/lp-stiff200/n-stiff10/force_0.00/rep0/

	rep_folder_list = [str(f.parent.absolute()) for f in Path(os.path.join(input_folder)).rglob(source_file_name)]

	for rep in rep_folder_list:

		params = {
			"persistence_length": 0,
			"stiff_persistence_length": 0,
			"n_stiff":0,
			'init_force': 0,
			"pattern":"1.0",
		}

		# rep_id = (rep.split("/")[-1][3:])
		params["replica_id"] = rep

		_pdf = pd.read_csv(os.path.join(rep, "parameters.dat"), sep = " ", header = None, names = [ "var_name" ,"value"], usecols = [1,3])
		_pdf = _pdf.set_index("var_name")
		for p in params:
			if p in _pdf.index:
				params[p] = _pdf.loc[p]["value"]

		output_filename = "eq_lp{:d}-{:d}_nstiff-{:d}_f-{:.2f}.dat".format(
			int(params["persistence_length"]),
			int(params["stiff_persistence_length"]),
			int(params["n_stiff"]),
			params["init_force"]
		)
		output_path = os.path.join(output_folder, "pat{}".format(params["pattern"]), output_filename)

		print(output_path)
		shutil.copy(os.path.join(rep, source_file_name), output_path)


if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("input_folder")
	ap.add_argument("output_folder")
	ap.add_argument("-n", "--source_file_name", type = str, default = "equil.dat")

	args = ap.parse_args()
	main(args.input_folder, args.output_folder, source_file_name = args.source_file_name)
