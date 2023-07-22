import argparse 
import subprocess
from pathlib import Path
import os

def parse(target_folders, log_file):

	DEFAULT_RATE_SMC = 100
	DEFAULT_RUN_DURATION = 100000
	DEFAULT_LPOL = 1000
	print_rate = 100
	expected_print_lines = DEFAULT_RUN_DURATION/print_rate

	with open(log_file, "w") as f:

		for target_folder in target_folders:
			# folder_id = target_folder.split("/")[-3] # manual, will break
			# folder_list = [(f.path, f.name) for f in os.scandir(target_folder) if f.is_dir()]

			rep_folder_list = [str(_.parent.absolute()) for _ in Path(os.path.join(target_folder)).rglob("smc_pos.txt")]

			for _, rep in enumerate(rep_folder_list):
				pos_file = os.path.join(rep, "smc_pos.txt")
				result = subprocess.check_output(['wc', '-l', pos_file, ])
				obs_lines = int(result.decode("utf-8").split(" ")[0])

				if obs_lines < expected_print_lines:
					print(rep)
					f.write(rep + "\n")

if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("-i","--input_folder", nargs = "+")
	ap.add_argument("-o","--output_file")

	args = ap.parse_args()
	
	parse(args.input_folder, args.output_file)




