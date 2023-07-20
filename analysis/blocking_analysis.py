import os
import numpy as np
import argparse
import matplotlib
import math
import matplotlib.pyplot as plt
import pandas as pd

def filter_generator(df, lps, get_lp_id = False):
	for c, lp in enumerate(lps):
		lp_filter = (df["n_stiff"] == lp)
		# get reps
		reps = df.loc[lp_filter, "replica_id"].unique()

		for rep in reps:
			if get_lp_id:
				yield c, (df["replica_id"] == rep) & lp_filter
			else:
				yield (df["replica_id"] == rep) & lp_filter

def main(input_file):
	"""
	blocking fraction starts from 100%

	case: immediate jump over, t = 0 at the last timepoint before jumping 
	case: stepwise reach stiff, t = 0 at reaching point before stiff
	
	want to shift the start time point
	"""

	lpol = 1000
	timestep = 100

	df = pd.read_csv(input_file)

	stiffs = df["n_stiff"].unique()

	output_file = os.path.join(os.path.dirname(input_file), "processed_" + os.path.basename(input_file))

	output_df = pd.DataFrame(columns = ["n_stiff", "blocked"])

	global_max_time = df["time"].max()

	for c, x in filter_generator(df, stiffs, get_lp_id = True):
		print(c)
		max_time = df.loc[x, "time"].max() 
		if max_time < global_max_time: # handle crashed cases due to fene bond...
			continue

		_stiff = stiffs[c]
		stiff_start = math.floor((lpol - _stiff)/2)
		stiff_end = stiff_start + _stiff

		_df = {"n_stiff": [], "blocked": [], "reached": [], "time": []}

		flag_blocked = 1
		flag_one_before = 0
		time_zero_candidate = 0

		for t in range(timestep, max_time + 1, timestep):
			tf = (df["time"] == t)

			_x = df.loc[x & tf, "x2"].values[0] 

			if (_x == stiff_start - 1) and (not flag_one_before):
				flag_one_before = 1 # stop considering
				time_zero_candidate = t

			if (_x > stiff_end):
				if flag_blocked:
					flag_blocked = 0
				if (not flag_one_before): # jumped over entirely
					prev_pos = df.loc[(x & df["time"] == last_time), "x2"].values[0]
					time_zero_candidate = df.loc[(x & df["x2"] == prev_pos), "time"].values[0]
					flag_one_before = 1 # stop considering

			_df["blocked"].append(flag_blocked)
			_df["n_stiff"].append(_stiff)
			_df["time"].append(t)
			_df["t_0"] = time_zero_candidate

			last_time = t

		_df = pd.DataFrame(_df)

		_df = _df[_df["time"] >= time_zero_candidate]
		output_df = pd.concat([output_df, _df], ignore_index = True)

	output_df.to_csv(output_file, index = False)

	return output_file


def plot_blocking(input_file, plot_folder, t_slice = 200000, xlog = False, uid = "") :
	df = pd.read_csv(input_file)

	fig, axs = plt.subplots(2)
	# first plot is blocking fraction conditioned on reaching
	# add error bars lol
	# second plot is N data points in bin

	# different colours

	cmap = matplotlib.colormaps['jet']

	norm = matplotlib.colors.Normalize(vmin = 0, vmax = 1)
	nstiff = df["n_stiff"].nunique()
	colors = [cmap(x) for x in np.linspace(0, 1, nstiff )]

	slice_y = []
	slice_yerr = []
	slice_n = []
	slice_s = []

	for c, (s, _df) in enumerate(df.groupby("n_stiff")):
		y = []
		yerr = []
		n = []
		times = []

		slice_s.append(s)

		for t, tf in _df.groupby("time"):
			_n = len(tf[tf["reached"] == 1])
			if _n:
				_y = tf.loc[tf["reached"] == 1, "blocked"].sum()/tf["reached"].sum()
				y.append(_y)
				yerr.append(1/np.sqrt(_n))
				n.append(_n/len(tf))

				if t == t_slice:
					slice_y.append(_y)
					slice_yerr.append(1/np.sqrt(_n))
					slice_n.append(_n/len(tf))
			else:
				y.append(np.nan)
				yerr.append(np.nan)
				n.append(0)


		times = _df["time"].unique()

		axs[0].errorbar(times, y, yerr, color = colors[c], label = str(s), marker = "s", capsize = 3)
		axs[0].set_ylim(0, 1)
		axs[0].set_ylabel("Blocking fraction conditioned on LEF reaching stiff")

		axs[1].plot(times, n, ".", color = colors[c], label = str(s))
		axs[1].set_xlabel("Simulation time")
		axs[1].set_ylabel("Fraction of LEF reaching stiff")
		if xlog: 
			axs[0].set_xscale("log")
			axs[1].set_xscale("log")

	plt.legend()
	fig.set_size_inches((8.6, 9))
	plt.tight_layout()

	plot_path = os.path.join(plot_folder, os.path.basename(input_file) + uid + ".png")

	plt.savefig(plot_path, dpi = 300)

	if t_slice > 0:
		fig, axs = plt.subplots(2)
		fig.suptitle(f"Slice at t = {t_slice}")

		axs[0].errorbar(slice_s, slice_y, slice_yerr, color = "k", capsize = 3, marker = ".")
		axs[0].set_xlabel("No. of stiff beads ($l_p$ = 50)") # make this automatic later
		axs[0].set_ylabel("Blocking fraction at $t$")
		axs[0].set_title("Errorbars: $n_{replicas}^{-1/2}$")

		axs[1].plot(slice_s, slice_n, "b.")
		axs[1].set_xlabel("No. of stiff beads ($l_p$ = 50)")
		axs[1].set_ylabel("Fraction of LEF reaching stiff")
		# axs[1].set_ylim(0, 1)
		if xlog: 
			axs[0].set_xscale("log")
			axs[1].set_xscale("log")

		fig.set_size_inches((8.6, 9))
		plt.tight_layout()

		plot_path = os.path.join(plot_folder, "time_slice-" + os.path.basename(input_file) + uid + ".png")

		plt.savefig(plot_path, dpi = 300)


if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("input_file")
	ap.add_argument("plot_folder")
	ap.add_argument("-n", "--no_compute", action = "store_false", default = True)
	ap.add_argument("-np", "--no_plot", action = "store_true")
	ap.add_argument("-t", "--t_slice", type = int, default = -1)
	ap.add_argument("-l", "--xlog", action = "store_true")
	ap.add_argument("-id", "--uid", default = "")

	args = ap.parse_args()

	if args.no_compute:
		analysis_file = main(args.input_file)
	else:
		analysis_file = os.path.join(os.path.dirname(args.input_file), "processed_" + os.path.basename(args.input_file))

	if not args.no_plot:
		plot_blocking(analysis_file, args.plot_folder, args.t_slice, args.xlog, args.uid)