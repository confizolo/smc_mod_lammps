import numpy as np
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

def main():

	input_file = "~/Documents/tap/smc-single-polymer/stiff_test-redo.csv"
	lpol = 1000

	df = pd.read_csv(input_file)

	stiffs = df["n_stiff"].unique()

	output_df = pd.DataFrame(columns = ["n_stiff", "blocked"])

	timestep = 10000

	for c, x in filter_generator(df, stiffs, get_lp_id = True):

		_stiff = stiffs[c]
		stiff_start = math.floor((lpol - _stiff)/2)
		stiff_end = stiff_start + _stiff

		max_time = df.loc[x, "time"].max() 

		_df = {"n_stiff": [], "blocked": [], "reached": [], "time": []}
		for t in range(timestep, max_time + 1, timestep):

			tf = (df["time"] == t)

			if df.loc[x & tf, "x2"].values[0] > stiff_end:
				_df["blocked"].append(0)
			else:
				_df["blocked"].append(1)

			if df.loc[x & tf, "x2"].values[0] < stiff_start:
				_df["reached"].append(0)
			else:
				_df["reached"].append(1)

			_df["n_stiff"].append(_stiff)
			_df["time"].append(t)

		_df = pd.DataFrame(_df)
		output_df = pd.concat([output_df, _df], ignore_index = True)

	output_df.to_csv("~/Documents/tap/smc-single-polymer/stiff_fix_analysis.csv", index = False)


def plot_blocking(input_file):
	input_file = "~/Documents/tap/smc-single-polymer/stiff_fix_analysis.csv"
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

	for c, (s, _df) in enumerate(df.groupby("n_stiff")):
		y = []
		yerr = []
		n = []
		for t, tf in _df.groupby("time"):
			_n = len(tf[tf["reached"] == 1])
			if _n:
				y.append(tf.loc[tf["reached"] == 1, "blocked"].sum()/tf["reached"].sum())
				yerr.append(1/np.sqrt(_n))
				n.append(_n/len(tf))
			else:
				y.append(None)
				yerr.append(None)
				n.append(0)

		times = _df["time"].unique()
		axs[0].errorbar(times, y, yerr, color = colors[c], label = str(s), marker = "s", capsize = 3)
		axs[0].set_ylim(0, 1)
		axs[0].set_ylabel("Blocking fraction conditioned on LEF reaching stiff")
		axs[1].set_ylabel("Fraction of LEF reaching stiff")
		axs[1].plot(times, n, ".", color = colors[c], label = str(s))
		axs[1].set_xlabel("Simulation time")
	plt.legend()
	fig.set_size_inches((8.6, 9))
	plt.tight_layout()
	plt.savefig("plots/stiff_fix_blocking_over_time.png", dpi = 300)


if __name__ == "__main__":
	# main()
	plot_blocking("")