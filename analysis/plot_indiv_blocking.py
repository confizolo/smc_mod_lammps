import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib

def plot_blocking(input_file, plot_folder, plot_log = False, uid = "", compare_key = "", other_key = "n_stiff", n_reps = 40) :
	df = pd.read_csv(input_file)
	cmap = matplotlib.colormaps['jet']

	norm = matplotlib.colors.Normalize(vmin = 0, vmax = 1)
	nstiff = df[other_key].unique()

	if len(compare_key):
		nc = df[compare_key].nunique()

	colors = [cmap(x) for x in np.linspace(0, 1, len(nstiff))]

	stiff_map = {nstiff[x]:x for x in range(len(nstiff))}

	MAX_TIME = df["run_duration"].values[0]
	TIMESTEP = df["ratesmc"].values[0]

	df = df[df["tangent_cutoff"]!=0]
	
	if len(compare_key):
		n_compare = df[compare_key].nunique()
		fig, axs = plt.subplots(n_compare, 2)
	else:
		n_compare = 1
		fig, axs = plt.subplots(2)

	df = df.sort_values(by = "dt")
	# n_reps = len(df) / 20
	# group by force too

	for c, (c_key, t_df) in enumerate(df.groupby(compare_key)):
		for _, (p, _df) in enumerate(t_df.groupby(other_key)):
			n_stiff = p

			if not len(_df):
				print("No data points for parameter set ", p)
				continue

			# we only need aggregate statistics, we don't care about individual behaviour / more conditions yet
			#
			p_df = _df[_df["dt"] > 0]
			dts = p_df["dt"].values

			y = (-np.arange(1, len(p_df) + 1) + len(_df))/(len(_df))

			y = np.insert(y, 0, 1)
			dts = np.insert(dts, 0, TIMESTEP)
			
			yerror = 1.96 * np.sqrt((y) * (1 - y)/n_reps)

			if len(compare_key) and nc > 1:
				idx_0 = (c, 0)
				idx_1 = (c, 1)
			else:
				idx_0 = 0
				idx_1 = 1

			axs[idx_0].plot(dts, y, color = colors[stiff_map[n_stiff]], label = str(n_stiff), lw = 2, alpha = 0.6)
			axs[idx_0].fill_between(dts, y + yerror, y - yerror, color = colors[stiff_map[n_stiff]], alpha = 0.4)
			# axs[1].plot(n_stiff, len(tdf), color = colors[stiff_map[n_stiff]])

			reaching_cdf = _df.sort_values("t0")["t0"].values
			reaching_y = np.array([_/len(reaching_cdf) for _ in range(len(reaching_cdf))])

			axs[idx_1].plot(reaching_cdf, reaching_y, color = colors[stiff_map[n_stiff]], label = str(n_stiff), lw = 2, alpha = 0.8)

			# axs[idx_1].errorbar(reaching_cdf, reaching_y, yerror, color = colors[stiff_map[n_stiff]], label = str(n_stiff), lw = 2, alpha = 0.8, capsize = 2)

		axs[idx_0].set_title("{} = {:.2f}".format(compare_key, c_key))
		axs[idx_1].set_title("{} = {:.2f}".format(compare_key, c_key))

		if plot_log:
			axs[idx_0].set_xscale("log")
			axs[idx_0].set_xlabel("Simulation time (log)")
			axs[idx_1].set_xscale("log")
			axs[idx_1].set_xlabel("Simulation time (log)")
			axs[idx_0].set_xlim(TIMESTEP, MAX_TIME)
			axs[idx_1].set_xlim(TIMESTEP, MAX_TIME)
		else:
			axs[idx_0].set_xlabel("Simulation time")
			axs[idx_1].set_xlabel("Simulation time")
			axs[idx_0].set_xlim(0, MAX_TIME)
			axs[idx_1].set_xlim(0, MAX_TIME)

		axs[idx_0].set_ylabel("Blocking fraction (expt)")
		axs[idx_0].legend(title = other_key, fancybox = True)

		axs[idx_1].set_ylabel("LEF reaching stiff CDF")
		axs[idx_1].legend(title = other_key, fancybox = True)

		axs[idx_0].set_ylim(0, 1)
		axs[idx_1].set_ylim(0, 1)

	fig.set_size_inches((16, 3 * n_compare))
	plt.tight_layout()

	plot_path = os.path.join(plot_folder, os.path.basename(input_file) + uid + ".png")

	plt.savefig(plot_path, dpi = 300)
		


		

	# 	y = []
	# 	yerr = []
	# 	n = []
	# 	times = []

	# 	slice_s.append(s)

	# 	for t, tf in _df.groupby("time"):
	# 		_n = len(tf[tf["reached"] == 1])
	# 		if _n:
	# 			_y = tf.loc[tf["reached"] == 1, "blocked"].sum()/tf["reached"].sum()
	# 			y.append(_y)
	# 			yerr.append(1/np.sqrt(_n))
	# 			n.append(_n/len(tf))

	# 			if t == t_slice:
	# 				slice_y.append(_y)
	# 				slice_yerr.append(1/np.sqrt(_n))
	# 				slice_n.append(_n/len(tf))
	# 		else:
	# 			y.append(np.nan)
	# 			yerr.append(np.nan)
	# 			n.append(0)


	# 	times = _df["time"].unique()

	# 	axs[0].errorbar(times, y, yerr, color = colors[c], label = str(s), marker = "s", capsize = 3)
	# 	axs[0].set_ylim(0, 1)
	# 	axs[0].set_ylabel("Blocking fraction conditioned on LEF reaching stiff")

	# 	axs[1].plot(times, n, ".", color = colors[c], label = str(s))
	# 	axs[1].set_xlabel("Simulation time")
	# 	axs[1].set_ylabel("Fraction of LEF reaching stiff")
	# 	if xlog: 
	# 		axs[0].set_xscale("log")
	# 		axs[1].set_xscale("log")

	# plt.legend()
	# fig.set_size_inches((8.6, 9))
	# plt.tight_layout()

	# plot_path = os.path.join(plot_folder, os.path.basename(input_file) + uid + ".png")

	# plt.savefig(plot_path, dpi = 300)

	# if t_slice > 0:
	# 	fig, axs = plt.subplots(2)
	# 	fig.suptitle(f"Slice at t = {t_slice}")

	# 	axs[0].errorbar(slice_s, slice_y, slice_yerr, color = "k", capsize = 3, marker = ".")
	# 	axs[0].set_xlabel("No. of stiff beads ($l_p$ = 50)") # make this automatic later
	# 	axs[0].set_ylabel("Blocking fraction at $t$")
	# 	axs[0].set_title("Errorbars: $n_{replicas}^{-1/2}$")

	# 	axs[1].plot(slice_s, slice_n, "b.")
	# 	axs[1].set_xlabel("No. of stiff beads ($l_p$ = 50)")
	# 	axs[1].set_ylabel("Fraction of LEF reaching stiff")
	# 	# axs[1].set_ylim(0, 1)
	# 	if xlog: 
	# 		axs[0].set_xscale("log")
	# 		axs[1].set_xscale("log")

	# 	fig.set_size_inches((8.6, 9))
	# 	plt.tight_layout()

	# 	plot_path = os.path.join(plot_folder, "time_slice-" + os.path.basename(input_file) + uid + ".png")

	# 	plt.savefig(plot_path, dpi = 300)


if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("input_file")
	ap.add_argument("plot_folder")
	ap.add_argument("-c", "--compare_key", default = "")
	ap.add_argument("-o", "--other_key", default = "n_stiff")
	ap.add_argument("-l", "--log", action = "store_true")
	ap.add_argument("-id", "--uid", default = "")
	ap.add_argument("-n", "--nreps", type = int, default = 40)

	args = ap.parse_args()
	plot_blocking(args.input_file, args.plot_folder, args.log, args.uid, args.compare_key, args.other_key, args.nreps)