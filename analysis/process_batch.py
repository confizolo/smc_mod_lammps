from pathlib import Path
from glob import glob
from scipy import signal
import numpy as np
import matplotlib
import os
from glob import glob
import pandas as pd
import matplotlib.pyplot as plt


def parse():

	DEFAULT_RATE_SMC = 100
	DEFAULT_RUN_DURATION = 100000
	DEFAULT_LPOL = 1000

	# what about the length of the polymer? assume that it is held constant?
	# future proofing for e.g. varying length and LP? 
	# in that case, make the hierarchy the same i.e. job_name/lpX_NY/stuff_inside

	target_folders = [
		"/home/zy/Documents/tap/smc-single-polymer/stiff-14-fix-edgecase/N1000/", 
	]


	df = pd.DataFrame(columns = ['persistence_length', 'replica_id', 'lpol', 'ratesmc', 'run_duration', 'time', 'smc_id', 'x1', 'x2', 'stiff_persistence_length', 'n_stiff', 'fixed_start_pos'])

	for target_folder in target_folders:
		# folder_id = target_folder.split("/")[-3] # manual, will break
		# folder_list = [(f.path, f.name) for f in os.scandir(target_folder) if f.is_dir()]

		rep_folder_list = [str(f.parent.absolute()) for f in Path(os.path.join(target_folder)).rglob("smc_pos.txt")]

		for _, rep in enumerate(rep_folder_list):

			params = {
				"ratesmc": DEFAULT_RATE_SMC,
				"run_duration": DEFAULT_RUN_DURATION,
				"lpol": DEFAULT_LPOL,
				"persistence_length": 0,
				"stiff_persistence_length": 0,
				"n_stiff":0,
				"start_position":0,
			}

			# rep_id = (rep.split("/")[-1][3:])
			params["replica_id"] = rep 

			_pdf = pd.read_csv(os.path.join(rep, "parameters.dat"), sep = " ", header = None, names = [ "var_name" ,"value"], usecols = [1,3])
			_pdf = _pdf.set_index("var_name")
			for p in params:
				if p in _pdf.index:
					params[p] = _pdf.loc[p]["value"]

			pos_file = os.path.join(rep, "smc_pos.txt")
			_df = pd.read_csv(pos_file, sep = " ", header = None, names = ['time', 'smc_id', 'x1', 'x2'])

			for p in params:
				_df[p] = params[p]
			df = pd.concat([df, _df], ignore_index = True)

	df.to_csv("~/Documents/tap/smc-single-polymer/stiff_test-redo.csv", index = False)


def shift_start(df):
	df = df.sort_values("time")
	df = df.reset_index()

	df["x1"] -= df.iloc[0]["x1"]
	df["x2"] -= df.iloc[0]["x2"]

	return df

def compute_velocity(df):

	# will need to compute a "distanced moved" once we implement async movement of the beads

	for c, row in df.iterrows():
		df.at[c, "hash"] = "{}_{}".format(int(row.x1), int(row.x2))

	df = df.drop_duplicates(subset = "hash", ignore_index = True)

	df.loc[:, "time_delta"] = df.loc[:, "time"] - df.loc[:, "time"].shift(1)

	df.loc[:, "v_1"] = 1 / df.loc[:, "time_delta"] 

	df["smoothed_v_1"] = df["v_1"].rolling(4).mean()

	for c, row in df.iterrows():
		if c == 0: pass
		else:
			df.at[c, "v_2"] = 1/(row.time_delta) + 0.5 * ((1/(row.time_delta)) - 1/(prev_dt))
		prev_dt = row.time_delta

	return df

def average_plotter():
	df = pd.read_csv("~/Documents/tap/smc-single-polymer/test_scan_lp/merged.csv")
	nlp = df["persistence_length"].nunique()

	cmap = matplotlib.colormaps['jet']
	norm = matplotlib.colors.Normalize(vmin = 0, vmax = 1)
	colors = [cmap(x) for x in np.linspace(0, 1, nlp)][::-1]
	# print(colors)

	fig, axs = plt.subplots(2, nlp)
	sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

	# are these expensive operations?? i have no idea

	reps = df["replica_id"].unique()
	lps = df["persistence_length"].unique()
	times = df["time"].unique()

	df = df.sort_values(["time","persistence_length", "replica_id"])

	df = df.drop(df[(df["persistence_length"] == 5) & (df["replica_id"] == 21 )].index)


	# directly modify the original dataframe so you
	# can slice it in different ways


	adf = pd.DataFrame(columns = ["time", "lp", "avg_pos", ])

	for lp in lps:
		lp_filter = (df["persistence_length"] == lp)
		for rep in reps:
			# manual lol
			if lp == 5 and rep == 21: continue
			rep_filter = (df["replica_id"] == rep)

			start_time = (df.loc[lp_filter & rep_filter]).reset_index().loc[0, "time"]

			x1_0, x2_0 =  (df.loc[lp_filter & rep_filter & (df["time"] == start_time), ["x1", "x2"]]).values[0]

			df.loc[lp_filter & rep_filter, "x1"] -= x1_0
			df.loc[lp_filter & rep_filter, "x2"] -= x2_0

	
	# for each lp, compute average position across all replicas for each point in time
		for t in times:
			time_filter = (df["time"] == t)
			mean_x2 = (df.loc[lp_filter & (time_filter), "x2"].mean())


			new_df = pd.DataFrame([{"time": t, "lp": lp, "avg_pos": mean_x2}])

			adf = pd.concat([adf, new_df], axis = 0, ignore_index = True)

		# for each lp, compute the velocity using the first order difference method
	
	adf = adf.sort_values(["lp", "time"])
	# continual motion, no need for hash

	fig, axs = plt.subplots(5)
	lps.sort()

	# try adding a second order velocity computation
	# and try using a rolling window to smooth

	for c, lp in enumerate(lps):
		adf_lp_filter = (adf["lp"] == lp)

		adf.loc[adf_lp_filter, "pos_delta"] = adf.loc[adf_lp_filter, "avg_pos"] - adf.loc[adf_lp_filter, "avg_pos"].shift(1)

		adf.loc[adf_lp_filter, "v_1"] = adf.loc[adf_lp_filter, "pos_delta"] / 1000 #... this is manual

		adf.loc[adf_lp_filter, "v_2"] = 1.5 * adf.loc[adf_lp_filter, "pos_delta"] - 0.5 * adf.loc[adf_lp_filter, "pos_delta"].shift(1)

		window_size = 4

		adf.loc[adf_lp_filter, "smoothed_v_1"] = adf.loc[adf_lp_filter, "v_1"].rolling(window_size).mean()

		adf.loc[adf_lp_filter, "smoothed_v_2"] = adf.loc[adf_lp_filter, "v_2"].rolling(window_size).mean()

		axs[0].plot(adf.loc[adf_lp_filter, "time"], adf.loc[adf_lp_filter, "avg_pos"], color = colors[c], label = "lp = {:.2f}".format(lp), alpha = 1)

		axs[1].plot(adf.loc[adf_lp_filter, "time"], adf.loc[adf_lp_filter,"v_1"], color = colors[c], label = "lp = {:.2f} 1st".format(lp), alpha = 1)
		# axs[1].plot(adf.loc[adf_lp_filter, "time"], adf.loc[adf_lp_filter,"v_2"], color = colors[c], label = "lp = {:.2f} 2nd".format(lp), alpha = 0.5, linestyle = "--")

		axs[2].plot(adf.loc[adf_lp_filter, "avg_pos"], adf.loc[adf_lp_filter, "v_1"], color = colors[c], label = "lp = {:.2f} 1st".format(lp), alpha = 1, linestyle = "-")
		# axs[2].plot(adf.loc[adf_lp_filter, "avg_pos"], adf.loc[adf_lp_filter, "v_2"], color = colors[c], label = "lp = {:.2f} 2nd".format(lp), alpha = 0.5, linestyle = "--")


		axs[3].plot(adf.loc[adf_lp_filter, "time"], adf.loc[adf_lp_filter, "smoothed_v_1"], color = colors[c], linestyle = "-", label = "lp = {:.2f} 1st".format(lp))

		# axs[3].plot(adf.loc[adf_lp_filter, "time"], adf.loc[adf_lp_filter, "smoothed_v_2"], color = colors[c], linestyle = "--", label = "lp = {:.2f} 2nd".format(lp))

		axs[4].plot(adf.loc[adf_lp_filter, "avg_pos"], adf.loc[adf_lp_filter, "smoothed_v_1"], color = colors[c], label = "lp = {:.2f} 1st".format(lp), alpha = 1, linestyle = "-")
		# axs[4].plot(adf.loc[adf_lp_filter, "avg_pos"], adf.loc[adf_lp_filter, "smoothed_v_2"], color = colors[c], label = "lp = {:.2f} 2nd".format(lp), alpha = 0.5, linestyle = "--")


		# df.loc[:, "v_1"] = 1 / df.loc[:, "time_delta"] 


	axs[3].set_title("Smoothing window = {}".format(window_size))
	axs[4].set_title("Smoothing window = {}".format(window_size))

	axs[0].set_xlabel("Timestep")
	axs[1].set_xlabel("Timestep")
	axs[2].set_xlabel("Average position")
	axs[3].set_xlabel("Timestep")
	axs[4].set_xlabel("Average position")

	axs[0].set_ylabel("Average position")
	axs[1].set_ylabel("First order instantaneous velocity")
	axs[2].set_ylabel("First order instantaneous velocity")
	axs[3].set_ylabel("First order instantaneous velocity")
	axs[4].set_ylabel("First order instantaneous velocity")

	axs[0].legend()
	axs[1].legend()
	axs[2].legend()
	axs[3].legend()
	axs[4].legend()

	fig.set_size_inches((8.6, 30))
	plt.suptitle("Average across each persistence length")
	plt.tight_layout()
	plt.savefig("plots/average_kymo.png", dpi = 150)
	plt.cla()
	plt.clf()

	# cumulative of ALL SMCS:

	tdf = pd.DataFrame(columns = ["time", "avg_pos",])

	for t in times:
		time_filter = (df["time"] == t)
		meanmean_x2 = (df.loc[(time_filter), "x2"].mean())
		new_df = pd.DataFrame([{"time": t, "avg_pos": meanmean_x2}])
		tdf = pd.concat([tdf, new_df], axis = 0, ignore_index = True)

	fig, axs = plt.subplots(3)
	tdf.loc[:, "v_1"] = (tdf.loc[:, "avg_pos"] - tdf.loc[:, "avg_pos"].shift(1))/1000
	tdf.loc[:, "smoothed_v_1"] = tdf.loc[:, "v_1"].rolling(window_size).mean()
	
	axs[0].plot(tdf["time"], tdf["avg_pos"], "k-")
	axs[1].plot(tdf["time"], tdf["v_1"], "r-", label = "unsmoothed v_1")
	axs[1].plot(tdf["time"], tdf["smoothed_v_1"], "b-", label = "smoothed v_1")
	axs[2].plot(tdf["avg_pos"], tdf["v_1"], "r-", label = "unsmoothed v_1")
	axs[2].plot(tdf["avg_pos"], tdf["smoothed_v_1"], "b-", label = "smoothed v_1")

	axs[0].set_xlabel("Timestep")
	axs[1].set_xlabel("Timestep")
	axs[2].set_xlabel("SMC position")

	axs[0].set_ylabel("SMC position")
	axs[1].set_ylabel("First order velocity")
	axs[2].set_ylabel("First order velocity")

	axs[1].legend()
	axs[2].legend()

	fig.set_size_inches((8.6, 12))
	plt.suptitle("Average across all SMC (across all persistence lengths)")
	plt.tight_layout()
	plt.savefig("plots/average_all_smc_kymo.png", dpi = 150)
	plt.cla()
	plt.clf()


	avg_dist = []
	for lp in lps:
		lp_filter = (adf["lp"] == lp) & (adf["time"] == 100000)
		avg_dist.append(adf.loc[lp_filter, "avg_pos"].values[0])
	
	plt.plot(lps, avg_dist, "k.")
	plt.xlabel("persistence length")
	plt.ylabel("average position over entire simulation")
	plt.gcf().set_size_inches((8.6, 6))
	plt.tight_layout()
	plt.savefig("plots/pos_vs_lp.png", dpi = 300)
	plt.cla()
	plt.clf()



def plotter():
	df = pd.read_csv("~/Documents/tap/smc-single-polymer/test_scan_lp/merged_2.csv")

	# get a few colours
	nlp = df["persistence_length"].nunique()
	cmap = matplotlib.cm.get_cmap('jet')
	norm = matplotlib.colors.Normalize(vmin = 0, vmax = 1)
	colors = [cmap(x) for x in np.linspace(0, 1, nlp )]

	fig, axs = plt.subplots(2, nlp)
	sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

	lp_counter = 0
	lp_ticks = []

	for lp, _ in df.groupby("persistence_length"):
		axs[0, lp_counter].set_title(f"lp = {lp}")
		lp_ticks.append("{:.2f}".format(lp))

		vel_dfs = []

		for rep, _df in _.groupby("replica_id"):
			_df = shift_start(_df)
			_vel_df = compute_velocity(_df)[["time", "smoothed_v_1"]]

			axs[0, lp_counter].plot(_df["time"], _df["x1"], color = colors[lp_counter], alpha = 0.4)
			axs[0, lp_counter].plot(_df["time"], _df["x2"], color = colors[lp_counter], alpha = 0.4)
			axs[0, lp_counter].set_ylim(-50, 50)

			vel_dfs.append(_vel_df)

		vel_df = pd.concat(vel_dfs, ignore_index = True).sort_values("time")

		win_len = max(2, len(vel_df)//50)

		smoothed = signal.savgol_filter(vel_df["smoothed_v_1"], win_len, 1)

		axs[1, lp_counter].plot(vel_df["time"],vel_df["smoothed_v_1"], ".", color = colors[lp_counter], alpha = 0.2)

		axs[1, lp_counter].plot(vel_df["time"], smoothed, "-", color = colors[lp_counter])

			# plt.plot(_df["x1"])
		

		lp_counter += 1


	# cbar = plt.colorbar(sm)
	# cbar.set_ticks(np.linspace(0, 1, nlp))
	# cbar.set_ticklabels(lp_ticks)
	fig.set_size_inches((18, 12))
	axs[0, 0].set_ylabel("position from start")
	axs[0, 2].set_xlabel("Simulation time")
	axs[1, 0].set_ylabel("first order velocity")
	axs[1, 2].set_xlabel("Simulation time")

	plt.tight_layout()
	plt.savefig("test_vary_lp.png", dpi = 300)
	plt.show()

	# then plot velocity curves
	# print(df)

def analysis():
	pass


if __name__ == "__main__":
	parse()
	# plotter()
	# average_plotter()