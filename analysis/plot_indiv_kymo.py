from scipy.stats import linregress
import matplotlib
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

def filter_generator(df, lps, get_lp_id = False):
	for c, lp in enumerate(lps):
		lp_filter = (df["persistence_length"] == lp)
		# get reps
		reps = df.loc[lp_filter, "replica_id"].unique()

		for rep in reps:
			if get_lp_id:
				yield c, (df["replica_id"] == rep) & lp_filter
			else:
				yield (df["replica_id"] == rep) & lp_filter


def compute_fixed_pos_velocity(input_file, prefix):

	df = pd.read_csv(input_file)
	
	# df = df.loc[df["replica_id"] == "test_scan_lp/lp05/7"]

	cmap = matplotlib.colormaps['jet']
	norm = matplotlib.colors.Normalize(vmin = 0, vmax = 1)
	nlp = df["persistence_length"].nunique()
	colors = [cmap(x) for x in np.linspace(0, 1, nlp )]

	fig, axs = plt.subplots(2, nlp)
	sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

	lps = df["persistence_length"].unique()
	# times = df["time"].unique()

	df = df.sort_values(["time","persistence_length", "replica_id"])

	# modify the dataframe in place to get an inst velocity 

	for x in filter_generator(df, lps):
		start_time = (df.loc[x].reset_index().loc[0, "time"])

		x1_0, x2_0 =  (df.loc[x & (df["time"] == start_time), ["x1", "x2"]]).values[0]

		df.loc[x, "x1"] -= x1_0
		df.loc[x, "x2"] -= x2_0
		df.loc[x, "ll"] = df.loc[x, "x2"] - df.loc[x, "x1"]

		for c, row in df.loc[x].iterrows():
			df.at[c, "hash"] = "{}_{}".format(row.replica_id, row.ll)
		
	vdf = df.drop_duplicates("hash")


	for x in filter_generator(vdf, lps):

		vdf.loc[x,"time_delta"] = vdf.loc[x, "time"] - vdf.loc[x, "time"].shift(1)

		vdf.loc[x, "v_1"] = 1/vdf.loc[x, "time_delta"]
		vdf.loc[x, "v_2"] = 1.5/vdf.loc[x, "time_delta"] - 0.5/vdf.loc[x, "time_delta"].shift(1)

	vdf.to_csv(f"~/Documents/tap/smc-single-polymer/test_scan_lp/{prefix}.indiv_inst_v.csv", index = False)

	# write back to original, fill forward
	# for each row in the reduced dataframe
	# take the match the hash and time to the original dataframe
	# then forward fill for each replica slice (could be funky if it's all in the same view)
	# also backfill with 0 for like. consistency i guess?

	# for c, row in vdf.iterrows():
	# 	print(c)
	# 	for y in ["v_1", "v_2",]:
	# 		# super expensive operation
	# 		df.loc[(df["time"] == row.time) & (df["hash"] == row.hash), y] = row[y]

	# for x in filter_generator(df, lps):
	# 	for y in ["v_1", "v_2",]:
	# 		df.loc[x, y] = df.loc[x, y].ffill()

	# # print(df)
	# df.to_csv(f"~/Documents/tap/smc-single-polymer/test_scan_lp/{prefix}.rewrite.indiv_inst_v.csv", index = False)
	



def plot_inst_velocities(input_file, prefix):
	df = pd.read_csv(input_file)
	lps = df["persistence_length"].unique()

	####
	# position-time scatterplot
	###
	# velocity-time scatterplot
	###
	# velocity-position scatterplot
	###
	# smoothed velocity-time scatterplot
	###
	# smoothed velocity-position scatterplot

	fig, axs = plt.subplots(len(lps), 4)
	# for lp_id, x in filter_generator(df, lps, get_lp_id = True):

	for lp_id, lp in enumerate(lps):
		x = (df["persistence_length"] == lp)

		axs[lp_id, 0].plot(df.loc[x, "time"], df.loc[x, "x2"], "k.", alpha = 0.05)
		axs[lp_id, 0].set_xlabel("Simulation timestep")
		axs[lp_id ,0].set_ylabel("Inst. Bead position")
		axs[lp_id, 0].set_title("lp = {:.2f} | Position vs time".format(lps[lp_id]))

		axs[lp_id, 1].plot(df.loc[x, "time"], df.loc[x, "v_1"], "k.", alpha = 0.05)
		axs[lp_id, 1].set_xlabel("Simulation timestep")
		axs[lp_id ,1].set_ylabel("Inst. velocity")
		axs[lp_id, 1].set_title("lp = {:.2f} | Velocity vs time".format(lps[lp_id]))


		# axs[lp_id, 2].hexbin(df.loc[x, "x2"], df.loc[x, "v_1"], gridsize = 20, cmap = 'jet')

		# slope, intercept, r, p, se = linregress(df.loc[x, "x2"], df.loc[x, "v_1"])

		# plot linear regress

		# _X = np.linspace(df.loc[x,"x2"].min(), df.loc[x,"x2"].max(), 2)
		# _Y = slope * _X + intercept

		# axs[lp_id, 2].plot(_X, _Y, "r-")


		axs[lp_id, 2].plot(df.loc[x, "x2"], df.loc[x, "v_1"], "b.", alpha = 0.02, markersize = 12)
		axs[lp_id, 2].set_xlabel("Inst. Bead position")
		axs[lp_id ,2].set_ylabel("Inst. velocity")
		axs[lp_id, 2].set_title("lp = {:.2f}".format(lps[lp_id]))

		# axs[lp_id, 2].set_title("lp = {:.2f} | slope = {:.2E}| R^2 = {:.2f}".format(lps[lp_id], slope, r**2))


		axs[lp_id, 3].hexbin(df.loc[x, "x2"], df.loc[x, "v_1"], cmap = "binary", gridsize = 5)
		axs[lp_id, 3].set_xlabel("Inst. Bead position")
		axs[lp_id, 3].set_ylabel("Inst. velocity")
		axs[lp_id, 3].set_title("lp = {:.2f}".format(lps[lp_id]))

	fig.set_size_inches((16, len(lps) * 3))
	plt.tight_layout()
	plt.savefig(f"plots/{prefix}.position-velocity-time.png", dpi = 300)
	plt.show()
	plt.clf()
	plt.cla()


	## plot histograms of inst velocity
	# is this aggregated?
	# as a function of position? time

def hist_inst_vel(input_file, prefix):
	df = pd.read_csv(input_file)
	lps = df["persistence_length"].unique()

	cmap = matplotlib.colormaps['viridis']
	nlp = df["persistence_length"].nunique()
	colors = [cmap(x) for x in np.linspace(0, 1, nlp )]

	fig, axs = plt.subplots(len(lps), 2)


	for x in filter_generator(df, lps):
		df.loc[x, "smoothed_v_1"] = df.loc[x, "v_1"].rolling(4).mean()

	for lp_id, lp in enumerate(lps):
		x = (df["persistence_length"] == lp)

		axs[lp_id, 0].hist(df.loc[x, "v_1"], color = colors[lp_id], edgecolor = "k")
		axs[lp_id, 1].hist(df.loc[x, "smoothed_v_1"], color = colors[lp_id], edgecolor = "k")

		axs[lp_id, 0].set_title(f"v_1 | lp = {lp:.2f}")
		axs[lp_id, 1].set_title(f"smoothed | v_1 lp = {lp:.2f}")

		axs[lp_id, 0].set_ylabel("No. of occurences")
		axs[lp_id, 0].set_xlabel("Instantaneous velocity")
		axs[lp_id, 1].set_ylabel("No. of occurences")
		axs[lp_id, 1].set_xlabel("Instantaneous velocity")


	fig.set_size_inches((12, len(lps) * 3))
	plt.tight_layout()
	plt.savefig(f"plots/{prefix}.inst-vel-hist.png", dpi = 300)
	plt.show()
	plt.clf()
	plt.cla()

	


if __name__ == "__main__":
	compute_fixed_pos_velocity("~/Documents/tap/smc-single-polymer/longer_run_96.csv", "longer_96")

	plot_inst_velocities("~/Documents/tap/smc-single-polymer/test_scan_lp/longer_96.indiv_inst_v.csv", "longer_96")

	# hist_inst_vel("~/Documents/tap/smc-single-polymer/test_scan_lp/longer_96.rewrite.indiv_inst_v.csv", "longer_96")