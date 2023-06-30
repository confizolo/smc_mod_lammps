import matplotlib
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

# plt.style.use('dark_background')

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

def group_average(input_file, prefix):
	df = pd.read_csv(input_file)


	lps = df["persistence_length"].unique()
	times = df["time"].unique()

	df = df.sort_values(["time","persistence_length", "replica_id"])

	# modify the dataframe in place to get an inst velocity 

	for x in filter_generator(df, lps):
		start_time = (df.loc[x].reset_index().loc[0, "time"])

		x1_0, x2_0 =  (df.loc[x & (df["time"] == start_time), ["x1", "x2"]]).values[0]

		df.loc[x, "x1"] -= x1_0
		df.loc[x, "x2"] -= x2_0

	adf = pd.DataFrame()

	for lp in lps:
		lp_filter = (df["persistence_length"] == lp)
		for t in times:
			time_filter = (df["time"] == t)
			x = lp_filter & time_filter

			mean_ll = (df.loc[x, "x2"] - df.loc[x, "x1"]).mean()
			std_ll = (df.loc[x, "x2"] - df.loc[x, "x1"]).std()

			new_df = pd.DataFrame([{
				"time":t,
				"lp":lp,
				"x_avg":mean_ll,
				"x_std":std_ll,
			}])
			adf = pd.concat([adf, new_df], axis = 0, ignore_index = True)

	adf = adf.sort_values(["lp", "time"])

	for c, lp in enumerate(lps):
		x = (adf["lp"] == lp)

		adf.loc[x, "pos_delta"] = adf.loc[x, "x_avg"] - adf.loc[x, "x_avg"].shift(1)

		adf.loc[x, "v_1"] = adf.loc[x, "pos_delta"] #... this is manual

		adf.loc[x, "v_2"] = 1.5 * adf.loc[x, "pos_delta"] - 0.5 * adf.loc[x, "pos_delta"].shift(1)

		window_size = 4

		adf.loc[x, "smoothed_v_1"] = adf.loc[x, "v_1"].rolling(window_size).mean()
		adf.loc[x, "smoothed_v_2"] = adf.loc[x, "v_2"].rolling(window_size).mean()


	adf.to_csv(f"/home/zy/Documents/tap/smc-single-polymer/{prefix}.avg.csv", index = False)


def plot_pos_vel(input_file, prefix):
	df = pd.read_csv(input_file)

	nlp = df["lp"].nunique()
	lps = df["lp"].unique()
	print(nlp)
	print(df)

	cmap = matplotlib.colormaps["viridis"].resampled(nlp)
	norm = matplotlib.colors.Normalize(vmin = min(lps), vmax = max(lps))

	colors = [cmap(x) for x in np.linspace(0, 1, nlp )]

	sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

	# aggregate plots first

	fig, axs = plt.subplots(2, 2)

	for c, lp in enumerate(lps):
		x = (df["lp"] == lp)
		MAX_TIME = df["time"].max()

		axs[0, 0].plot(df.loc[x, "time"], df.loc[x, "x_avg"], "-", color = colors[c])
		axs[0, 0].set_title("position-time")
		axs[0, 0].set_xlabel("Simulation timestep")
		axs[0, 0].set_ylabel("Average replica position")

		axs[1, 0].plot(df.loc[x, "time"], df.loc[x, "x_std"], "-", color = colors[c])
		axs[1, 0].set_title("stdev-time")
		axs[1, 0].set_xlabel("Simulation timestep")
		axs[1, 0].set_ylabel("Stdev replica position")

		axs[0, 1].plot(lp, df.loc[((x) & (df["time"] == MAX_TIME)), "x_avg"].values[0], ".", color = colors[c])
		axs[0, 1].set_xlabel("Persistence length")
		axs[0, 1].set_ylabel("Position at end of simulation")

		axs[1,1].plot(lp, df.loc[((x) & (df["time"] == MAX_TIME)), "x_std"].values[0], ".", color = colors[c])
		axs[1,1].set_xlabel("Persistence length")
		axs[1,1].set_ylabel("Stdev at end of simulation")

	fig.colorbar(sm, location = "top", ax = axs[0, 0])
	fig.colorbar(sm, location = "top", ax = axs[0, 1])

	# refer to https://stackoverflow.com/questions/13784201/how-to-have-one-colorbar-for-all-subplots

	fig.set_size_inches((14, 10))
	plt.tight_layout()
	plt.savefig(f"plots/{prefix}.avg.position-time.png", dpi = 300)
	plt.cla(); plt.clf()

	# then plot velocity time, first and second orders
	# and underneath, smoothed velocity time

	fig, axs = plt.subplots(2, 2)

	for c, lp in enumerate(lps):
		x = (df["lp"] == lp)

		axs[0, 0].plot(df.loc[x, "time"], df.loc[x, "v_1"], "-", color = colors[c], alpha = 0.9)
		axs[0, 0].set_title("1st order velocity-time")
		axs[0, 0].set_xlabel("Simulation timestep")
		axs[0, 0].set_ylabel("Inst. velocity of mean replica position")

		axs[0, 1].plot(df.loc[x, "time"], df.loc[x, "v_2"], "-", color = colors[c], alpha = 0.9)
		axs[0, 1].set_title("2nd order velocity-time")
		axs[0, 1].set_xlabel("Simulation timestep")
		axs[0, 1].set_ylabel("Inst. velocity of mean replica position")
	
		axs[1, 0].plot(df.loc[x, "time"], df.loc[x, "smoothed_v_1"], "-", color = colors[c], alpha = 0.9)
		axs[1, 0].set_title("Smoothed 1st order velocity-time")
		axs[1, 0].set_xlabel("Simulation timestep")
		axs[1, 0].set_ylabel("Inst. velocity of mean replica position")

		axs[1, 1].plot(df.loc[x, "time"], df.loc[x, "smoothed_v_2"], "-", color = colors[c], alpha = 0.9)
		axs[1, 1].set_title("Smoothed 2nd order velocity-time")
		axs[1, 1].set_xlabel("Simulation timestep")
		axs[1, 1].set_ylabel("Inst. velocity of mean replica position")
	fig.colorbar(sm, location = "top", ax = axs[0, 0])
	fig.colorbar(sm, location = "top", ax = axs[0, 1])
	fig.set_size_inches((14, 10))

	plt.tight_layout()
	plt.savefig(f"plots/{prefix}.avg.velocity-time.png", dpi = 300)
	plt.cla(); plt.clf()

	# aggregated velocity position scatter plot
	fig, axs = plt.subplots(4)

	axs[0].hexbin(df["x_avg"], df["v_1"], gridsize = 20, cmap = "inferno")
	axs[1].hexbin(df["x_avg"], df["smoothed_v_1"], gridsize = 20, cmap = "inferno")
	axs[2].hexbin(df["x_avg"], df["v_2"], gridsize = 20, cmap = "inferno")
	axs[3].hexbin(df["x_avg"], df["smoothed_v_2"], gridsize = 20, cmap = "inferno")

	axs[0].set_title("First order velocity")
	axs[1].set_title("Smoothed First order velocity w = 4")
	axs[2].set_title("Second order velocity")
	axs[3].set_title("Smoothed Second order velocity w = 4")

	for ax in axs:
		ax.set_ylim(0, 0.4)

	fig.suptitle("Aggregated velocity-position heatmap")
	fig.set_size_inches((4, 16))

	plt.tight_layout()
	plt.savefig(f"plots/{prefix}.agg.velocity-pos.png", dpi = 300)
	plt.cla(); plt.clf()


	# individual velocity position scatter plot below
	# plot velocity vs pos std because.. reasons
	fig, axs = plt.subplots(nlp, 2)

	for c, lp in enumerate(lps):
		x = (df["lp"] == lp)

		axs[c, 0].plot(df.loc[x, "x_avg"], df.loc[x, "smoothed_v_1"], ".", color = colors[c])
		axs[c, 0].set_xlabel("Average position")
		axs[c, 0].set_ylabel("Smoothed v_1")
		axs[c, 0].set_title(f"lp = {lp:.2f}")

		axs[c, 1].plot(df.loc[x, "x_std"], df.loc[x, "smoothed_v_1"], ".", color = colors[c])
		axs[c, 1].set_xlabel("Stdev position")
		axs[c, 1].set_ylabel("Smoothed v_1")
		axs[c, 1].set_title(f"lp = {lp:.2f}")

	# fig.suptitle("loopsize = 2 * avg_position")
	fig.set_size_inches((12 , nlp * 3))


	plt.tight_layout()
	plt.savefig(f"plots/{prefix}.avg.velocity-pos.png", dpi = 300)
	plt.cla(); plt.clf()


	fig, axs = plt.subplots(2)

	for c, lp in enumerate(lps):
		x = (df["lp"] == lp)
		axs[0].plot(df.loc[x, "x_avg"], df.loc[x, "smoothed_v_1"], ".", color = colors[c])
		axs[0].set_xlabel("Average position")
		axs[0].set_ylabel("Smoothed v_1")
		# axs[0, 0].set_title(f"lp = {lp:.2f}")

		axs[1].plot(df.loc[x, "x_std"], df.loc[x, "smoothed_v_1"], ".", color = colors[c])
		axs[1].set_xlabel("Stdev position")
		axs[1].set_ylabel("Smoothed v_1")
		# axs[0, 1].set_title(f"lp = {lp:.2f}")
	fig.set_size_inches((8.6, 9))
	# fig.suptitle("loopsize = 2 * avg_position")
	fig.colorbar(sm, location = "top", ax = axs[0])
	# fig.colorbar(sm, location = "top", ax = axs[0, 1])

	plt.tight_layout()
	plt.savefig(f"plots/{prefix}.avg.merged.velocity-pos.png", dpi = 300)
	plt.cla(); plt.clf()

	# also plot histograms

	fig, axs = plt.subplots(nlp, 2)

	for c, lp in enumerate(lps):
		x = (df["lp"] == lp)
		axs[c, 0].hist(df.loc[x, "v_1"], color = colors[c], edgecolor = "k")
		axs[c, 1].hist(df.loc[x, "smoothed_v_1"], color = colors[c], edgecolor = "k")

		axs[c, 0].set_xlabel("v_1")
		axs[c, 0].set_ylabel("No. of occurences")
		axs[c, 0].set_title(f"v_1 | lp = {lp:.2f}")

		axs[c, 1].set_xlabel("Smoothed v_1")
		axs[c, 1].set_ylabel("No. of occurences")
		axs[c, 1].set_title(f"Smoothed v_1 w=4 | lp = {lp:.2f}")

	# fig.suptitle("loopsize = 2 * avg_position")
	fig.set_size_inches((16 , nlp * 3))

	plt.tight_layout()
	plt.savefig(f"plots/{prefix}.hist.avg-pos.velocity.png", dpi = 300)
	plt.cla(); plt.clf()



if __name__ == "__main__":
	# group_average("/home/zy/Documents/tap/smc-single-polymer/test_scan_lp/merged_2.csv", "merged_2")
	# plot_pos_vel("/home/zy/Documents/tap/smc-single-polymer/test_scan_lp/merged_2.avg.csv", "merged_2")
	# group_average("/home/zy/Documents/tap/smc-single-polymer/longer_run_96.csv", "longer_96")
	plot_pos_vel("/home/zy/Documents/tap/smc-single-polymer/longer_96.avg.csv", "longer_96")