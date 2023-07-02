import os

def generate_stiff(sl: int, molecule_folder, force = False):
	MOLECULE_FILE = "/home/zy/Documents/tap/smc-lammps/data/In_conf.Nb1000.RW4.fixed.dat"
	NEW_ATOM = 4
	NEW_ANGLE = 2
	assert sl
	# sl = 500 # this is variable

	text = []
	with open(MOLECULE_FILE, "r") as f:
		for line in f:
			text.append(line)	

	text[3] = "4 atom types\n"
	text[7] = "2 angle types\n"
	text.insert(18, "4 1\n")

	if not os.path.exists(molecule_folder):
		os.makedirs(molecule_folder)

	output_path = os.path.join(molecule_folder, f"stiff_{sl}.dat")
	lpol = 1000

	if os.path.exists(output_path) and not force:
		return os.path.abspath(output_path)

	# centre the stiff region
	stiff_start = int((lpol - sl)/2 + 1)
	stiff_end = int(stiff_start + sl)

	read_flag = False
	atom_data = []
	angle_data = []

	for c, line in enumerate(text):
		if "Velocities" in line: 
			read_flag = False
			atom_text_end = c 

		if read_flag and len(line.strip()):
			_data = line.split(" ")
			if int(_data[0]) in range(stiff_start, stiff_end):
				_data[2] = str(NEW_ATOM) 

			atom_data.append(" ".join(_data))

		if "Atoms" in line: 
			read_flag = True
			atom_text_start = c + 1
	
	# not needed but make it explicit
	read_flag = False
	for c, line in enumerate(text):
		if read_flag and len(line.strip()):
			_data = line.split(" ")
			if all([int(_data[j].strip()) in range(stiff_start, stiff_end) for j in [-1, -2, -3]]):
				_data[1] = str(NEW_ANGLE)
			angle_data.append(" ".join(_data))

			pass
		if "Angles" in line: 
			read_flag = True
			angle_text_start = c + 1

	with open(output_path, "w") as f:
		for i in range(atom_text_start):
			f.write(text[i])
		f.write("\n")
		for x in atom_data: 
			f.write(x)
		f.write("\n")
		for i in range(atom_text_end, angle_text_start):
			f.write(text[i])
		f.write("\n")
		for x in angle_data:
			f.write(x)

	return os.path.abspath(output_path)

if __name__ == "__main__":
	generate_stiff(100)