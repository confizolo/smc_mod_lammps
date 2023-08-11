import os

def generate_stiff(sl: int, molecule_folder, molecule_file, force = False, extend_boundary = False, pattern = "",):
	"""
	format: 1.1 1 bead, 1 blank
	count such that  there is one bead at the end (fencepost counting)
	"""

	def parse_pattern(sl, x):
		bunch, gap = x.split(".")
		bunch, gap = int(bunch), int(gap)

		_sl = sl + ((sl/bunch)-1) * gap
		_start = int((lpol - _sl)/2 + 1) # if sl = 0, this is 501
		_end = int(stiff_start + _sl) # if sl = 0, this is 501

		_build_str = ""

		while _sl > 0:
			_build_str += "o" * bunch
			_sl -= bunch
			_build_str += "x" * gap

		return _start, _end, _build_str
	
	# MOLECULE_FILE = "/home/zy/Documents/tap/smc-lammps/data/In_conf.Nb1000.RW4.fixed.dat"
	MOLECULE_FILE = molecule_file
	NEW_ATOM = 4
	NEW_ANGLE = 2
	# sl = 500 # this is variable

	text = []
	with open(MOLECULE_FILE, "r") as f:
		for line in f:
			text.append(line)	

	text[3] = "4 atom types\n"
	text[7] = "2 angle types\n"

	if extend_boundary:
		text[9] = "-{:.2f} {:.2f} xlo xhi\n".format(extend_boundary, extend_boundary)

	text.insert(18, "4 1\n")

	if not os.path.exists(molecule_folder):
		os.makedirs(molecule_folder)

	output_path = os.path.join(molecule_folder, f"stiff_{sl}{pattern}.dat")
	lpol = 1000

	if os.path.exists(output_path) and not force:
		return os.path.abspath(output_path)

	# centre the stiff region
	if len(pattern):
		stiff_start, stiff_end, build_str = parse_pattern(sl, pattern)	

	else:
		stiff_start = int((lpol - sl)/2 + 1) # if sl = 0, this is 501
		stiff_end = int(stiff_start + sl) # if sl = 0, this is 501
		build_str = "o" * sl

	read_flag = False
	atom_data = []
	angle_data = []

	for c, line in enumerate(text):
		if "Velocities" in line: 
			read_flag = False
			atom_text_end = c 

		if read_flag and len(line.strip()):
			_data = line.split(" ")
			_id = int(_data[0])
			if _id in range(stiff_start, stiff_end):
				if build_str[_id - stiff_start] == "o":
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
			if all([build_str[int(_data[j].strip())-stiff_start] == "o" for j in [-1, -2, -3]]):
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