"""Minimal helpers to read LAMMPS trajectory and bond dump files."""

from typing import Dict, List, Optional, Union

import numpy as np


FrameDict = Dict[str, Union[np.ndarray, int, float]]


def read_lammpstrj_frame(f_in) -> Optional[FrameDict]:
    """Read one ``.lammpstrj`` frame from an opened file object."""
    
    line = f_in.readline()
    if line == "":
        return None

    timestep = int(f_in.readline())
    f_in.readline()
    n_atoms = int(f_in.readline())
    f_in.readline()
    for _ in range(3):
        lohi = f_in.readline().split(" ")
        lo, hi = float(lohi[0]), float(lohi[1])
        L = abs(lo) + abs(hi)
    f_in.readline()

    indexes: List[int] = []
    pos: List[List[float]] = []
    pos_uwr: List[List[float]] = []
    types: List[int] = []
    mols: List[int] = []

    for _ in range(n_atoms):
        line = f_in.readline().split(" ")
        indexes.append(int(line[0]))
        mols.append(int(line[1]))
        types.append(int(line[2]))
        xs, ys, zs = float(line[3]), float(line[4]), float(line[5])
        ix, iy, iz = int(line[6]), int(line[7]), int(line[8])
        pos.append([xs + ix * L, ys + iy * L, zs + iz * L])
        pos_uwr.append([xs, ys, zs])

    sorting = np.argsort(indexes)

    pos = np.array(pos)[sorting, :]
    pos_uwr = np.array(pos_uwr)[sorting, :]

    output: FrameDict = {}
    output["L"] = L
    output["types"] = np.array(types)[sorting]
    output["index"] = np.array(indexes)[sorting]
    output["timestep"] = timestep
    output["mols"] = np.array(mols)[sorting]
    output["xyz"] = pos
    output["xyz_uwr"] = pos_uwr

    return output


def read_lammpstrj(fn_in: str) -> List[FrameDict]:
    """Read the whole LAMMPS trajectory from a file."""

    frames: List[Dict[str, np.ndarray]] = []
    with open(fn_in) as f_in:
        while True:
            frame = read_lammpstrj_frame(f_in)
            if not frame:
                break
            frames.append(frame)
    return frames


def read_bonds_frame(f_in) -> Optional[List[np.ndarray]]:
    """Read one ``.bonds`` frame from an opened file object."""

    line = f_in.readline()
    if line == "":
        return None

    timestep = int(f_in.readline())
    f_in.readline()
    n_bonds = int(f_in.readline())
    f_in.readline()
    for _ in range(3):
        lohi = f_in.readline().split(" ")
        lo, hi = float(lohi[0]), float(lohi[1])
        L = abs(lo) + abs(hi)
    f_in.readline()

    indexes: List[int] = []
    batom1: List[int] = []
    batom2: List[int] = []
    btype: List[int] = []

    for _ in range(n_bonds):
        line = f_in.readline().split(" ")
        indexes.append(int(line[0]))
        batom1.append(int(line[1]))
        batom2.append(int(line[2]))
        btype.append(int(line[3]))

    batoms = np.hstack(
        [
            np.array(batom1).reshape(-1, 1),
            np.array(batom2).reshape(-1, 1),
            np.array(btype).reshape(-1, 1),
        ]
    )
    fil = batoms[:, 1] < batoms[:, 0]
    batoms_old = batoms.copy()
    batoms[fil, 0] = batoms_old[fil, 1]
    batoms[fil, 1] = batoms_old[fil, 0]
    batoms = np.sort(batoms.view("i8,i8,i8"), order=["f1"], axis=0).view(np.int64)
    batoms = np.sort(batoms.view("i8,i8,i8"), order=["f2"], axis=0).view(np.int64)

    return [batoms]


def read_bonds(fn_in: str) -> List[List[np.ndarray]]:
    """Read the whole LAMMPS bonds dump from a file."""

    frames: List[List[np.ndarray]] = []
    with open(fn_in) as f_in:
        while True:
            frame = read_bonds_frame(f_in)
            if not frame:
                break
            frames.append(frame)
    return frames
