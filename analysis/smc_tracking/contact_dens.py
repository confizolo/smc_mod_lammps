
from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple

import numpy as np
from tqdm import tqdm

from numba import njit  
from read_files import read_lammpstrj


@njit
def count_contacts(frame: np.ndarray, cutoff: float = 0.2) -> Tuple[int, int]:
    """Count intra- and inter-polymer contacts for one frame.

    Parameters
    ----------
    frame:
        Array of shape ``(N, 5)`` containing ``x, y, z, type, mol`` columns.
    cutoff:
        Euclidean distance threshold that defines a contact.

    Returns
    -------
    (int, int)
        ``(intra_contacts, inter_contacts)``
    """

    intra = 0
    inter = 0
    for i in range(frame.shape[0]):
        for j in range(i):
            dist = np.linalg.norm(frame[i, :3] - frame[j, :3])
            if dist < cutoff:
                if frame[i, 4] == frame[j, 4]:
                    intra += 1
                else:
                    inter += 1
    return intra, inter

def analyse_contacts(
    trajectory: str,
    output_dir: str,
    *,
    cutoff: float = 0.2,
    bead_types: Optional[Sequence[int]] = None,
    bead_per_polymer: int = 1,
    frame_start: Optional[int] = None,
    frame_stop: Optional[int] = None,
    frame_step: int = 1,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute intra- and inter-polymer contact counts for a trajectory.
    Parameters
    ----------  
    trajectory:
        Path to LAMMPS trajectory file.
    output_dir:
        Directory to save output contact count files.
    cutoff:
        Euclidean distance threshold that defines a contact.
    bead_types:
        If provided, only consider beads of these types for contact counting.
    bead_per_polymer:
        Number of beads per polymer; used to determine polymer indices (SMC filtering).
    frame_start:
        First frame to analyze (default: first frame).
    frame_stop:
        Last frame to analyze (default: last frame).
    frame_step:
        Step size between frames to analyze.    
    """

    frames = read_lammpstrj(trajectory)
    if not frames:
        raise RuntimeError(f"no frames found in {trajectory}")

    selected = frames[frame_start:frame_stop:frame_step]
    intra_counts = np.zeros(len(selected), dtype=int)
    inter_counts = np.zeros(len(selected), dtype=int)

    for idx, frame in enumerate(tqdm(selected, desc="frames")):
        coords = frame["xyz"]
        types = frame["types"]
        mols = frame["mols"] // bead_per_polymer

        data = np.column_stack((coords, types, mols))
        if bead_types is not None:
            mask = np.isin(types, np.asarray(bead_types))
            data = data[mask]

        intra, inter = count_contacts(data.astype(np.float64), cutoff)
        intra_counts[idx] = intra
        inter_counts[idx] = inter

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savetxt(out_dir / "intra_cnt.txt", intra_counts, fmt="%d")
    np.savetxt(out_dir / "inter_cnt.txt", inter_counts, fmt="%d")

    return intra_counts, inter_counts

