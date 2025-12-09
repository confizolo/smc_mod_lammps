# Test Files

This directory contains two sample scripts used to run PPA analysis on simulation of large-scale solutions of polymers.

## Algorithm Overview

Primitive Path Analysis (PPA) follows the protocol popularized by Everaers et al. to reveal the entanglement network:

1. **Freeze chain termini** – anchor beads are immobilized to lock chain length.
2. **Remove intrachain excluded-volume forces** – only bonds and angular constraints remain so chains can collapse onto their primitive paths.
3. **Integrate with strong damping** – a Langevin thermostat at low target temperature drives the system toward the primitive-path state.
4. **Equilibrate under constrained dynamics** – a short run (here 5×10⁵ steps) relaxes chains while keeping the constraints active.

### Script: `run.digest.ppa.lam`

Designed for bulk systems containing many polymers. Key steps:

- **Group definition** – terminal bead indices are selected for every polymer (e.g., every 1000th bead in a long chain) and assigned to `startbead` / `endbead`.
- **Freezing termini** – the `fix setforce` command (named `ppafix`) nulls forces on terminal beads, immobilizing them in space.

### Script: `run.digest.ppa.noloops.lam`

Targeted to systems where loops (marked as type 2) must be removed before PPA. Recent refactor makes the file suitable for two 10-bead polymers used in validation tests.

- **Automated indexing** – variables compute bead indices for each polymer so the input remains valid when chain length or polymer count changes.
- **Loop removal** – `delete_atoms group loops bond yes` removes type-2 atoms and associated bonds, leaving only backbone beads.
- **Termini locking** – the first and last bead of each 10-bead polymer are grouped and immobilized using `fix setforce`, mirroring the behaviour of the larger production script.

## Practical Notes

- Generate restart/data files with atom types that distinguish loops (type 2) from backbone (type 1) before running the `noloops` variant.
- Both scripts expect an outer driver to define variables such as `outfold`, `rname`, `rep`, `time`, and random seeds.

## References
- Ralf Everaers et al., Rheology and Microscopic Topology of Entangled Polymeric Liquids. Science303,823-826(2004).DOI:10.1126/science.1091215