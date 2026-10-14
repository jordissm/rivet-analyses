# Rivet analyses

Run from `rivet-analyses/` with Rivet and `rivet-build` available on PATH:

```bash
bash build-analyses.sh electra
```

The script builds each `.cc` into `build/electra/Rivet_<source-name>.so`.
Use a second argument to choose the output root. The `pdg/` directory contains
additional analyses that can be built in the same way.

The profile-specific runners require the SMASH event directory, the aggregated
eHIJING metadata JSONL file and the matching plugin library. They keep outputs
separate by cross-section profile. Read each runner's arguments and binning
defaults before use; the supplied runners target Rivet 3.1.7.

See the [manual's Rivet guide](../docs/analysis/rivet.md) for an example and
the [SLURM guide](../docs/running/slurm.md) for the separate analysis submitter.

