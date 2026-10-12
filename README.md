# psst: PBS script submission tool
This repository houses the `psst.py` Python tool for coordinating the configuration and submission of job scripts on the QUT HPC. It is created as a personal utility, but if it helps anyone else, that is great.

## What is the point of psst?
The author has a lot of PBS scripts floating around in the [Various_scripts](https://github.com/zkstewart/Various_scripts). Finding anything is a small challenge, and it may reduce mental burden to have a single launching point for all these various scripts. Additionally, the scripts could use some modernisation to make them (1) better written and hence more understandable, and (2) compatible with the new HPC requirements of submitting jobs with Research Project ID (RPID) attached to them. `psst` will make that easier.

## Setup
`psst.py` is executable, so you can symlink it into a directory on your `PATH` to run it as `psst`:

```
ln -s /path/to/psst/psst.py ~/bin/psst
```

Default program, database and reference genome locations, as well as RPID shortcuts, are read from `default_config.json`. To use your own values, copy it to `my_config.json` in the same directory and edit that copy; `my_config.json` is ignored by git and is used in preference to `default_config.json` whenever it exists. It must contain the same sections and keys as `default_config.json`, except for `rpid_shortcuts` which can contain any shortcut names you like.

## How to use
Each job type is a subcommand, and `-h` at any level lists what is available:

```
psst -h
psst assembly -h
psst assembly hifiasm -h
```

### Only `--rpid` is required
Every job requires `--rpid`, which accepts either a full RPID or one of the shortcuts from your config file (case-insensitive). Everything else is optional:

* Arguments with a default value (e.g., program locations, CPUs, memory, walltime) are filled in with that default.
* Arguments without a default (e.g., an input FASTA file or output prefix) are written into the job script as `@@variable@@` placeholders. psst will tell you which placeholders remain, and you should edit the script to replace them before submitting it.

This means that, for example, `psst blast run mmseqs2 --rpid ABCDEF1234` gives you a job script in its default layout, ready for manual editing.

### Submitting with `--qsub`
`--qsub` submits the job script immediately after writing it. It can only be used when every variable has a value, either from the command line or from a default; if any placeholders would remain, psst stops with an error before writing anything.

Jobs are submitted from your current directory, which becomes the job's working directory (`$PBS_O_WORKDIR`). Some jobs expect to run within a particular directory layout (e.g., `assembly hifiasm` reads its input reads from `../fastq_corrected/output`), so run psst from the directory the job should run in. If you submit a script yourself instead, run `qsub` from that same directory.

### Other behaviours to be aware of
* **Existing files are never overwritten.** If a job script already exists at the location psst would write to, psst stops with an error; delete or move the old script first, or use `--script` (or `--workdir` for `basecalling dorado`) to choose a different location.
* **Script locations.** Job scripts are written to `<outprefix>.<job>.sh` if the job takes an output prefix (`-o`) and you gave one, otherwise to `<job>.sh` in the current directory, unless `--script` is given.

## For transparency
This repository will be almost entirely coded using the Claude AI assistant. At the time of writing this README (11th October 2026) I am using Opus 5.5.
