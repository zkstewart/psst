# psst: notes for working on this repository

psst (`psst.py`) generates PBS job scripts for the QUT HPC from bash templates, and can submit them with `qsub`. It consolidates older scripts from the [Various_scripts](https://github.com/zkstewart/Various_scripts) repository, several of which psst jobs still call (so changes may need to land there first; commit Various_scripts before psst when both change).

## Layout
- `psst.py`: the whole CLI (nested argparse subparsers), shared argument helpers, and per-mode dispatch (`<x>main()` functions calling validators).
- `modules/validation.py`: all validation, as `validate_<abbrev>` (parent mode) and `validate_<abbrev>_<submode>` functions, plus shared helpers (`validate_input_file`, `validate_program_dir`, `validate_reads`, `validate_pbs_resources`, `parse_mem_gb`, ...).
- `modules/<mode>.py`: one module per parent mode (e.g., `mapping.py` for `map`); each leaf function fills a template and calls `finalise_job()`.
- `modules/pbs.py`: template filling, PBS header, writing and submitting scripts.
- `modules/reads.py`: read file discovery and pair checking (`get_read_files`).
- `modules/config.py` / `modules/locations.py`: program/database locations from `default_config.json`, or the untracked `my_config.json` if present (must have the same keys; `rpid_shortcuts` is freeform).
- `templates/pbs_header.sh` (shared) and `templates/<mode>/*.sh` (one directory per parent mode).

## Adding a mode
1. Parser in `psst.py`: `add_parser(..., parents=[p], add_help=False)`, `set_defaults(func=..., scriptSuffix=...)`, then `add_pbs_arguments(parser, jobName, walltime, mem, cpus)` last. New parent modes get a subparser with `dest="<name>Mode"`, a `validate_<abbrev>()` (often `pass`), a dispatch block in `main()` and a `<abbrev>main()` function; a parent mode without submodes is itself the leaf parser, and its dispatch block calls `validate_<abbrev>()` directly.
2. Validator in `modules/validation.py`, imported and called from the dispatch function.
3. Job function in `modules/<mode>.py` passing a dict of every template field to `finalise_job(args, "<mode>/<template>.sh", values[, arraySize=...])`.
4. Template in `templates/<mode>/`.
5. New program locations: add a key to `default_config.json` and a constant to `modules/locations.py`; remind the user to add it to their HPC `my_config.json`.

## Argument conventions (also summarised at the top of `psst.py`)
- Only `--rpid` is required (full RPID or a config shortcut). Everything else has a default or defaults to `None`, which leaves an `@@field@@` placeholder in the script; `--qsub` is refused while placeholders remain.
- Help text uses `%(default)s`; placeholder defaults say "default == left as a placeholder".
- Location arguments use `action=UserSet`: validation only checks values the user typed (`args.userSet`), never default HPC paths. Exceptions are modes that run immediately rather than writing a job (e.g., `basespace prepare`).
- Derived defaults are set in validation (e.g., an index directory from `-r`, a prefix from the genome name).
- Dests must not reuse routing names (`mode`, subparser dests, `func`, `scriptSuffix`, `pbsMode`, `userSet`); `check_reserved_dests()` enforces this.
- Shared argument sets have helpers (`add_reads_arguments`, `add_bwa_arguments`, `add_varscripts_argument`, ...); reuse them, and give a mode's validation a matching shared validator.
- Walltime is capped at 48h. Mapping modes that need an index warn rather than fail on a missing index when `--afterok` is given, so jobs can be chained.

## Templates and job scripts
- Templates are plain bash with `@@field@@` placeholders (not `str.format`, so bash braces are safe). `fill_template` errors if template fields and supplied keys differ; values are `shlex`-quoted, so pass flags as booleans (`"true"`/`"false"`) and assemble argument arrays in bash, not as one string. An empty string disables an optional input.
- Layout: a `####` block of commented variable assignments, a closing `####`, then `# STEP N:` sections. Use `set -eo pipefail` (after any `module load`/`conda activate`).
- Never overwrite outputs: psst refuses existing scripts/outputs, and jobs check again at runtime. Write outputs under a temporary `.incomplete` name and `mv` on success so partial output never looks complete; per-sample jobs skip samples whose final output exists.
- Per-sample work is a PBS array job: `get_array_size(numSamples)` gives `#PBS -J 1-N` (no array for one sample, a `@@numSamples@@` placeholder if reads weren't given); templates use `${PBS_ARRAY_INDEX:-1}` and rediscover read files with the same rules as `modules/reads.py` (paired reads end in `1`/`2` before the suffix; the sample name strips a trailing `[._-]R?`).
- Multi-job modes (`basecalling dorado`, `call *`) use `add_pbs_resource_arguments(step=...)` per job plus `add_pbs_submission_arguments()`, render every script before writing any, and with `--qsub` chain later jobs via `submit_job(..., afterok=...)`. Genome-wide work is split into chunks of whole contigs (`call`: a `.fai`-based awk partition in genome order, with `.empty` markers for chunks without contigs).
- Index directories made by psst use fixed file names that later jobs rely on (e.g., `map index giraffe` writes `pangenome.gbz`, `reference.fa`, `sites.vcf.gz`, ...; validation checks them via `GIRAFFE_MAP_FILES`/`GIRAFFE_CALL_FILES` in `mapping.py`), so later modes take just `-i <index dir>`.
- Prefer streaming pipelines over intermediate files, and give tools memory limits derived from `--mem` where they would otherwise size themselves from the whole node (e.g., salmon `--ramLimit`, STAR `--limitGenomeGenerateRAM`).

## Testing
There is no test suite. Verify changes by generating scripts in the scratchpad (placeholders, every validation error path, `--qsub` refusal; a fake `qsub` on PATH tests submission) and, where possible, running the generated job end to end with real tools (a scratch conda env from bioconda works; stub `module` for module-loaded tools) on small simulated data; export `PBS_O_WORKDIR` to the test directory first, since the header's `cd $PBS_O_WORKDIR` otherwise runs jobs in `$HOME`. Afterwards, confirm every existing mode still generates with only `--rpid`.

## Repository practices
- Commit and push only when asked; commit messages describe the change in a short paragraph and end with the attribution line.
- Keep `README.md` brief; don't document niche options in depth there.
