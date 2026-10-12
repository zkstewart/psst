# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os, shlex, subprocess, sys

from .pbs import finalise_job

def basespace_download(args):
    '''
    Generates (and optionally submits) a PBS job script which downloads a BaseSpace
    project with the bs command line tool.

    Parameters:
        args -- an argparse.Namespace containing validated "basespace download" arguments.
    '''
    finalise_job(args, "run_basespace_download.sh", {
        "bsDir": args.bsDir,
        "projectID": args.projectID,
        "outputDir": args.outputDir
    })

def basespace_prepare(args):
    '''
    Runs Various_scripts/Fasta_related/concatenate_trimmed_lanes.py, which writes a
    shell script to concatenate reads split across sequencing lanes. Unlike other
    modes, this does not write a job script from a template; it runs immediately.

    Parameters:
        args -- an argparse.Namespace containing validated "basespace prepare" arguments.
    '''
    cmd = [sys.executable, os.path.join(args.varScriptDir, "Fasta_related", "concatenate_trimmed_lanes.py"),
           "-d", args.readsDir,
           "-rs", args.readsSuffix,
           "-o", args.outputFileName,
           "--laneIdentifier", args.laneIdentifier,
           "--rpid", args.rpid]
    if not args.symlink:
        cmd.append("--nosymlink")
    if args.singleEnd:
        cmd.append("--singleEnd")

    # Create the output directory, which concatenate_trimmed_lanes.py requires to exist
    os.makedirs(os.path.dirname(args.outputFileName), exist_ok=True)

    print("Running: " + shlex.join(cmd), flush=True) # flush so output stays in order with the script's
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise RuntimeError(f"concatenate_trimmed_lanes.py failed with exit code {result.returncode}")

def basespace_archive(args):
    '''
    Generates (and optionally submits) a PBS job script which archives a directory
    of reads into a .tar file with an accompanying md5 checksum.

    Parameters:
        args -- an argparse.Namespace containing validated "basespace archive" arguments.
    '''
    finalise_job(args, "archive_download.sh", {
        "readsDir": args.readsDir,
        "prefix": args.prefix
    })
