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

from .pbs import finalise_job

def clustering_mmseqs2(args):
    '''
    Generates (and optionally submits) a PBS job script which clusters sequences with
    MMseqs2 (cascaded clustering or linclust) via Various_scripts/run_clustering.py.

    Parameters:
        args -- an argparse.Namespace containing validated "clustering mmseqs2" arguments.
    '''
    finalise_job(args, "run_mmseqs_clust.sh", {
        "varScriptDir": args.varScriptDir,
        "mmseqsDir": args.mmseqsDir,
        "fasta": args.fastaFile,
        "molecule": args.molecule,
        "outputFile": args.outputFile,
        "tmpDir": args.tmpDir,
        "program": args.program,
        "cpus": args.cpus,
        "identity": args.identity,
        "coverage": args.coverage,
        "evalue": args.evalue,
        "mode": args.clusterMode,
        "sensitivity": args.sensitivity,
        "steps": args.steps
    })

def clustering_cdhit(args):
    '''
    Generates (and optionally submits) a PBS job script which clusters sequences with
    CD-HIT via Various_scripts/run_clustering.py.

    Parameters:
        args -- an argparse.Namespace containing validated "clustering cdhit" arguments.
    '''
    finalise_job(args, "run_cdhit_clust.sh", {
        "varScriptDir": args.varScriptDir,
        "cdhitDir": args.cdhitDir,
        "fasta": args.fastaFile,
        "molecule": args.molecule,
        "outputFile": args.outputFile,
        "cpus": args.cpus,
        "identity": args.identity,
        "cdhitMem": args.cdhitMem,
        "align": args.align,
        "shorterCov": args.shorterCov,
        "longerCov": args.longerCov
    })
