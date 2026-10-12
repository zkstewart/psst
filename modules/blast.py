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

def blast_run_mmseqs2(args):
    '''
    Generates (and optionally submits) a PBS job script which runs an MMseqs2
    search via Various_scripts/run_mmseqs2.py.

    Parameters:
        args -- an argparse.Namespace containing validated "blast run mmseqs2" arguments.
    '''
    finalise_job(args, "run_mms2.sh", {
        "mmseqsDir": args.mmseqsDir,
        "varScriptDir": args.varScriptDir,
        "query": args.query,
        "database": args.database,
        "outputPrefix": args.outputPrefix,
        "tmpDir": args.tmpDir,
        "cpus": args.cpus,
        "iterations": args.iterations,
        "sensitivity": f"{args.sensitivity:g}" if args.sensitivity != None else None,
        "searchType": args.searchType
    })
