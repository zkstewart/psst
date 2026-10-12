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

def annotable_generate(args):
    '''
    Generates (and optionally submits) a PBS job script which runs an MMseqs2
    search and then produces an annotation table from its results.

    Parameters:
        args -- an argparse.Namespace containing validated "annotable generate" arguments.
    '''
    finalise_job(args, "mms2_table_generation.sh", {
        "mmseqsDir": args.mmseqsDir,
        "varScriptDir": args.varScriptDir,
        "query": args.query,
        "databaseTag": args.databaseTag,
        "database": args.database,
        "idmappingFile": args.idmappingFile,
        "oboFile": args.oboFile,
        "cpus": args.cpus,
        "numHits": args.numHits,
        "iterations": args.iterations,
        "outputPrefix": args.outputPrefix,
        "tmpDir": args.tmpDir
    })
