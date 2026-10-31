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

def map_wgs(args):
    '''
    Generates (and optionally submits) a PBS array job script which maps short DNA reads
    with bwa mem, piping its output through samtools to sort and mark duplicates. Each
    subjob maps one sample; the array size is a placeholder if no reads directory was
    given, and no array is requested if there is only one sample.

    Parameters:
        args -- an argparse.Namespace containing validated "map wgs" arguments.
    '''
    if args.numSamples == None:
        arraySize = "@@numSamples@@" # left for the user to fill in
    else:
        arraySize = args.numSamples if args.numSamples > 1 else None
    finalise_job(args, "map/run_bwa_wgs.sh", {
        "genome": args.genome,
        "readsDir": args.readsDir,
        "readsSuffix": args.readsSuffix,
        "paired": "false" if args.singleEnd else "true",
        "outputDir": args.outputDir,
        "library": args.library,
        "cpus": args.cpus,
        "sortMem": args.sortMem
    }, arraySize=arraySize)
