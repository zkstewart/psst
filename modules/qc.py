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

from .mapping import get_array_size
from .pbs import finalise_job

def qc_library_trim(args):
    '''
    Generates (and optionally submits) a PBS array job script which trims adapters and
    low-quality ends from reads with bbduk. Each subjob trims one sample; the array size
    is a placeholder if no reads directory was given, and no array is requested if there
    is only one sample.

    Parameters:
        args -- an argparse.Namespace containing validated "qc library trim" arguments.
    '''
    finalise_job(args, "qc/run_bbduk_trim.sh", {
        "condaEnv": args.condaEnv,
        "readsDir": args.readsDir,
        "readsSuffix": args.readsSuffix,
        "paired": "false" if args.singleEnd else "true",
        "outputDir": args.outputDir,
        "adapters": args.adapters,
        "trimq": args.trimq,
        "minLength": args.minLength,
        "cpus": args.cpus,
        "xmx": args.xmx
    }, arraySize=get_array_size(args.numSamples))
