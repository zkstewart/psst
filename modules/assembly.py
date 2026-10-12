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

def assembly_hifiasm(args):
    '''
    Generates (and optionally submits) a PBS job script which runs a hifiasm genome
    assembly, followed by FASTA conversion, statistics and completeness scoring. The
    Hi-C template is used if --hic was given (with or without files).

    Parameters:
        args -- an argparse.Namespace containing validated "assembly hifiasm" arguments.
    '''
    values = {
        "genScriptDir": args.genScriptDir,
        "compleasmDir": args.compleasmDir,
        "cpus": args.cpus,
        "lineage": args.lineage
    }
    if args.hicFiles == None:
        finalise_job(args, "run_hifiasm_assembly.sh", values)
    else:
        hic1, hic2 = args.hicFiles if args.hicFiles != [] else (None, None)
        finalise_job(args, "run_hifiasm_assembly.hic.sh", {
            **values,
            "ploidy": args.ploidy,
            "hic1": hic1,
            "hic2": hic2
        })

def assembly_revcomp(args):
    '''
    Generates (and optionally submits) a PBS job script which reverse complements
    scaffolded haplotypes to match the orientation of a reference genome.

    Parameters:
        args -- an argparse.Namespace containing validated "assembly revcomp" arguments.
    '''
    finalise_job(args, "run_revcomp.sh", {
        "varScriptDir": args.varScriptDir,
        "refGenome": args.refGenome,
        "cpus": args.cpus,
        "outputDir": args.outputDir,
        "preset": args.preset
    })

def assembly_scaffolding(args):
    '''
    Generates (and optionally submits) a PBS job script which scaffolds haplotype
    assemblies against a diploid reference genome.

    Parameters:
        args -- an argparse.Namespace containing validated "assembly scaffolding" arguments.
    '''
    finalise_job(args, "run_scaffolding_pipeline.sh", {
        "varScriptDir": args.varScriptDir,
        "paf2dotplotDir": args.paf2dotplotDir,
        "refHap1": args.refHap1,
        "refHap2": args.refHap2,
        "cpus": args.cpus,
        "outputDir": args.outputDir,
        "preset": args.preset
    })
