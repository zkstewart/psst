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
        finalise_job(args, "assembly/run_hifiasm_assembly.sh", values)
    else:
        hic1, hic2 = args.hicFiles if args.hicFiles != [] else (None, None)
        finalise_job(args, "assembly/run_hifiasm_assembly.hic.sh", {
            **values,
            "ploidy": args.ploidy,
            "hic1": hic1,
            "hic2": hic2
        })

def phase_values(args):
    '''
    Parameters:
        args -- an argparse.Namespace containing validated arguments of an "assembly" mode
                which phases contigs with annotarium (see validation.validate_phasing()).
    Returns:
        values -- a dict of the template values shared by these modes.
    '''
    return {
        "annotariumDir": args.annotariumDir,
        "minimap2Dir": args.minimap2Dir,
        "contigFiles": args.contigFiles,
        "refGenome": args.refGenome,
        "phaseDir": args.phaseDir if args.phaseDir != None else "", # empty phases the contigs
        "ploidy": args.ploidy,
        "preset": args.preset,
        "outputDir": args.outputDir,
        "cpus": args.cpus
    }

def assembly_revcomp(args):
    '''
    Generates (and optionally submits) a PBS job script which phases contigs into haplotype
    partitions against a reference genome with annotarium (unless given an earlier phasing),
    then writes each partition's contigs paired 1:1 with the reference contigs, renamed and
    reverse complemented to match them.

    Parameters:
        args -- an argparse.Namespace containing validated "assembly revcomp" arguments.
    '''
    finalise_job(args, "assembly/run_revcomp.sh", phase_values(args))

def assembly_scaffolding(args):
    '''
    Generates (and optionally submits) a PBS job script which phases contigs into haplotype
    partitions against a reference genome with annotarium (unless given an earlier phasing),
    then scaffolds each partition against the reference with RagTag and dotplots it.

    Parameters:
        args -- an argparse.Namespace containing validated "assembly scaffolding" arguments.
    '''
    finalise_job(args, "assembly/run_scaffolding.sh", {
        **phase_values(args),
        "paf2dotplotDir": args.paf2dotplotDir
    })
