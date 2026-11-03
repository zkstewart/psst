#!/usr/bin/env python3
# A stand-in for annotarium which lets the "assembly" job scripts be tested without real
# data. Every invocation is logged to $STUB_DIR/annotarium.log as one line of arguments
# joined by '|' (so that arguments containing spaces can be told apart), then:
#   assembly phase -- writes phasing.tsv, assigning contig 'c<N>' to each partition N up
#                     to --ploidy, and coverage.tsv to its -o directory
#   assembly homology_pair -- writes '>partition<N>' to its -o file, or fails if N is
#                             $STUB_FAIL_PARTITION
#   assembly scaffold -- writes partition<N>.fasta to its -o directory for each partition
#                        in its -i phasing.tsv

import os, sys

args = sys.argv[1:]
mode = args[1]
def value(flag):
    return args[args.index(flag) + 1]

with open(os.path.join(os.environ["STUB_DIR"], "annotarium.log"), "a") as fileOut:
    fileOut.write("|".join(args) + "\n")

def assigned_partitions(phasingFile):
    rows = [ line.split("\t") for line in open(phasingFile) if not line.startswith("#") ]
    return sorted(set([ row[2] for row in rows if row[3] == "assigned" ]), key=int)

if mode == "phase":
    outputDir = value("-o")
    os.makedirs(outputDir, exist_ok=True)
    with open(os.path.join(outputDir, "phasing.tsv"), "w") as fileOut:
        fileOut.write("#contig\tlength\tpartition\tstatus\taligned_bp\tunique_bp\treference_regions\treference_bp\treference_strand\n")
        for partition in range(1, int(value("--ploidy")) + 1):
            fileOut.write(f"c{partition}\t100\t{partition}\tassigned\t100\t100\tr:1-100\t100\t+\n")
        fileOut.write("u\t100\t.\tunaligned\t0\t0\t.\t.\t.\n")
    with open(os.path.join(outputDir, "coverage.tsv"), "w") as fileOut:
        fileOut.write("#reference\tlength\n")
elif mode == "homology_pair":
    partition = value("--partition")
    if partition == os.environ.get("STUB_FAIL_PARTITION"):
        sys.exit(1)
    with open(value("-o"), "w") as fileOut:
        fileOut.write(f">partition{partition}\n")
elif mode == "scaffold":
    for partition in assigned_partitions(value("-i")):
        with open(os.path.join(value("-o"), f"partition{partition}.fasta"), "w") as fileOut:
            fileOut.write(f">partition{partition}\n")
