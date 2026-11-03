# Simulation of small genomes, variants and reads for psst's integration tests. All
# simulation is seeded so that test data are identical between runs.

import gzip, os, random, subprocess

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")

def other_base(rng, base, exclude=()):
    return rng.choice([ b for b in BASES if b != base and b not in exclude ])

def revcomp(seq):
    return seq[::-1].translate(COMPLEMENT)

def write_fasta(fileName, seqs):
    '''
    Parameters:
        fileName -- a string of the FASTA file to write.
        seqs -- a dict of sequence IDs to sequence strings.
    '''
    with open(fileName, "w") as fileOut:
        for seqID, seq in seqs.items():
            fileOut.write(f">{seqID}\n")
            for i in range(0, len(seq), 60):
                fileOut.write(seq[i:i+60] + "\n")

def random_genome(rng, sizes, homopolymers=False):
    '''
    Parameters:
        rng -- a random.Random instance.
        sizes -- a dict of contig IDs to lengths.
        homopolymers -- OPTIONAL; a boolean indicating whether to plant 8 bp A runs
                        every 20 kb in contigs of at least 40 kb, so that indels within
                        them test left-alignment.
    Returns:
        genome -- a dict of contig IDs to lists of bases.
    '''
    genome = { c: [ rng.choice(BASES) for _ in range(n) ] for c, n in sizes.items() }
    if homopolymers:
        for c, n in sizes.items():
            if n < 40000:
                continue
            for p in range(5000, n - 5000, 20000):
                genome[c][p-1:p+9] = list("CAAAAAAAAG")
    return genome

def apply_variants(genome, variants):
    '''
    Parameters:
        genome -- a dict of contig IDs to lists of bases.
        variants -- a list of (contig, 0-based position, ref length, alt string) tuples
                    which must not overlap.
    Returns:
        haplotype -- a dict of contig IDs to sequence strings.
    '''
    haplotype = {}
    for c, seq in genome.items():
        seq = seq[:]
        for _, p, refLength, alt in sorted([ v for v in variants if v[0] == c ], key=lambda v: -v[1]):
            seq[p:p+refLength] = list(alt)
        haplotype[c] = "".join(seq)
    return haplotype

def simulate_population(rng, genome, samples):
    '''
    Simulates diploid samples carrying SNPs, multiallelic SNPs, MNPs, deletions and
    insertions (including single base deletions within homopolymer runs).

    Parameters:
        rng -- a random.Random instance.
        genome -- a dict of contig IDs to lists of bases.
        samples -- a list of sample names.
    Returns:
        haplotypes -- a dict of sample names to a pair of haplotypes (dicts of contig
                      IDs to sequence strings).
    '''
    edits = { (s, h): [] for s in samples for h in (0, 1) }
    for c, seq in genome.items():
        for p in range(1000, len(seq) - 1000, 1500):
            r = seq[p]
            k = rng.random()
            if k < 0.6: # SNP with random genotypes
                alt = other_base(rng, r)
                for s in samples:
                    gt = rng.choice([(0, 0), (0, 1), (1, 1), (0, 1)])
                    for h in (0, 1):
                        if gt[h]:
                            edits[(s, h)].append((c, p, 1, alt))
            elif k < 0.7: # multiallelic SNP
                a1 = other_base(rng, r)
                a2 = other_base(rng, r, (a1,))
                for i, s in enumerate(samples):
                    edits[(s, i % 2)].append((c, p, 1, a1 if i % 2 == 0 else a2))
            elif k < 0.8: # MNP (two adjacent substitutions in phase)
                alt = other_base(rng, r) + other_base(rng, seq[p+1])
                for s in samples[:2]:
                    edits[(s, 0)].append((c, p, 2, alt))
            elif k < 0.9: # 3 bp deletion
                for s in samples:
                    if rng.random() < 0.7:
                        edits[(s, rng.randint(0, 1))].append((c, p, 4, r))
            else: # 4 bp insertion
                ins = "".join(rng.choice(BASES) for _ in range(4))
                for s in samples:
                    if rng.random() < 0.7:
                        edits[(s, rng.randint(0, 1))].append((c, p, 1, r + ins))
        for p in range(5000, len(seq) - 5000, 20000): # deletions within homopolymers
            if seq[p:p+8] == list("AAAAAAAA"):
                for s in samples[:2]:
                    edits[(s, 1)].append((c, p + 6, 1, ""))
    return { s: (apply_variants(genome, edits[(s, 0)]), apply_variants(genome, edits[(s, 1)]))
             for s in samples }

def add_errors(rng, seq, subRate, indelRate):
    out = []
    for b in seq:
        x = rng.random()
        if x < indelRate / 2: # deletion
            continue
        if x < indelRate: # insertion
            out.append(b)
            out.append(rng.choice(BASES))
            continue
        out.append(other_base(rng, b) if x < indelRate + subRate else b)
    return "".join(out)

def simulate_paired_reads(rng, haplotypes, prefix, coverage, readLength=150, insert=400,
                          hic=False, duplicateFraction=0.0):
    '''
    Writes paired reads from a set of haplotypes as <prefix>_1.fq.gz and <prefix>_2.fq.gz.

    Parameters:
        rng -- a random.Random instance.
        haplotypes -- a list of haplotypes (dicts of contig IDs to sequence strings).
        prefix -- a string of the output file prefix.
        coverage -- a float of the read coverage per haplotype.
        hic -- OPTIONAL; a boolean indicating whether to pair each read with a mate from
               a random distant location, as in Hi-C.
        duplicateFraction -- OPTIONAL; a float of the fraction of pairs to duplicate.
    '''
    pairs = []
    for haplotype in haplotypes:
        contigs = list(haplotype)
        weights = [ len(haplotype[c]) for c in contigs ]
        numPairs = int(sum(weights) * coverage / (2 * readLength))
        for _ in range(numPairs):
            c = rng.choices(contigs, weights)[0]
            seq = haplotype[c]
            start = rng.randint(0, len(seq) - insert)
            r1 = seq[start:start+readLength]
            if hic:
                c2 = rng.choices(contigs, weights)[0]
                start2 = rng.randint(0, len(haplotype[c2]) - readLength)
                r2 = revcomp(haplotype[c2][start2:start2+readLength])
            else:
                r2 = revcomp(seq[start+insert-readLength:start+insert])
            if rng.random() < 0.5:
                r1, r2 = r2, r1
            pairs.append((add_errors(rng, r1, 0.002, 0), add_errors(rng, r2, 0.002, 0)))
    pairs += [ pairs[i] for i in range(int(len(pairs) * duplicateFraction)) ]
    with gzip.open(f"{prefix}_1.fq.gz", "wt") as out1, gzip.open(f"{prefix}_2.fq.gz", "wt") as out2:
        for i, (r1, r2) in enumerate(pairs):
            out1.write(f"@p{i}\n{r1}\n+\n{'I' * len(r1)}\n")
            out2.write(f"@p{i}\n{r2}\n+\n{'I' * len(r2)}\n")

def simulate_long_reads(rng, haplotypes, fileName, coverage):
    '''
    Writes nanopore-like reads (5-15 kb with ~2% errors, from both strands) from a set
    of haplotypes to a gzipped FASTQ file.
    '''
    with gzip.open(fileName, "wt") as fileOut:
        i = 0
        for haplotype in haplotypes:
            contigs = list(haplotype)
            weights = [ len(haplotype[c]) for c in contigs ]
            for _ in range(int(sum(weights) * coverage / 10000)):
                seq = haplotype[rng.choices(contigs, weights)[0]]
                length = min(len(seq), rng.randint(5000, 15000))
                start = rng.randint(0, len(seq) - length)
                read = seq[start:start+length]
                if rng.random() < 0.5:
                    read = revcomp(read)
                read = add_errors(rng, read, 0.004, 0.016)
                i += 1
                fileOut.write(f"@r{i}\n{read}\n+\n{'5' * len(read)}\n")

def simulate_pangenome(rng, genome):
    '''
    Simulates the variants of a pangenome with two diploid samples (P1 and P2) relative
    to a reference, including multiallelic SNPs, MNPs, complex alleles, and SNPs within
    another haplotype's deletion.

    Returns:
        records -- a list of (contig, 0-based position, REF, [ALTs], {haplotype: allele})
                   tuples, where haplotypes 0 and 1 belong to P1, and 2 and 3 to P2.
    '''
    records = []
    for c, seq in genome.items():
        p = 500
        while p < len(seq) - 500:
            r = seq[p]
            k = rng.random()
            haps = set(rng.sample(range(4), rng.choice([1, 1, 2, 2, 3])))
            carriers = { h: 1 for h in haps }
            if k < 0.45: # SNP
                records.append((c, p, r, [other_base(rng, r)], carriers))
            elif k < 0.55: # multiallelic SNP
                a1 = other_base(rng, r)
                records.append((c, p, r, [a1, other_base(rng, r, (a1,))], { 0: 1, 2: 2 }))
            elif k < 0.65: # MNP
                ref = "".join(seq[p:p+3])
                records.append((c, p, ref, [other_base(rng, ref[0]) + ref[1] + other_base(rng, ref[2])], carriers))
            elif k < 0.75: # deletion
                ref = "".join(seq[p:p+rng.randint(2, 15)+1])
                records.append((c, p, ref, [ref[0]], carriers))
            elif k < 0.85: # insertion
                records.append((c, p, r, [r + "".join(rng.choice(BASES) for _ in range(rng.randint(1, 12)))], carriers))
            elif k < 0.92: # complex
                ref = "".join(seq[p:p+5])
                records.append((c, p, ref, [other_base(rng, ref[0]) + "".join(rng.choice(BASES) for _ in range(2))], carriers))
            else: # a deletion on one haplotype containing a SNP on another
                ref = "".join(seq[p:p+8])
                records.append((c, p, ref, [ref[0]], { 0: 1 }))
                records.append((c, p + 4, seq[p+4], [other_base(rng, seq[p+4])], { 2: 1 }))
            p += rng.randint(250, 450)
    return records

def write_pangenome_vcf(fileName, genome, records):
    '''
    Writes pangenome records as a bgzipped, indexed VCF like cactus-pangenome's (i.e.,
    contigs named without the reference sample prefix, and phased P1 and P2 genotypes).
    '''
    plain = fileName[:-3] if fileName.endswith(".gz") else fileName
    with open(plain, "w") as fileOut:
        fileOut.write("##fileformat=VCFv4.2\n")
        for c, seq in genome.items():
            fileOut.write(f"##contig=<ID={c},length={len(seq)}>\n")
        fileOut.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fileOut.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP1\tP2\n")
        for c, p, ref, alts, carriers in records:
            gt = [ str(carriers.get(h, 0)) for h in range(4) ]
            fileOut.write(f"{c}\t{p+1}\t.\t{ref}\t{','.join(alts)}\t60\tPASS\t.\tGT\t{gt[0]}|{gt[1]}\t{gt[2]}|{gt[3]}\n")
    subprocess.run(["bgzip", "-f", plain], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", fileName], check=True)

def pangenome_haplotypes(genome, records):
    '''
    Returns:
        haplotypes -- a list of the four pangenome haplotypes (P1 then P2), each a dict
                      of contig IDs to sequence strings.
    '''
    return [ apply_variants(genome, [ (c, p, len(ref), alts[carriers[h]-1])
                                      for c, p, ref, alts, carriers in records if carriers.get(h, 0) ])
             for h in range(4) ]

def build_pangenome_gbz(workDir, referenceFasta, vcfFile, gbzFile, refSample="reference"):
    '''
    Builds a GBZ graph like cactus-pangenome's from a reference and phased VCF: the VCF's
    samples become haplotype paths, and the reference becomes reference paths of sample
    refSample (e.g., 'reference#0#chr1').
    '''
    vg = lambda *args, **kwargs: subprocess.run(["vg", *args], check=True, cwd=workDir, **kwargs)
    with open(os.path.join(workDir, "graph.vg"), "w") as fileOut:
        vg("construct", "-r", referenceFasta, "-v", vcfFile, "-a", stdout=fileOut, stderr=subprocess.DEVNULL)
    vg("gbwt", "-x", "graph.vg", "-v", vcfFile, "-o", "haplotypes.gbwt", stderr=subprocess.DEVNULL)
    vg("gbwt", "-x", "graph.vg", "-E", "-o", "references.gbwt", stderr=subprocess.DEVNULL)
    vg("gbwt", "-m", "-o", "merged.gbwt", "haplotypes.gbwt", "references.gbwt", stderr=subprocess.DEVNULL)
    vg("gbwt", "-x", "graph.vg", "merged.gbwt", "--gbz-format", "-g", "graph.gbz", stderr=subprocess.DEVNULL)
    gfa = subprocess.run(["vg", "convert", "-f", "graph.gbz"], check=True, cwd=workDir,
                         capture_output=True, text=True).stdout

    # Rewrite the reference paths as walks of the reference sample
    segmentLengths = {}
    for line in gfa.split("\n"):
        if line.startswith("S\t"):
            fields = line.split("\t")
            segmentLengths[fields[1]] = len(fields[2])
    out = []
    for line in gfa.split("\n"):
        if line == "" or (line.startswith("W\t") and "_gbwt_ref" in line):
            continue
        if line.startswith("H"):
            out.append(f"H\tVN:Z:1.1\tRS:Z:{refSample}")
        elif line.startswith("P\t"):
            fields = line.split("\t")
            steps = fields[2].split(",")
            walk = "".join((">" if s[-1] == "+" else "<") + s[:-1] for s in steps)
            length = sum(segmentLengths[s[:-1]] for s in steps)
            out.append(f"W\t{refSample}\t0\t{fields[1]}\t0\t{length}\t{walk}")
        else:
            out.append(line)
    with open(os.path.join(workDir, "pangenome.gfa"), "w") as fileOut:
        fileOut.write("\n".join(out) + "\n")
    vg("gbwt", "-G", "pangenome.gfa", "--gbz-format", "-g", gbzFile, stderr=subprocess.DEVNULL)
