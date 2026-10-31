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

import os, re

# Separator (and optional 'R') left at the end of a paired read prefix once the 1/2 is
# removed (e.g., 'sample_R' from 'sample_R1.fq.gz'); job templates which discover read
# files must strip this with the same regex (i.e., sed -E 's/[._-]R?$//')
PAIR_SEPARATOR_REGEX = re.compile(r"[._-]R?$")

def get_read_files(readsDir, readsSuffix, singleEnd):
    '''
    Locates read files within a directory and checks that paired reads are complete.
    Paired read files must have a 1 or 2 immediately preceding readsSuffix (e.g.,
    'sample_1.fq.gz' and 'sample_2.fq.gz' with readsSuffix '.fq.gz').

    Parameters:
        readsDir -- a string of the directory containing read files.
        readsSuffix -- a string of the suffix which identifies read files.
        singleEnd -- a boolean indicating whether reads are single-end rather than paired.
    Returns:
        samples -- a list of tuples like (sampleName, forwardFile, reverseFile), sorted
                   by forward file name; reverseFile is None if singleEnd is True.
    '''
    readFiles = sorted([ f for f in os.listdir(readsDir) if f.endswith(readsSuffix) ])
    if readFiles == []:
        raise FileNotFoundError(f"No files ending with '{readsSuffix}' were found in {readsDir}")

    samples = []
    if singleEnd:
        for file in readFiles:
            samples.append((file[:-len(readsSuffix)], os.path.join(readsDir, file), None))
    else:
        forwardReads = [ f for f in readFiles if f.endswith(f"1{readsSuffix}") ]
        reverseReads = set([ f for f in readFiles if f.endswith(f"2{readsSuffix}") ])
        for file in readFiles:
            if file not in forwardReads and file not in reverseReads:
                raise ValueError(f"{file} ends with the expected suffix '{readsSuffix}' but is not preceded by a 1 or 2!")

        # Pair each forward read file with its reverse read file
        for file in forwardReads:
            prefix = file[:-len(f"1{readsSuffix}")]
            pairFile = f"{prefix}2{readsSuffix}"
            if pairFile not in reverseReads:
                raise FileNotFoundError(f"{file} has no matching reverse read file ({pairFile}) in {readsDir}")
            reverseReads.remove(pairFile)
            samples.append((PAIR_SEPARATOR_REGEX.sub("", prefix), os.path.join(readsDir, file),
                            os.path.join(readsDir, pairFile)))
        if reverseReads != set():
            raise FileNotFoundError("Reverse read files have no matching forward read file: " +
                                    ", ".join(sorted(reverseReads)))

    # Validate sample names, which name output files
    sampleNames = [ s[0] for s in samples ]
    for sampleName in sampleNames:
        if sampleName == "":
            raise ValueError(f"A read file in {readsDir} has no sample name before '{readsSuffix}'")
        if sampleNames.count(sampleName) > 1:
            raise ValueError(f"Sample name '{sampleName}' is shared by multiple read files in {readsDir}")
    return samples
