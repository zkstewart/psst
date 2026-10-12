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

import os

from .config import CONFIG

# Location of job script templates within this repository
TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "templates")

# Locations below are read from my_config.json or default_config.json (see config.py)

# Program locations on the HPC
MMSEQS_DIR = CONFIG["programs"]["mmseqs_dir"]
VARSCRIPT_DIR = CONFIG["programs"]["various_scripts_dir"]
GENSCRIPT_DIR = CONFIG["programs"]["genome_analysis_scripts_dir"]
COMPLEASM_DIR = CONFIG["programs"]["compleasm_dir"]
DORADO_DIR = CONFIG["programs"]["dorado_dir"]

# Database locations on the HPC
REFSEQ_PLANT_FASTA = CONFIG["databases"]["refseq_plant_fasta"]
UNIREF90_FASTA = CONFIG["databases"]["uniref90_fasta"]
IDMAPPING_FILE = CONFIG["databases"]["uniref_idmapping_file"]
OBO_FILE = CONFIG["databases"]["go_obo_file"]

# Reference genome locations on the HPC
REFERENCE_HAP1_FASTA = CONFIG["references"]["reference_hap1_fasta"]
REFERENCE_HAP2_FASTA = CONFIG["references"]["reference_hap2_fasta"]
