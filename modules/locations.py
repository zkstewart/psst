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

# Location of job script templates within this repository
TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "templates")

# Program locations on the HPC
MMSEQS_DIR = "/home/stewarz2/various_programs/mmseqs/bin"
VARSCRIPT_DIR = "/home/stewarz2/scripts/Various_scripts"

# Database locations on the HPC
REFSEQ_PLANT_FASTA = "/home/stewarz2/various_programs/refseq_db/refseq_plant.protein.faa"
UNIREF90_FASTA = "/home/stewarz2/various_programs/uniref_db/uniref90.fasta"
IDMAPPING_FILE = "/home/stewarz2/various_programs/uniref_db/idmapping_selected.tab"
OBO_FILE = "/home/stewarz2/various_programs/uniref_db/go.obo"
