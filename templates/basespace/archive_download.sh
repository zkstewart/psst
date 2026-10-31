####

# Specify the directory to archive
READSDIR=@@readsDir@@

# Specify the prefix of the archive file
PREFIX=@@prefix@@

####

tar --create --dereference --file=${PREFIX}.tar ${READSDIR}
md5sum ${PREFIX}.tar > ${PREFIX}.tar.md5
