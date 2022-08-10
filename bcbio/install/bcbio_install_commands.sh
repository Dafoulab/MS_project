wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py

python3 bcbio_nextgen_install.py [bcbio_path] --tooldir=[bcbio_tools_path] --nodata

export PATH="[bcbio_path]/anaconda/bin":"[bcbio_tools_path]/bin":$PATH

bcbio_nextgen.py upgrade -u skip \
--genomes hg38 \
--aligners bwa \
--aligners star \
--distribution ubuntu

bcbio_nextgen.py upgrade -u development --tools --data










