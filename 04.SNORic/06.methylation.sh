#########################################################################
# File Name: 06.methylation.sh
# Author: C.J. Liu
# Mail: chunjie.sam.liu@gmail.com
# Created Time: Tue 14 Mar 2017 12:40:17 AM CDT
#########################################################################
#!/bin/bash

Rscript 06.methylation.r

awk '{print NR"\t"$0}' /home/cliu18/liucj/projects/4.SNORic/database/methylation/methylation.tsv >/home/cliu18/liucj/projects/4.SNORic/database/methylation/methylation.tsv.n

sshpass -p 'u201012670' rsync -rvz  -e 'ssh -p 3000' --progress /home/cliu18/liucj/projects/4.SNORic/database/methylation/methylation.tsv.n liucj@211.69.207.247:/home/liucj/tmp/snoric/