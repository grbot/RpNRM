#!/bin/sh


for files in *.nwk
do



#    /home/gerrit/miniconda3/envs/py3/bin/HYPHYMPI likeli.bf $files
   /usr/local/bin/hyphy likeli.bf $files 




done


