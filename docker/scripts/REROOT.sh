#!/bin/sh


for files in /srv/shiny-server/RpNRM/*.nwk
do



#    /home/gerrit/miniconda3/envs/py3/bin/HYPHYMPI likeli.bf $files
    #/usr/local/bin/hyphy likeli.bf $files 
     /usr/local/envs/hyphy  /srv/shiny-server/RpNRM/scripts/likeli.bf $files 
   




done


