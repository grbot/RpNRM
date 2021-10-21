#!/bin/sh

rm data.csv
rm data1.nwk
cd /srv/shiny-server/rpnrm
for files in /srv/shiny-server/rpnrm/*.nwk
do

  #cd /srv/shiny-server/rpnrm
  /usr/local/envs/hyphy/bin/hyphy  /srv/shiny-server/rpnrm/likeli.bf $files 
   




done


