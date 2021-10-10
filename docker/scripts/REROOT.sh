#!/bin/sh


for files in /srv/shiny-server/RpNRM/data/*.nwk
do

  cd /srv/shiny-server/RpNRM/data
  /usr/local/envs/hyphy/bin/hyphy  /srv/shiny-server/RpNRM/data/likeli.bf $files 
   




done


