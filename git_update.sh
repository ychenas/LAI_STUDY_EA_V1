#!/bin/sh

#init git 
git init

#add files 
git add *.R *.rda spei* *.sh

#commit files 
git commit -m "update sh-files for version 1"

#set the origin, only for the first time
#git remote add origin git@github.com:ychenatsinca/LAI_STUDY_EA_V1.git
#add branch name, here is main 
#git branch -M main
#push commit files to the server/origin as master or branch 
git push -u origin main


