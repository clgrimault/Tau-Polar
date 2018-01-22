#!/usr/bin/python
import os
import sys
import re


print " =========  Starting submission on CRAB ========"
filelist=[]
for file in os.listdir("."):
    if sys.argv[1] in file and file.endswith("root"):
        filelist.append(file)
        


for file in filelist:
    os.system('hadd Combined.root file')
#print filelist
        
