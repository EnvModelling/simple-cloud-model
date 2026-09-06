# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 23:20:02 2020

@author: mccikpc2
"""

import os
import tempfile
from subprocess import check_output

import getpass
username=getpass.getuser()

def do_model_run():
    from runsDefine import outputDir

    inputFile=os.getcwd()+'/../namelist.pamm'
    
    dumpFileObj=tempfile.NamedTemporaryFile(delete=False)
    dumpFile=dumpFileObj.name
    
    tmpFileObj=tempfile.NamedTemporaryFile(delete=False)
    tmpFile=tmpFileObj.name
    
    print(tmpFile)
    print(dumpFile)
    
    changeFile(inputFile,dumpFile,'/tmp/output.nc','/tmp/' + username + '/output.nc')
    
    if not os.path.exists('/tmp/' + username):
        os.mkdir('/tmp/' + username)
    
    changeFile(dumpFile,tmpFile,\
            'outputfile = \'' + '/tmp/' + username + '/output.nc\'',\
            'outputfile = \'' + outputDir + '/' + username + \
             '/output.nc\'')
    
    str1='./main.exe ' + tmpFile
    
    result = check_output(str1, shell=True, cwd='../').decode()

    #print(result)

"""
    0. function to change file
"""
def changeFile(inFile,outFile,inString,outString):
    fin = open(inFile,"rt")
    
    lines=[]
    for line in fin:
        lines.append(line)
        
    fin.close()


    fout = open(outFile,"wt")

    for line in lines:
        fout.write(line.replace(inString,outString))
    

    fout.close()
    
if __name__=="__main__":
    do_model_run()