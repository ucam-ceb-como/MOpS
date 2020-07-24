# *****************************************************************************
#
# File:                 mops.py
# Project:              python3
# Author(s):            Eric J. Bringley (eb656)
#
# Contact:
#   Mr. Eric J. Bringley
#   Department of Chemical Engineering
#   University of Cambridge
#   New Museums Site
#   Pembroke Street
#   Cambridge
#   CB2 3RA
#   United Kingdom
#
#   Email:   eb656@cam.ac.uk
#   Website: como.cheng.cam.ac.uk
#
# Purpose:
#   This python module provides two classes (case and exec) that make reading 
#   and working with kinetics input and output files a bit less cumbersome.
# *****************************************************************************

import sys
import pandas as pd
import numpy as np 
import re
import os
import xml.etree.ElementTree as ET
import subprocess

#define case class:
class case():
    #Initialization
    def __init__(self,path="WorkingDir"):
        # check that path is a string, else convert to stirng
        if(not isinstance(path, str)):
            pathstr=str(path)
            path=pathstr
        # check directory exists
        if(os.path.isdir(path)):
            # print("{} exists".format(path))
            self.directoryPath=path
            self.inputFile=os.path.join(self.directoryPath, "InputParams.xml")
        else:
            raise ValueError("{} is not a valid directory".format(path))
        self._caseName=path
        print("Reading in case from {}".format(self.directoryPath))
        self.files=[]
        self.fileDict={}
        # self._readInputFile()
        # self._readCatalogueFile()

    def readOutputFiles(self):
        "Reads the output files of kinetics"
        read=self._readOutputFiles()
        return read

    #return full case path:
    def caseName(self):
        "Returns the absolute case path"
        return self._caseName

    #return full case path:
    def path(self):
        "Returns the absolute case path"
        return os.path.abspath(self.directoryPath)

    #return file keys:
    def filenames(self):
        return self.fileDict.keys()

    #return InputParamsElementTree
    def inputfile(self):
        return self.inputParamsXML

    #Read Input file
    def _readInputFile(self):
        "reads the kinetics input file"
        return None
    
    #Read Catalogue file
    def _readOutputFiles(self):
        "reads the output files"
        catFileName = self.outputBaseName+"Catalogue.xml"
        catFile = ""

        if(os.path.isfile(catFile)):
            return True
        else:
            return False

            # raise ValueError("{} is not a valid file. Check that the simulation has completed running.".format(catFile))

    # #overload len function
    # def __len__(self):
    #     return(len(self.fileDict))

    # #overload [] operator
    # def __getitem__(self,key):
    #     return self.fileDict[key]

    #Destructor:
    def __del__(self):
        class_name = self.__class__.__name__
        print("{} {} deleted from python memory".format(class_name,self.directoryPath))

# define executable wrapper class
class exec():
    "A class that is a wrapper for a MOPS executable"


    #Initialization
    def __init__(self,executable,path="",args=""):
        self._execName=executable
        self._execPath=path
        self._argumentString=args
        self._args=self._argumentString.split()
        # check that path is a string, else convert to stirng
        if(not isinstance(self._execPath, str)):
            pathstr=str(path)
            path=pathstr
        try:
            completed=subprocess.run([self._execName], capture_output=True, text=True)
        except FileNotFoundError:
            if(os.path.isdir(self._execPath)):
                print("{} is not in your default path... trying again with full path".format(self._execName))
                try:
                    completed=subprocess.run([os.path.join(self._execPath,self._execName)],capture_output=True, text=True)
                    print("Found executable {}".format(self._execName))
                    self._execPath=os.path.abspath(self._execPath)
                except FileNotFoundError:
                    print("Cannot execute {}... please check the executable exists and the path is correct.".format(self._execName))
            else:
                print("{} is not in your default path... please provide a full path as an argument".format(self._execName))
                raise FileNotFoundError


    def run(self, caseobj):
        if(not isinstance(caseobj, case)):
            raise TypeError("Argument must be of type mops.case")
            completed=False
        else:
            #append filesep to casename
            rundir=caseobj.path()
            print("Running case: {}".format(rundir))
            command = ['cd', rundir, ";"] + [os.path.join(self._execPath,self._execName)] + self._args
            try:
                completed=subprocess.run(command, text=True,cwd=rundir)
            except FileNotFoundError:
                try:
                    command = [os.path.join(self._execPath,self._execName)] + self._args
                    completed=subprocess.run(command,text=True,cwd=rundir)
                except FileNotFoundError:
                    print("Cannot execute {}... please check the executable exists and the path is correct.".format(self._execName))
                    completed=False
        return completed

    def setArgs(self, argumentString):
        "Sets arguements for running MOpS executable (-p --flamepp -g 'gasphase.csv')"
        if(not isinstance(argumentString, str)):
            argumentString=str(argumentString)
        self._argumentString=argumentString
        self._args=self._argumentString.split()
        return None


    #Destructor:
    def __del__(self):
        class_name = self.__class__.__name__
        print("Wrapper for {} {} deleted from python memory".format(class_name,self._execName))
