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


# define case class:
class case():
    # Initialization
    def __init__(self, path="WorkingDir",
                 inputFile="mops.inx", sweepFile="sweep.xml",
                 gasphaseFile="gasphase.csv"
                 ):
        # Early declaration of object properties:
        self._sweepFilePath = ""
        self._sweepFile = None
        self._inputFilePath = ""
        self._inputFile = None
        self._caseName = ""
        self._outputFileName = ""
        self.files = []
        self.fileDict = {}
        self.timearray = []
        self.directoryPath = ""

        # check that path is a string, else convert to stirng
        if(not isinstance(path, str)):
            path = str(path)
        # check directory exists
        if(os.path.isdir(path)):
            self.directoryPath = path
        else:
            raise ValueError("{} is not a valid directory".format(path))
        self._inputFilePath = os.path.join(self.directoryPath, inputFile)
        # Read input file; returns true on success.
        if(not self._readInputFile()):
            raise ValueError("{} is not a file".format(self._inputFile))
        # Check for sweep file:
        self._sweepFilePath = os.path.join(self.directoryPath, sweepFile)
        if(not self._readSweepFile()):
            raise ValueError("{} is not a file".format(self._sweepFilePath))        
        self._caseName = path
        print("Reading in case from {}".format(self.directoryPath))
        # self._readOutputFiles()

    # read Input file
    def _readInputFile(self):
        """
            Reads the mops input file (*.inx)
            The inx file is an xml file. Use ElementTree module to read
            the input file and extract information
        """
        # Check that self._inputFile is a valid file
        if(not os.path.isfile(self._inputFilePath)):
            print("The input file {} does not exist".format(self._inputFile))
            return False
        # Try to read input file.
        self._inputFile = ET.parse(self._inputFilePath)
        self._inputFileRoot = self._inputFile.getroot()
        element_filename = self._inputFileRoot.find("./output/filename")
        self._outputFileName = element_filename.text
        # element_time_start = self._inputFile.find("./timeintervals/start")
        # self.timearray.append(element_time_start.text)
        times = self._inputFileRoot.findall("./timeintervals/time")
        for time in times:
            # psl files have 6 sig fig in their names, truncate to match:
            t = time.text
            tstr = (f"{float(t):.6}")
            self.timearray.append(tstr)
        return True

    # re-write Input file to disk
    def _writeInputFile(self, newTree=None, newName=""):
        if not newTree:
            newTree = self._inputFile
        if not newName:
            newName=self._inputFilePath
        newTree.write(newName)
        return None


    # read Sweep file
    def _readSweepFile(self):
        """
            Reads the particle input file (sweep.xml)
            Use ElementTree module to read the input file and extract information
        """
        # Check that self._inputFile is a valid file
        if(not os.path.isfile(self._sweepFilePath)):
            print("The input file {} does not exist".format(self._sweepFilePath))
            return False
        # Try to read input file.
        tree = ET.parse(self._sweepFilePath)
        self._sweepFile = tree.getroot()
        # Do nothing currently with sweep file
        return True

    def readOutputFiles(self):
        "Reads the output files of kinetics"
        # ensure latest input file is read:
        self._readInputFile()
        read = self._readOutputFiles()
        return read

    # return full case path:
    def caseName(self):
        "Returns the absolute case path"
        return self._caseName

    # return full case path:
    def path(self):
        "Returns the absolute case path"
        return os.path.abspath(self.directoryPath)

    # get M0 
    def getM0(self):
        "Gets the value of M0 in the input file"
        M0 = -1
        element_M0 = self._inputFileRoot.find("./maxm0")
        M0 = float(element_M0.text)
        return M0

    # Set M0 
    def setM0(self, newM0):
        "Sets the value of M0 in the input file"
        element_M0 = self._inputFileRoot.find("./maxm0")
        element_M0.text = str(newM0)
        return None       

    # Get pcount
    def getpcount(self):
        "Gets the value of pcount in the input file"
        pcount = -1
        element_pcount = self._inputFileRoot.find("./pcount")
        pcount = float(element_pcount.text)
        return pcount

    # Set pcount
    def setpcount(self, newpcount):
        "Sets the value of pcount in the input file"
        element_pcount = self._inputFileRoot.find("./pcount")
        element_pcount.text = str(newpcount)
        return None  

    # get filename_string 
    def getFilename(self):
        "Gets the value of output/filename"
        element = self._inputFileRoot.find("./output/filename")
        fname = element.text
        return fname

    # Set M0 
    def setFilename(self, newfname):
        "Sets the value of M0 in the input file"
        element = self._inputFileRoot.find("./output/filename")
        element.text = newfname
        return None  

    # Save input file
    def saveInputFile(self, newName=""):
        self._writeInputFile(newName=newName)
        return None

    # return file keys:
    def filenames(self):
        return self.fileDict.keys()

    # return time array
    def times(self):
        return self.timearray

    # return fine time
    def tfinal(self):
        return max(self.timearray)

    # return InputParamsElementTree
    def inputfile(self):
        return self._inputFileRoot

    # Read Catalogue file
    def _readOutputFiles(self):
        "reads the output files"
        # Read Particle Size List in timearray
        for time in self.timearray:
            # check if PSL file exists:
            pslFileName = "psl("+str(time)+"s)"
            self._readFile(pslFileName)
        # Read Aggregate Data files
        self._readFile("part")
        self._readFile("primary")
        self._readFile("primary-nodes")
        return True

    def _readFile(self, fileName):
        """
            Reads an output file given a generic file name
            Adds it to the fileDict
        """
        fullFileName = self._outputFileName+"-"+fileName+".csv"
        if(os.path.isfile(os.path.join(self.directoryPath, fullFileName))):
            # print("File {} found".format(fullFileName))
            key = fileName
            df = pd.read_csv(os.path.join(self.directoryPath, fullFileName))
            self.fileDict[key] = df
        else:
            print("File {} not found".format(fullFileName))
        return None

    # overload len function
    def __len__(self):
        return(len(self.fileDict))

    # overload [] operator
    def __getitem__(self, key):
        return self.fileDict[key]

    # Destructor:
    def __del__(self):
        class_name = self.__class__.__name__
        print("{} {} deleted from python memory".format(class_name,
                                                        self.directoryPath))


# define executable wrapper class
class exec():
    "A class that is a wrapper for a MOPS executable"

    # Initialization
    def __init__(self, executable, path="", args=""):
        self._execName = executable
        self._execPath = path
        self._argumentString = args
        self._args = self._argumentString.split()
        # check that path is a string, else convert to stirng
        if(not isinstance(self._execPath, str)):
            pathstr = str(path)
            path = pathstr
        try:
            completed = subprocess.run([self._execName], capture_output=True, text=True)
        except FileNotFoundError:
            if(os.path.isdir(self._execPath)):
                print("{} is not in your default path... trying again with full path".format(self._execName))
                try:
                    completed = subprocess.run([os.path.join(self._execPath, self._execName)], capture_output=True, text=True)
                    print("Found executable {}".format(self._execName))
                    self._execPath = os.path.abspath(self._execPath)
                except FileNotFoundError:
                    print("Cannot execute {}... please check the executable exists and the path is correct.".format(self._execName))
            else:
                print("{} is not in your default path... please provide a full path as an argument".format(self._execName))
                raise FileNotFoundError

    def run(self, caseobj):
        if(not isinstance(caseobj, case)):
            raise TypeError("Argument must be of type mops.case")
            completed = False
        else:
            # append filesep to casename
            rundir = caseobj.path()
            print("Running case: {}".format(rundir))
            command = [os.path.join(self._execPath, self._execName)] + self._args
            try:
                completed = subprocess.run(command, text=True,cwd=rundir)
            except FileNotFoundError:
                try:
                    command = [os.path.join(self._execPath, self._execName)] + self._args
                    completed = subprocess.run(command, text=True, cwd=rundir)
                except FileNotFoundError:
                    print("Cannot execute {}... please check the executable exists and the path is correct.".format(self._execName))
                    completed = False
        return completed

    def setArgs(self, argumentString):
        "Sets arguements for running MOpS executable (-p --flamepp -g 'gasphase.csv')"
        if(not isinstance(argumentString, str)):
            argumentString = str(argumentString)
        self._argumentString = argumentString
        self._args = self._argumentString.split()
        return None

    # Destructor:
    def __del__(self):
        class_name = self.__class__.__name__
        print("Wrapper for {} {} deleted from python memory".format(class_name, self._execName))
