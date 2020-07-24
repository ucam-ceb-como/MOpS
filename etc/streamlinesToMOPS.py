# *****************************************************************************
#
# File:                 streamlinesForMOPS.py
# Project:              python3
# Author(s):            Eric J. Bringley (eb656)
#
# Contact:
#   Mr. Eric J. Bringley
#   Department of Chemical Engineering
#   University of Cambridge
#   West Cambridge Site
#   Philippa Fawcett Drive
#   Cambridge
#   CB3 0AS
#   United Kingdom
#
#   Email:   eb656@cam.ac.uk
#   Website: como.cheng.cam.ac.uk
#
# Purpose:
#   This script takes streamlines from OpenFOAM and formats them as inputs for
#   the TiO2 (TTIP) particle models in MOPS.
#
# *****************************************************************************

import sys
import pandas as pd
import numpy as np 
import re
import os
import xml.etree.ElementTree as ET
