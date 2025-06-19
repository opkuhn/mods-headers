#!/usr/bin/env python
# /lbt/mods_runtime/anaconda/bin/python
# /opt/anaconda/bin/python

###############################################################################
# Script to correct headers which have multiple keyword entries. All of
# the mods1r data taken on 20230526 had many duplicate keyword/value entries. 
# This was fixed by cycling the IC program on M1.RC. 
# But while the data were ingested into the archive, the incorrect headers
# led to their landing in the wrong UTdate subdirectories under Repository and,
# because TELRA and TELDEC were among the duplicate keywords, having incorrect 
# WCS information. Also - incorrect airmass information, objname, etc.
#
# This deletes the first instance of the keyword/value; reads the second value,
# which is the correct one, and then writes out the keyword/value/comment 
# string in the same place where it appeared originally.
#
# The index numbers were determined from the output of fitsverify.
#
# 20230527 opk/lbto
###############################################################################

import numpy as np
from astropy.io import fits
import argparse

parser = argparse.ArgumentParser(description = 'fix headers')
parser.add_argument('filename', type=str, help = 'file which needs fixing')
args = parser.parse_args()
filename= args.filename

# indices and keywords copied from output of fitsverify, but see below - astropy headers are 0-indexed. 

keyls = {120, 103, 109, 95, 100, 93, 99, 98, 94, 118, 106, 107, 87, 92, 85, 91, 90, 86, 25, 115, 108, 102, 113, 114, 111, 110, 105, 104, 119} 

keywd = {120 : ("AIRMASS" , "Airmass (secZD) at start of obs"), 
103 : ("DATE-OBS", " UTC Date at start of obs"), 
109 : ("EQUINOX", "Equinox of coordinates"), 
95 : ("GUIDEC", "Guide Star DEC"), 
100 : ("GUIEPOCH", "Guide Star Epoch"), 
93 : ("GUINAME", "Guide Star Name"), 
99 : ("GUIPMDEC", "Guide Star Dec proper motion [mas per yr]"), 
98 : ("GUIPMRA", "Guide Star RA ,proper motion [mas per yr]"), 
94 : ("GUIRA", "Guide Star RA"), 
118 : ("HA", "Hour Angle at start of obs" ), 
106 : ("LST-OBS", "Local Siderial Time at start of obs"), 
107 : ("MJD-OBS", "Modified JD=JD-2400000.5 at start of obs" ), 
87 : ("OBJDEC", "Target DEC" ), 
92 : ("OBJEPOCH", "Target Epoch" ), 
85 : ("OBJNAME", "Target Name" ), 
91 : ("OBJPMDEC", "Target Dec proper motion [mas per yr]" ), 
90 : ("OBJPMRA", "Target RA proper motion [mas per yr]" ), 
86 : ("OBJRA", "Target RA" ), 
25 : ("OBSERVAT", "Observatory Site"), 
115 : ("PARANGLE", " Parallactic Angle at start of obs [deg]"), 
108 : ("RADECSYS", "Coordinate System" ), 
102 : ("TCSLINK", "TCS Communications Link Status" ), 
113 : ("TELALT", "Telescope Altitude at start of obs [deg]" ), 
114 : ("TELAZ", "Telescope Azimuth at start of obs [deg]" ),
111 : ("TELDEC", "Telescope DEC" ), 
110 : ("TELRA", "Telescope RA" ), 
105 : ("TIMESYS", "Time System" ), 
104 : ("UTC-OBS", "UTC Time at start of obs" ), 
119 : ("ZD", "Zenith Distance at start of obs [deg]" )} 

f = fits.open(filename,mode='update')
h = f[0].header

for i in keyls:

      # delete all of the duplicate keywords, python uses 0-indexed
      del h[i-1]
      fitskey = keywd[i][0]
      val = h[fitskey]
      print("Deleted the first occurence of %s, at card %d. Now %s has correct value %s" % (keywd[i][0],i-1,keywd[i][0],val))
      print("Deleting the second occurence of %s" % (keywd[i][0]))
      del h[keywd[i][0]]
      print("Inserting the keyword/value/comment %s/%s/%s at position %d" % (keywd[i][0],val,keywd[i][1],i-1))
      h.insert(i-1,(keywd[i][0],val,keywd[i][1]))

# write this back to the original file  
print("Updating the file.")
f.flush()
