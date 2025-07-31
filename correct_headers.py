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

keyls = {119, 102, 108, 92, 93, 94, 97, 98, 99, 117, 105, 106, 86, 91, 84, 90, 89, 85, 24, 114, 107, 101, 112, 113, 110, 109, 104, 103, 118} 

keywd = {119 : ("AIRMASS" , "Airmass (secZD) at start of obs"), 
102 : ("DATE-OBS", " UTC Date at start of obs"), 
108 : ("EQUINOX", "Equinox of coordinates"), 
92 : ("GUINAME", "Guide Star Name"), 
93 : ("GUIRA", "Guide Star RA"), 
94 : ("GUIDEC", "Guide Star DEC"), 
97 : ("GUIPMRA", "Guide Star RA ,proper motion [mas per yr]"), 
98 : ("GUIPMDEC", "Guide Star Dec proper motion [mas per yr]"), 
117 : ("HA", "Hour Angle at start of obs" ), 
99 : ("GUIEPOCH", "Guide Star Epoch"), 
105 : ("LST-OBS", "Local Siderial Time at start of obs"), 
106 : ("MJD-OBS", "Modified JD=JD-2400000.5 at start of obs" ), 
86 : ("OBJDEC", "Target DEC" ), 
91 : ("OBJEPOCH", "Target Epoch" ), 
84 : ("OBJNAME", "Target Name" ), 
90 : ("OBJPMDEC", "Target Dec proper motion [mas per yr]" ), 
89 : ("OBJPMRA", "Target RA proper motion [mas per yr]" ), 
85 : ("OBJRA", "Target RA" ), 
24 : ("OBSERVAT", "Observatory Site"), 
114 : ("PARANGLE", " Parallactic Angle at start of obs [deg]"), 
107 : ("RADECSYS", "Coordinate System" ), 
101 : ("TCSLINK", "TCS Communications Link Status" ), 
112 : ("TELALT", "Telescope Altitude at start of obs [deg]" ), 
113 : ("TELAZ", "Telescope Azimuth at start of obs [deg]" ),
110 : ("TELDEC", "Telescope DEC" ), 
109 : ("TELRA", "Telescope RA" ), 
104 : ("TIMESYS", "Time System" ), 
103 : ("UTC-OBS", "UTC Time at start of obs" ), 
118 : ("ZD", "Zenith Distance at start of obs [deg]" )} 

f = fits.open(filename,mode='update')
h = f[0].header

for i in keyls:

      # delete all of the duplicate keywords
      del h[i]
      fitskey = keywd[i][0]
      val = h[fitskey]
      print("Deleted the first occurrence of %s, at card \#%d. Now %s has correct value %s" % (keywd[i][0],i,keywd[i][0],val))
      print("Deleting the second occurrence of %s" % (keywd[i][0]))
      del h[keywd[i][0]]
      print("Inserting the keyword/value/comment %s/%s/%s at position %d" % (keywd[i][0],val,keywd[i][1],i))
      h.insert(i,(keywd[i][0],val,keywd[i][1]))

# write this back to the original file  
print("Updating the file.")
f.flush()
