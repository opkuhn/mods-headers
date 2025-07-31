#!/usr/bin/env python

"""

tcs2hdr.py -

2025-Jun-13 1st version, opk/lbto
2025-Jun-17   - make a list of duplicates, use it to identify duplicates since when TCS info is missing, there may be added entries that
                 are not duplicates but also do not have a comment.
2025-Jul-19 - generalized to both blue and red channels of both MODS. (It had been just MODS2)
2025-Jul-29 - caught and corrected a corner case where some DATE-OBS and MJD-OBS cards had the equal sign but not value. 
                     
This script fills in the values of missing TCS keywords, using the TCSTATUS information from the ISIS log. 

Headers of some MODS images, particularly MODS2R images, are often missing TCS keywords: DATE-OBS, TELRA, TELDEC, guider and temperature data, etc. 
This prevents the files from flowing into the archive.


In addition, when TCS info is missing, there are also duplicate keywords added onto the end of the header. The presence of these also 
prevents data flow to the archive, as these images do not pass fitsverify.

In this script, header cards with missing information are identified because they lack an "=". The keyword and comment exist, but the
"=" and the value are missing. 
After inspection of the ISIS logs, it appears that:

(1) missing TCS info is indeed recorded in the ISIS log, in the TCSTATUS line that corresponds to the start of the exposure.
(2) missing TCS info often follows an IE timeout
(3) In the set of duplicate keywords, the AGWX position is often 1, i.e. at a limit. This may be related to an issue with sometimes non-sensical 
  AGW positions --- unsure whether part of the cause or just the effect?


"""

import numpy as np
from astropy.io import fits
from astropy.time import Time
import argparse
import os
import subprocess
import shlex
import sys
from datetime import datetime,timedelta



def find_dups(outfile):
   hdrcom = ['COMMENT', 'HISTORY', 'END']
   with open(outfile,"r") as f:
         kwd = []
         duplist = []
         for line in f:
            newkey = line.split(" ")[0].strip("'=','\n'")
            if newkey in kwd and newkey not in hdrcom:
               duplist.append(newkey)
               print("duplist: %s" % ( newkey))
            kwd.append(newkey)
   print (duplist)
   return duplist

def is_duplicate(line):
  # duplicates are signaled by no comment after the "/" and not in simplekeys
  llen = len(line.split("/"))
  l0 = len(line.split("/")[0].split("="))
  if l0 == 2:
     k = (line.split("/")[0].split("=")[0]).strip()
     kv = line.split("/")[0].split("=")[1]
  if llen==2:
     kc = (line.strip()).split("/")[1]
  if llen==3:
     kc = (line.strip()).split("/")[2]
  if (llen==2 and kc=='' and k not in simplekeys) or (llen==3 and kc=='' and k not in simplekeys):
         return True

def is_number(istr):
   try:
      float(istr)
      return True
   except ValueError:
      return False

def int_or_float(num):
   if (int(float(num)) == float(num)):
      return int(float(num))
   else:
      return float(num)

def tcsstatus(inImage,isislog,chan,L,minst):
   G = open(isislog,"r")
   isis = G.readlines()
   ti = -1
   c = str.upper(chan)

   # find the line reporting that the input image was written
   # the missing TCS information can be found in the preceding TCSTATUS line, which should have a timestamp 
   # that is exptime before the image was written.
   #srchstr = "DONE: Wrote LASTFILE=/lhome/data/"+inImage  # for mods2
   #srchstr = "STATUS: Wrote LASTFILE=/lhome/data/"+inImage # for mods1
   srchstr = "Wrote LASTFILE=/lhome/data/"+inImage
   fl = [i for i in range(len(isis)) if srchstr in isis[i]]
   if len(fl) == 1: 
        fi = int(fl[0])
        print(fi, isis[fi])
        L.write("%d %s\n" % (fi,isis[fi]))
        if minst == "2" and c ==  "R":
               tall = [i for i in range(len(isis)) if "M2.TC>M2.RC DONE: TCSTATUS" in isis[i]]
               iall = [i for i in range(len(isis)) if "M2.IE>M2.RC DONE: ISTATUS" in isis[i]]
               eall = [i for i in range(len(isis)) if "M2.RC>MC2 WARNING: GO IE Timeout with ERROR" in isis[i]]
        if minst == "2" and c ==  "B":
               tall = [i for i in range(len(isis)) if "M2.TC>M2.BC DONE: TCSTATUS" in isis[i]]
               iall = [i for i in range(len(isis)) if "M2.IE>M2.BC DONE: ISTATUS" in isis[i]]
               eall = [i for i in range(len(isis)) if "M2.BC>MC2 WARNING: GO IE Timeout with ERROR" in isis[i]]
        if minst == "1" and c ==  "R":
               tall = [i for i in range(len(isis)) if "TC>M1.RC DONE: TCSTATUS" in isis[i]]
               iall = [i for i in range(len(isis)) if "M1.IE>M1.RC DONE: ISTATUS" in isis[i]]
               eall = [i for i in range(len(isis)) if "M1.RC>MC1 WARNING: GO IE Timeout with ERROR" in isis[i]]
        if minst == "1" and c ==  "B":
               tall = [i for i in range(len(isis)) if "TC>M1.BC DONE: TCSTATUS" in isis[i]]
               iall = [i for i in range(len(isis)) if "M1.IE>M1.BC DONE: ISTATUS" in isis[i]]
               eall = [i for i in range(len(isis)) if "M1.BC>MC1 WARNING: GO IE Timeout with ERROR" in isis[i]]

        print (tall)
        # if there is no tcstatus line that corresponds to the start of exposure
        # in this case, len(fl) ==1 but ti = -1
        if len(tall) > 0:
           ti = max(filter(lambda nnn: nnn<fi, tall),default=0)
           if ti > 0:
               print(ti, isis[ti])
               L.write("%d %s\n" % (ti,isis[ti]))
           else:
               print(isis[fi].strip(" ")[0])
               utend = isis[fi].split(" ")[0]
               ti = -1 * Time(utend,format='isot',scale='utc').mjd
               L.write("no TCSTATUS: %s \n" % (isis[fi]))
        else: 
           print(isis[fi].strip(" ")[0])
           utend = isis[fi].split(" ")[0]
           ti = -1 * Time(utend,format='isot',scale='utc').mjd
           L.write("no TCSTATUS: %s \n" % (isis[fi]))
        if len(iall) > 0:
           ii = max(filter(lambda nnn: nnn<fi, iall))
           print(ii, isis[ii])
           L.write("%d %s\n" % (ii,isis[ii]))
        if len(eall) > 0:
           ei = max(filter(lambda nnn: nnn<fi, eall))
           print(ei, isis[ei])
           L.write("%d %s\n" % (ei,isis[ei]))
   elif len(fl) != 0: 
        print("Error: len(fl)=",len(fl),"****",fl)
        L.write("Error: len(fl)=",len(fl),"****\n",fl)

   G.close()
   if len(fl) == 1 and ti>0: 
        return isis[ti]
   elif len(fl) == 1 and ti<=0: 
        return ti
   elif len(fl) != 1: 
        return len(fl)

def readhead(inImage):
   # This is necessary for recognizing duplicate keywords 

   outfile = "header.dat"
   text_file = open(outfile,"w")

   cmd = "imhead " + str(inImage)

   try:
      process_output = subprocess.run(shlex.split(cmd),text=True,capture_output=True)
   except subprocess.CalledProcessError:
      print ("Error detected while executing command %s" % (cmd))

   text_file.write("%s" % process_output.stdout)
   text_file.close()
   return outfile
   
# the MODS header keywords that have a value but do not have a comment 
simplekeys = ['SIMPLE', 'BITPIX',  'NAXIS', 'NAXIS1',  'NAXIS2', 'BZERO', 'OVERSCNX', 'OVERSCNY', 'CCDROI']
hdrcom = ['COMMENT', 'HISTORY', 'END']

#isisdir = "../ISIS/"
note = "***"

parser = argparse.ArgumentParser(description = 'missing TCS info and duplicate header entries')

group = parser.add_mutually_exclusive_group()
group.add_argument('--image', type=str, dest='inImage', help = 'input image')
group.add_argument('--list', type=str, dest='inList', help = 'input list of images')

args = parser.parse_args()

if args.inImage:
   inImage = args.inImage
   nstr = len(str.split(inImage,"."))
   nimgs = 1

if args.inList:
   inList = args.inList
   print ("input image list is %s" % (inList))
   nstr = len(str.split(inList,"."))

   imglist = np.genfromtxt(inList,dtype=str)
   nimgs = len(imglist)


# Open a log file to make notes. 

logfile = datetime.strftime(datetime.now(),"%Y-%m-%dT%H:%M:%S") + '_addtcs2hdr.log'
L = open(logfile,"w")

L.write(" ------- tcs2hdr.py logfile  -------         \n\n") 

#### loop over the list of, or just the single, image(s)

for e in range(nimgs):

   if nimgs>1:
      inImage = imglist[e]

   L.write("#################################################\n\n")
   L.write("Working on %s ... \n\n" % (inImage))

   ####### ISIS logs ######
   # use the image filename to determine which ISIS log we need to open; relies on the date in the filename
   # being correct, which isn't always the case as the observer sometimes forgets to refresh it. 
   yymmdd = inImage.split(".")[1]
   instr = inImage.split(".")[0]
   chan = instr[5] # 'b' or 'r'
   minst = instr[4] # '1' or '2'

   fdate = datetime.strptime(yymmdd,"%Y%m%d")
   day = timedelta(days=1)
   fdate = [datetime.strftime((fdate-day),"%Y%m%d"), datetime.strftime((fdate+day),"%Y%m%d")]

   if minst == '1':
       isisdir = "/home/lbto/scratch/okuhn/dataCleanup_20250718/mods1data/ISIS/ISIS/"
   elif minst == '2':
       isisdir = "/home/lbto/scratch/okuhn/dataCleanup_20250718/mods2data/ISIS/ISIS/"

   isislog = isisdir + "isis." + str(yymmdd) + ".log"

   # if tisis == 0, then search in the ISIS log from the day before or after
   ntries = 0
   while ((tisis := tcsstatus(inImage,isislog,chan,L,minst)) == 0 and ntries<2):
        print(isislog)
        isislog = isisdir + "isis." + str(fdate[ntries]) + ".log"
        ntries = ntries+1
   if ntries == 2:
        print ("Error: Reached max number of tries")
   if is_number(tisis) and tisis < 0:
        print ("No TCS status line from which to draw missing TCS keyword values.")

   #convert tisis into a list of arrays 
   L.write(("Using values from %s\n" % (isislog)))

   # use shlex to split on spaces _except_ when these are enclosed by quotes
   if not is_number(tisis): 
        ilist = (shlex.split(tisis.strip()))
        ia=[]
        for i in range(len(ilist)):
           ik = ilist[i].split("=")[0].strip()
           if len(ilist[i].split("="))>1:
              #iv = (ilist[i].split("=")[1].strip())
              # remove added quotes around objname and guiname
              iv = ilist[i].split("=")[1].strip("'")
           else: 
              iv = 0
           ia.append(np.array([ik,iv]))

   ###### Read the header into a text file ######

   outfile = readhead(inImage)

   # find the duplicate keywords
   duplist = find_dups(outfile)

   # read the header text file line by line
   # pick out the MISSING and DUPLICATE entries 
   #
   # MISSING entries will have to be added to the header (keyword,value,comment) 
   #
   # DUPLICATE entries will be checked - to confirm the 1st entry is the one to keep, 
   # and, if so, delete the 2nd entries. 

   wv=[]
   dv=[]
   m=[]
   d=[]
   duplist=[]
   recnum = 0
   idup = 0
   istart = 0
   iend = 0
  

   duplist = find_dups(outfile)

   isodate = '0'

   with open(outfile,"r") as f:
         for line in f:
            # is there an equal sign in the keyword/value part?
            if "=" in line.split("/")[0]:
               l = (line.split("/")[0]).split("=")
               # and if not a duplicate
               if not is_duplicate(line):
                      wv.append(np.array([l[0].strip(),l[1].split("/")[0]]))
               # missing telescope info --- for date-obs and mjd-obs there may be an "=" but no value.
               if "DATE-OBS" in line:  
                      if shlex.split(line.split("/")[0],"=")[1].strip() == "":
                          kyd = shlex.split(line.split("/")[0],"=")[0].strip("=")
                          com = line.split("/")[1].strip()
                          linemod = kyd + "    /   " + com + '\n'
                          m.append(linemod)
                          missing = com + note
                          wv.append(np.array([line.split("/")[0].strip(),missing]))
               if "MJD-OBS" in line:  
                      ls = len(shlex.split(line.split("/")[0],"="))
                      if shlex.split(line.split("/")[0],"=")[ls-1].strip() == "":
                          com = line.split("/")[1].strip()
                          kyd = shlex.split(line.split("/")[0],"=")[0].strip("=")
                          linemod = kyd + "    /   " + com + '\n'
                          m.append(linemod)
                          missing = com + note
                          wv.append(np.array([line.split("/")[0].strip(),missing]))

            if "COMMENT" in line:
               wv.append(np.array(["COMMENT","VALUE"]))
            if "LDGROT" in line:
               wv.append(np.array(["LDGROT","VALUE"]))
            if "END" in line:
               wv.append(np.array(["END","VALUE"]))
            if "ISODATE" in line:
               iso = shlex.split(l[1].split("/")[0])
               isodate = iso[0]
               wv.append(np.array([l[0].strip(),isodate]))
            # missing telescope info --- no "=" and not a comment or end
            if "=" not in line.split("/")[0] and "COMMENT" not in line and "END" not in line:
               #kwd = line.split("/")[0].strip()
               com = line.split("/")[1].strip()
               m.append(line)
               missing = com + note
               wv.append(np.array([line.split("/")[0].strip(),missing]))
            # duplicates are signaled by no comment after the "/" and not in simplekeys
            if "=" in line.split("/")[0] and is_duplicate(line):
               k=(line.split("/")[0]).split("=")[0]
               kv=(line.split("/")[0]).split("=")[1]
               if idup == 0: #this is the first duplicate encountered
                      istart = recnum
                      iend = recnum
               elif idup > 0:
                      iend = recnum
               d.append(k)
               dv.append(np.array([k,kv]))
               wv.append(np.array([k,kv]))
               idup = idup + 1
            recnum = recnum+1
    
   if iend > istart:  
            print(wv[istart],wv[iend])

# print the value of the duplicate keywords to check and confirm that, in these
# cases, the first one if correct.

   L.write("\n\nDuplicate keywords start...\n")
   for i in range(len(dv)):
            L.write("%s %s\n" % (dv[i][0],dv[i][1]))
            for j in range(len(wv)):
               if dv[i][0].strip() == wv[j][0].strip():
                    if dv[i][1] != wv[j][1]:
                        print (inImage, wv[j][0], wv[j][1], dv[i][1])
                        L.write("%s: %d %s %s %s\n" % (inImage, j, wv[j][0], wv[j][1], dv[i][1]))
   print("\n")                     
   L.write("\n\nDuplicate keywords end.\n\n\n")
#

   # Now open the FITS image and store the header as hdr.

   hdul = fits.open(inImage,mode='update')
   hdr = hdul[0].header
   data = hdul[0].data

   # Provided that the duplicates are stale and we want to keep the originals,
   # then delete header lines between minid and iwcs-1
   #print ("OK to delete the following duplicate keyword/value pairs?")
   #print ("after(including) %s... %s %s and before (including) %s but excluding... %s %s \n" % (hdr[minid], w[minid], v[minid], hdr[iadd-1], w[iadd], v[iadd]))
   if idup>0:
      print ("deleting keywords between %d and %d" % (istart,iend))
      del hdr[istart:iend+1]

# loop through all of the keywords with missing values
# search for each one in the isis tcsstatus string immediately prior to the writing of the file (ti)
# result is the set: keyword, (missing)value, comment


   for i in range(len(m)):
        mk = m[i].split("/")[0].strip()
        mc = m[i].split("/")[1].strip()
        missing = mc + note
        del hdr[mk]
        mk_in_isis = 0
        if is_number(tisis) and tisis < 0:
            mval = " " # by default, empty
            #add the bare minimum, date-obs, mjd-obs, equinox (=2000) and radecsys = ("FK5")
            # date-obs from isodate or, failing that, the file-write time
            if str.lower(mk[len(mk)-3:len(mk)]) == "obs":
               if len(isodate)>1:
                  if str.lower(mk) == "date-obs":
                      mval = isodate
                  if str.lower(mk) == "mjd-obs":
                      mval = Time(isodate,format='isot',scale='utc').mjd
               else:
                  # get date from file write time
                  missing = "***  time is when the file was written"
                  ut = Time((tisis*-1),format='mjd',scale='utc').isot
                  if str.lower(mk) == "date-obs":
                      mval = ut 
                  if str.lower(mk) == "mjd-obs":
                      mval = -1*tisis
            if str.lower(mk) == "equinox":
               mval = 2000
            if str.lower(mk) == "radecsys":
               mval = "FK5"
            print("***",mk,mval,missing)
            L.write("Adding %s %s %s \n" % (mk,str(mval),missing))
            hdr[mk] = (mval,missing)
        if not is_number(tisis):
            for j in range(len(ia)):
               if mk == ia[j][0]:
                  mk_in_isis = 1
                  if is_number(ia[j][1]):
                      mval = int_or_float(ia[j][1])
                  else:
                      mval = ia[j][1]
                  print("***",mk,mval,missing)
                  L.write("Adding %s %s %s \n" % (mk,str(mval),missing))
                  hdr[mk] = (mval,missing)
            if not mk_in_isis:
                  print ("###",mk," is not in the ISIS TCSTATUS")
                  L.write("### %s is not in the ISIS TCSTATUS\n" % (mk))


   hdr['HISTORY'] = "*** indicates TCS keywords missing from original header."
   hdr['HISTORY'] = "These header cards were added with a script that looked up"
   hdr['HISTORY'] = "their values in the TCSTATUS line that was written into the"
   hdr['HISTORY'] = "MODS ISIS log at a time corresponding to the exposure start."
   hdr['HISTORY'] = "If there is no such TCSTATUS line in the log, most values are"
   hdr['HISTORY'] = "left empty, and times are when the file was written."
   print ("Updating the file")
   hdul.flush()
 
   L.write("Header updated.\n\n\n\n\n")

L.close()

