import sys
import os

fdir = "coi_compare"

if len(sys.argv) < 2:
	print "Usage: python copy.py <tag>"
	sys.exit()

cmd1 = "cp IO/TempOut.dat IO/" + fdir + "/TempOut_" + sys.argv[1] + ".dat"
cmd2 = "cp IO/clima_allout.tab IO/" + fdir + "/clima_allout_" + sys.argv[1] + ".tab"
#cmd3 = "cp IO/TempOut.dat IO/TempIn.dat"

os.popen4(cmd1)
os.popen4(cmd2)
#os.popen4(cmd3)

#print "Tag:  ", sys.argv[1], "mb"
print cmd1
print cmd2
#print cmd3
