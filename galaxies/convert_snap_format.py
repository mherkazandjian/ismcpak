# - converts fi snapthosts from version zero to versions 1
# - copy this script into the dir containing all the snapshots and run it
#
#------------------------------------------------------
import subprocess
import time

#snaps = range(21)
snaps = range(5)
#snaps = [18]

#------------------------------------------------------
def gen_runinfo(snap_idx):
    
    fObjR = open('runinfo_template','r')
    fObjW = open('runinfo', 'w')
    
    for line in fObjR:
        
        if 'INDEX' in line:
            lineNew = line.replace('INDEX', '%06d' % snap_idx )
            
        else:
            lineNew = line
            
        fObjW.write(lineNew)

    fObjR.close()
    fObjW.close()
    
for snap in snaps:
    print 'converting snapshot %d' % snap
    #getting the runinfo file from the template for this snapshot
    gen_runinfo(snap)
    
    #running fi to convert the snapshot
    subprocess.Popen(["rm","-fvr","simlog"])
    p = subprocess.Popen("./fi")
    p.wait()
    
    time.sleep(0.1)
    subprocess.Popen(["mv","tmp.000000", "fiout.%06d" % snap])
    
