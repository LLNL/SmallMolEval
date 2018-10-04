import sys,os,glob

import os,glob
l2=[]
for e in l:
  n = e.split('/')[-1]
  l2.append(n)

for e in l2:
  os.system('python ../atomwise/remove_AVE_bias2.py -activeMols ../DUDE/'+e+'/actives_final.sdf.gz -inactiveMols ../DUDE/'+e+'/decoys_final.sdf.gz -outDir '+e)
