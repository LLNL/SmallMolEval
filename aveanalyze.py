import os, sys, glob

#directory of dataset files split by AVE bias
dir=sys.argv[1]

l=glob.glob(dir+'/*')
n=dir.split('/')[-1]
print(l)

for e in l:
  e2=e.split('/')[-1]
  os.system('python ~/DataSets/atomwise/analyze_AVE_bias.py -activeMolsTraining '+e+'/actives.T.smi -inactiveMolsTraining '+e+'/inactives.T.smi -activeMolsTesting '+e+'/actives.V.smi -inactiveMolsTesting '+e+'/inactives.V.smi -outFile '+n+e2)
