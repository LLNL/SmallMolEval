import glob,os,sys

arg1 = '240'
arg2 = '180'

if len(sys.argv) == 3:
  arg1 = sys.argv[1]
  arg2 = sys.argv[2]

l=glob.glob('*.out.txt')

nl=[]

for e in l:
  nl.append(e.split('.')[0])

for e in nl:
  os.system("~/DataSets/scripts/gf.plot '"+e+"' '"+arg1+"' '"+arg2+"'")
