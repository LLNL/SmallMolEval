from rdkit import Chem
import gzip, sys, glob
from rdkit.Chem import rdqueries
from rdkit.Chem import Descriptors
from rdkit.Chem import rdmolops

from math import*
import decimal

#get list of molecules from gzipped sdf
def get_mols(fhandle):
  inf = gzip.open(fhandle)
  gzsuppl = Chem.ForwardSDMolSupplier(inf)
  ms = [x for x in gzsuppl if x is not None]
  inf.close()
  return ms

#for m in ms:
#  print(m.GetProp('_Name'))

#get simple features for molecule list
#returns dictionary of descriptor with molecule name as key 
#and list of unique vectors, for multiple entries of same molecule name 
def get_simple(ms):
  errors = 0
  descr = {}
  for m in ms:
    try:
      name = m.GetProp('_Name')
      if name not in descr:
        descr[name]=[]

      #get descriptor vector
      dv = []
      dv.append(m.GetNumHeavyAtoms())
 
      a_nums = [5,35,6,17,9,53,7,8,15,16]
      for n in a_nums:
        q = rdqueries.AtomNumEqualsQueryAtom(n)
        dv.append(len(m.GetAtomsMatchingQuery(q)))

      dv.append(Descriptors.NumHAcceptors(m))
      dv.append(Descriptors.NumHDonors(m))
      dv.append(Descriptors.MolLogP(m))
      dv.append(Descriptors.RingCount(m))

      #print(rdmolops.AssignAtomChiralTagsFromStructure(m))

      descr[name].append(dv)
    except ValueError:
      if len(descr[name]) == 0:
        del descr[name]
      errors = errors +1
  print(str(errors)+' ValueError(s) has(have) occured')
  return descr

#normalize descriptor vectors
#def norm(descr)

#return number of entries in descriptor dictionary
def get_dsize(descr):
  total = 0
  for id in descr:
    total = total + len(descr[id])
  return total

# remove duplicate description vectors
def remove_ddups(descr):
  for id in descr:
    newList = []
    for i in range(len(descr[id])):
      if descr[id][i] not in newList:
        newList.append(descr[id][i])
    descr[id] = newList 

#I is 1 if distance is less than t
def I(d,t):
  if d < t:
    return 1
  else:
    return 0

#nearest neighbor and empty space functions
#pass active to active distance list for nearest neighbor
#pass decoy to active distance list for empty space
def GF(dist,t):
  n = len(dist)
  sumI = 0.0
  for i in range(n):
    d = min(dist[i])
    #print(d,I(d,t))
    sumI = sumI + I(d,t)
  return sumI/n

#write data for graph and get all Sum values
def writeS(a2a,d2a,trange, fhandle):
  sumG = 0
  sumF = 0
  sumS = 0
  f=open(fhandle, 'w')

  for t in drange(0,trange,.1):
    G = GF(a2a,t)
    F = GF(d2a,t)
    sumG = sumG + G
    sumF = sumF + F
    sumS = sumS + (F-G)
    f.write(str(t)+'\t'+str(G)+'\t'+str(F)+'\n')
  print(sumG,sumF,sumS)
  f.close()

#write data for graph and get all Sum values
def writeG(a2a,trange, fhandle):
  sumG = 0
  f=open(fhandle, 'w')

  for t in drange(0,trange,.1):
    G = GF(a2a,t)
    sumG = sumG + G
    f.write(str(t)+'\t'+str(G)+'\n')
  print(sumG)
  f.close()



#    i = euclidean_distance(i,j)
 
#get distances
#2 molecule lists
#All if every representation used, 
#default is Fasle and only first representation is used
def dist(ms1,ms2,same=True,all=False):
  #first array dimension
  dist = []
  #for each molecule name in the dictionary
  m = 0
  for m1 in ms1:
    #if all=True get length of entries, else only use first entry
    if all: i = len(ms1[m1])
    else: i = 1
    #for each vector
    for j in range(i):
      dist.append([])
      n = 0
      for m2 in ms2:
        if all: k = len(ms2[m2])
        else: k = 1
        #nnd = float('inf')
        for l in range(k):
          if same and m == n:
            dist[m].append(float('inf'))
          else:
            #print(m,n)
            dist[m].append(euclidean_distance(ms1[m1][j],ms2[m2][l]))
          n = n + 1
      m = m + 1
  return dist

#get distance array for large decoy sets
#to be used when comparing to a large decoy set
#decoy files split and pass directory to decoys
#only use first active
def dist_alt(ms_a,decoy_d,file=True):
  #list of minimum distances between each decoy and the active set
  dist = []

  #decoy num
  j = 0

  if file == False:
    #get list of decoy files
    dlist = glob.glob(decoy_d+'/*.sdf.gz')
  else:
    dlist = []
    dlist.append(decoy_d)
  for df in dlist:
    ms_d = get_mols(df)
    d_descr = get_simple(ms_d)
    #initialize list to infinity for each decoy
    for i in range(len(d_descr)):
      dist.append([float('inf')])
    for d in d_descr:
      for a in ms_a:
        if euclidean_distance(ms_a[a][0],d_descr[d][0]) < dist[j][0]:
          dist[j][0] = euclidean_distance(ms_a[a][0],d_descr[d][0])
      j = j + 1

  return dist


def euclidean_distance(x,y):
  return sqrt(sum(pow(a-b,2) for a, b in zip(x, y)))
  
def drange(x, y, jump):
  while x < y:
    yield float(x)
    x += decimal.Decimal(jump)

#get similarities of vectors, 
#i=1 is use only first vector for each compound, otherwise use all
# s is similarity measure
# autoscale, on or off
#def get_sim(descr,i=1,s=,autoscale=on):


## Verified that atom count and heavy ato count always gives same
#for id in descr:
#  for i in range(len(descr[id])):
#    if descr[id][i][0] != descr[id][i][1]:
#      print(id)
