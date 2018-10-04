import DescriptorSets as DS
import sys

a1=DS.get_mols(sys.argv[1])

a1_simple = DS.get_simple(a1)

d1=DS.get_mols(sys.argv[2])

d1_simple= DS.get_simple(d1)

#print(DS.get_dsize(a1_simple))

DS.remove_ddups(a1_simple)
DS.remove_ddups(d1_simple)

#print(DS.get_dsize(a1_simple))

#print(DS.I([1,1,1],[1,-1,1],2))

#active to active distances
a2a = DS.dist(a1_simple,a1_simple, same=True)

#active to decoy distances
d2a = DS.dist(d1_simple,a1_simple,same=False)


DS.writeS(a2a,d2a,5,sys.argv[3]+'.out.txt')
