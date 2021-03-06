import DescriptorSets as DS
import sys

a1=DS.get_mols(sys.argv[1])

a1_simple = DS.get_simple(a1)


DS.remove_ddups(a1_simple)

#active to decoy distances
#d2a = DS.dist_alt(a1_simple, sys.argv[2], file=True)

print(DS.get_dsize(a1_simple))


#active to active distances
a2a = DS.dist(a1_simple,a1_simple, same=True)



DS.writeG(a2a,5,sys.argv[2]+'.out.txt')
