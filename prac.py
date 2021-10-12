import Bio.PDB
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

phi_psi = ([0,0])
phi_psi = np.array(phi_psi)
pdb1 ='/home/ashwini/python_project/1fat.pdb'

for model in Bio.PDB.PDBParser().get_structure('1fat',pdb1) :
    for chain in model :
        polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
        for poly_index, poly in enumerate(polypeptides) :
            print "Model %s Chain %s" % (str(model.id), str(chain.id)),
            print "(part %i of %i)" % (poly_index+1, len(polypeptides)),
            print "length %i" % (len(poly)),
            print "from %s%i" % (poly[0].resname, poly[0].id[1]),
            print "to %s%i" % (poly[-1].resname, poly[-1].id[1])
            phi_psi = poly.get_phi_psi_list()
            for res_index, residue in enumerate(poly) :
                #res_name = "%s%i" % (residue.resname, residue.id[1])
                #print res_name, phi_psi[res_index]
                phi_psi = np.vstack([phi_psi \
                ,np.asarray(phi_psi[res_index])]).astype(np.float) 
                #np.float - conversion to float array from object

phi, psi = np.transpose(phi_psi)

phi = np.degrees(phi)
psi = np.degrees(psi)

phi = phi[~np.isnan(phi)] # avoiding nan
psi = psi[~np.isnan(psi)]

f,ax = plt.subplots(1)
plt.title('RAMACHANDRAN PLOT')
plt.xlabel('$\phi^o$', size=20,fontsize=15)
plt.ylabel('$\psi^o$ ', size=20,fontsize=15)

h=ax.hexbin(phi, psi,  extent=[-200,200,-200,200],cmap=plt.cm.Greens)
#h=ax.hexbin(phi, psi,gridsize=100,  extent=[-180,180,-180,180],cmap=plt.cm.Blues)

f.colorbar(h)
plt.grid()
plt.show()