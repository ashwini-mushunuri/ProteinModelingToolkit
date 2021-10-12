#Authors: Ashwini Mushunuri & Sagarika Chakraborty
# -*- coding: utf-8 -*-
from prody import *
from pylab import *
from Tkinter import *
import tkMessageBox
from tkFileDialog   import askopenfilename
import tkFileDialog as fileDialog
from Bio.PDB import PDBParser
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import Bio.PDB
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from Bio import SeqIO
import os



#################################################################################################



def OpenFile():
	global opt
	opt = fileDialog.askopenfilename()
	'''FileDialog is a module with open and save dialog functions.
    .askopenfilename: To open file: Dialog that requests selection of an existing file.'''
   	global ext
   	global handle
   	handle = open(opt, 'r')
	global file1
	file1 = open("file.txt","a")


def fasta():

	for record in SeqIO.parse(handle, "pdb-seqres") :  #SeqIO.parse(): takes a file handle (or filename) and format name, and returns a SeqRecord iterator.
		x = ">" + record.id + "\n" + record.seq     #pdb-seqres:  primary sequence of the polymeric molecules present in the entry
		label = Label(root, text= str(x))
		seq1 = str(record.seq)
		file1.write(seq1)
		print x
		label.pack()

'''The ProteinAnalysis class takes one argument, the protein sequence as a string and builds
a sequence object using the Bio.Seq module'''

def mol_weight():
	f = open('file.txt', 'r')
	my_seq = f.read()
	analysed_seq = ProteinAnalysis(my_seq)
	weight = analysed_seq.molecular_weight()  #Calculates the molecular weight of a protein.
	tkMessageBox.showinfo("molecular weight", weight)


def amino_per():
	f = open('file.txt', 'r')
	my_seq = f.read()
	analysed_seq = ProteinAnalysis(my_seq)
	amino =  analysed_seq.get_amino_acids_percent()   #returns the number of amino acid in percentage of entire sequence
	amino = str(dict.items(amino))
	tkMessageBox.showinfo("amino acid percentage", amino)

def gravity():
	f = open('file.txt', 'r')
	my_seq = f.read()
	analysed_seq = ProteinAnalysis(my_seq)
	gravy = analysed_seq.gravy()     #the sum of hydropathy values of all the amino acids, divided by the number of residues in the sequence.
	tkMessageBox.showinfo("gravy", gravy)

def fraction():
	f = open('file.txt', 'r')
	my_seq = f.read()
	analysed_seq = ProteinAnalysis(my_seq)   #This methods returns a list of the fraction of amino acids which tend to be in helix, turn or sheet.
	fraction = analysed_seq.secondary_structure_fraction()
	tkMessageBox.showinfo("secondary structure fraction [Helix, Turn, Sheet].", fraction)

	
def isoelectric():
	f = open('file.txt', 'r')
	my_seq = f.read()
	analysed_seq = ProteinAnalysis(my_seq)
	iso = analysed_seq.isoelectric_point()   #calculates pl of a protein structure
	tkMessageBox.showinfo("isoelectric point", iso)


def instability():
	f = open('file.txt', 'r')
	my_seq = f.read()
	analysed_seq = ProteinAnalysis(my_seq)
	index =  analysed_seq.instability_index()   #tests a protein for stability.
	tkMessageBox.showinfo("instabilty index", index)

def aromaticity():
	f = open('file.txt', 'r')
	my_seq = f.read()
	analysed_seq = ProteinAnalysis(my_seq)
	aro = analysed_seq.aromaticity()	    #Calculates the aromaticity value of a protein 
	tkMessageBox.showinfo("aromaticity", aro)


def viz():
	base=os.path.basename(opt)
	inp = os.path.splitext(base)[0]
	prot = parsePDB(inp) 
	showProtein(prot)
	raw_input("press enter to continue...")

def rama():
	phi_psi = ([0,0])
	phi_psi = np.array(phi_psi)
	pdb1 = opt

	for model in Bio.PDB.PDBParser().get_structure('XYZ',pdb1) :
    		for chain in model:
        		polypeptides = Bio.PDB.PPBuilder().build_peptides(chain) # extract polypeptides.
        		for poly_index, poly in enumerate(polypeptides) :  #enumerate() : allows us to loop over something and have an automatic counter.
            			phi_psi = poly.get_phi_psi_list() # Return the list of phi/psi dihedral angles.
            			for res_index, residue in enumerate(poly) : 
                			#res_name = "%s%i" % (residue.resname, residue.id[1])
                			#print res_name, phi_psi[res_index]
                			phi_psi = np.vstack([phi_psi \
                			,np.asarray(phi_psi[res_index])]).astype(np.float)  #vstack() : Stack arrays in sequence vertically (row wise)
                			#np.float - conversion to float array from object  #astype(): Copy of the array, cast to a specified type.

	phi, psi = np.transpose(phi_psi)  #Returns a view of the array with axes transposed.

	phi = np.degrees(phi)  #converts angle x from radians to degrees.
	psi = np.degrees(psi)

	phi = phi[~np.isnan(phi)] # math library function, avoids nan, used to determine whether a given parameter is a valid number or not.
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

def About():
	tkMessageBox.showinfo("About","Protein analysis Toolkit\nProtein fractions:\nAmino acids in helix: V, I, Y, F, W, L.\nAmino acids in turn: N, P, G, S.\nAmino acids in sheet: E, M, A, L.")


def NewFile():
    print "New File!"
    global handle
    handle.close()
    opt = fileDialog.askopenfilename()   
    handle = open(opt, 'r')




#######################################################################################


root = Tk() #creates a GUI interface
'''Creating an instance of Tk initializes this interpreter and creates the root window.'''
root.geometry('1000x1000') # sets dimensions of the window
menu = Menu(root) #creates a menubar on the window
root.config(menu=menu)
filemenu = Menu(menu)
menu.add_cascade(label="File", menu=filemenu)  #add.cascade() : adds new drop down menu
filemenu.add_command(label="New", command=NewFile)
filemenu.add_command(label="Open...", command=OpenFile)
filemenu.add_separator()
filemenu.add_command(label="Exit", command=root.quit)

operations = Menu(menu)
menu.add_cascade(label="operations", menu=operations)
operations.add_command(label="fasta sequence", command=fasta)
operations.add_command(label="molecular weight", command=mol_weight)
operations.add_command(label="amino acid percentage", command=amino_per)
operations.add_command(label="gravity", command=gravity)
operations.add_command(label="secondary structure fraction", command=fraction)
operations.add_command(label="isoelectricity point", command=isoelectric)
operations.add_command(label="instability index", command=instability)
operations.add_command(label="aromaticity", command=aromaticity)

plots = Menu(menu)
menu.add_cascade(label="plots", menu=plots)
plots.add_command(label="visualize protein", command=viz)
plots.add_command(label="ramachandran plot", command=rama)

helpmenu = Menu(menu)
menu.add_cascade(label="Help", menu=helpmenu)
helpmenu.add_command(label="About...", command=About)

#button = Button(root, text="Print Me", command=fasta)
#button.pack()


widget = Label(font=("Times New Roman", 40), text='PROTEIN ANALYSIS TOOLKIT', fg="Blue") #label(): prints text on GUI interface
widget.pack()   
mainloop()
