#######################################################################################################################
# Python script to construct and run simulations of intrinsically disordered proteins
# and polynucleotides using the extended HPS model[1] 
# [1]  Sequence dependent co-phase separation of protein-polynucleotide mixtures elucidated using molecular simulations
#     RM Regy, GL Dignon, W Zheng, YC Kim, J Mittal - Nucleic Acids Research
#
# Authors: Roshan M Regy, Wenwei Zheng, Gregory L. Dignon, Jeetain Mittal 
#######################################################################################################################
# Prerequisities: This code was tested using HOOMD v2.9.3 with azplugins installed and GSD v2.2.0 
#
# Usage: python polynucleotide_protein_1chain.py One_letter_amino_acid_sequence_file length_of_polyA
#
# This code has the following steps:
#
# (1) Create a single chain initial configuration of the protein and polynucleotide as linear polymer chains --> generates start.gsd 
#
# (2) Run a NVT simulation at user defined temperature and time step and for user defined number of steps
#   --> generates production.dcd, restart_tmp1.gsd, restart_tmp2.gsd, production.log, stress.log  
#######################################################################################################################
import sys,os,numpy as np
import hoomd, hoomd.md as md
import hoomd.deprecated as old
from hoomd import azplugins
import gsd, gsd.hoomd, gsd.pygsd 
aaparams='''#AA     Mass    Charge  Sigma   Lambda
ALA     71.08   0.00    5.040   0.730
ARG     156.20  1.00    6.560   0.000
ASN     114.10  0.00    5.680   0.432
ASP     115.10  -1.00   5.580   0.378
CYS     103.10  0.00    5.480   0.595
GLN     128.10  0.00    6.020   0.514
GLU     129.10  -1.00   5.920   0.459
GLY     57.05   0.00    4.500   0.649
HIS     137.10  0.5    6.080   0.514
ILE     113.20  0.00    6.180   0.973
LEU     113.20  0.00    6.180   0.973
LYS     128.20  1.00    6.360   0.514
MET     131.20  0.00    6.180   0.838
PHE     147.20  0.00    6.360   1.000
PRO     97.12   0.00    5.560   1.000
SER     87.08   0.00    5.180   0.595
THR     101.10  0.00    5.620   0.676
TRP     186.20  0.00    6.780   0.946
TYR     163.20  0.00    6.460   0.865
VAL     99.07   0.00    5.860   0.892'''
# Length of polynucleotide
polyAlength = int(sys.argv[2])
# PRODUCTION RUN PARAMETERS
production_dt=0.01 # Time step for production run in picoseconds
production_steps=100000000 # Total number of steps 
production_T=190 # Temperature for production run in Kelvin



seq={'R':'ARG','H':'HIS','K':'LYS','D':'ASP','E':'GLU',
     'S':'SER','T':'THR','N':'ASN','Q':'GLN','C':'CYS',
     'U':'SEC','G':'GLY','P':'PRO','A':'ALA','V':'VAL',
     'I':'ILE','L':'LEU','M':'MET','F':'PHE','Y':'TYR',
     'W':'TRP'}

# ##### Read one letter amino acid sequence from file
filein = sys.argv[1]
fileout='%s_seq3.dat'%(filein)
nline=1
count=0
fout=open(fileout,'w')
with open(filein,'r') as fid:
    for i in fid:
        if i[0]!='#':
            for j in i:
                if j in seq:
                    fout.write(' %s'%seq[j])
                    count+=1
                    if count==nline:
                        fout.write('\n')
                        count=0
fout.close()


# #### 1.2 Read sequence and force field parameters
# ##### Input parameters for all the amino acids (force field)
ff_para = 'stats_module.dat'
aalist={}
with open(ff_para,'r') as fid:
    for i in fid:
        if i[0]!='#':
            tmp=i.rsplit()
            aalist[tmp[0]]=np.loadtxt(tmp[1:],dtype=float)
aakeys=list(aalist.keys())
# This translates each amino acid type into a number, which will be used in HOOMD
# For example, GLY is with an ID of 10
aamass=[]
aacharge=[]
aaradius=[]
aahps=[]
print ('aakeys')
print (aakeys)
print ('aalist[i][1]')
print (aalist[aakeys[1]][1])
for i in aakeys:
    aamass.append(aalist[i][0])
    aacharge.append(aalist[i][1])
    aaradius.append(aalist[i][2])
    aahps.append(aalist[i][3])
# Now we can translate the entire sequence into a number code according to the order in 'aakeys'
chain_id=[]
chain_mass=[]
chain_charge=[]
with open(fileout,'r') as fid:
    for i in fid:
        iname=i.rsplit()[0]
        chain_id.append(aakeys.index(iname))
        chain_mass.append(aalist[iname][0])
        chain_charge.append(aalist[iname][1])

pbond_length=0.38
chain_length=len(chain_id)
# Add RNA chains 
rbond_length=0.5
for i in range(polyAlength):
    chain_id.append(len(aakeys)-1)
    chain_mass.append(329.2)
    chain_charge.append(-1)

#box_length=np.max([pbond_length*chain_length,rbond_length*polyAlength])+10
box_length=200
print ('box length')
print (box_length)
# #### 1.3 Now we can build HOOMD data structure for one single frame
s=gsd.hoomd.Snapshot()
s.particles.N = chain_length+polyAlength
print ('aakeys')
print (aakeys)
s.particles.types = aakeys
s.particles.typeid = chain_id
s.particles.mass = chain_mass
s.particles.charge = chain_charge

#  Build initial position as a linear chain
pos=[]
for i in range(chain_length):
    # Change the z-coordinate to have a linear chain
    pos.append((0,0,(i-int(chain_length/2))*pbond_length))
    print ((i-int(chain_length/2))*pbond_length)

for i in range(polyAlength):
    pos.append((5,0,(i-int(polyAlength/2))*rbond_length))
    print ((i-int(polyAlength/2))*rbond_length)
pos=np.array(pos)
s.particles.position= pos

# Initialize bond
nbonds=chain_length-1+polyAlength-1
s.bonds.N = nbonds
s.bonds.types = ['AA_bond','NT_bond']
s.bonds.typeid = [0]*(chain_length-1)+[1]*(polyAlength-1)
bond_pairs=np.zeros((nbonds,2),dtype=int)
for i in range(0,chain_length-1):
    print ('%s-%s-A'%(i,i+1))
    bond_pairs[i,:] = np.array([i,i+1])
for cnt,i in enumerate(range(chain_length,nbonds+1)):
    print ('%s-%s-B'%(i,i+1))
    bond_pairs[cnt+chain_length-1,:] = np.array([i,i+1])
s.bonds.group = bond_pairs
print (bond_pairs)
# Box size
s.configuration.dimensions=3
s.configuration.box=[box_length,box_length,box_length,0,0,0]
s.configuration.step=0

# #### 1.4 Write intial singe chain gsd file
f = gsd.hoomd.open(name='start.gsd', mode='wb')
f.append(s)
f.close()
#################################################################################################
# ## start.gsd contains one single chain of the given protein 
# ----------------------------------------------------------------------------------------------
# ## 4.0. Run a simulation using start.gsd created in the previous step
################################################################################################
hoomd.context.initialize()
system = hoomd.init.read_gsd('start.gsd')

n_steps = production_steps # 1 microseconds

fileroot = 'Production'
nl = hoomd.md.nlist.cell()

## Bonds
harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('AA_bond', k=8360, r0=0.381)
harmonic.bond_coeff.set('NT_bond', k=8360, r0=0.5)
## Nonbonded
nl.reset_exclusions(exclusions=['1-2', 'body'])
nb = azplugins.pair.ashbaugh(r_cut=0, nlist=nl)
for i in aakeys:
    for j in aakeys:
        nb.pair_coeff.set(i,j,lam=(aalist[i][3]+aalist[j][3])/2.,
                          epsilon=0.8368, sigma=(aalist[i][2]+aalist[j][2])/10./2.,r_cut=2.0)    

## Electrostatics
yukawa = hoomd.md.pair.yukawa(r_cut=0.0, nlist=nl)
for i,atom1 in enumerate(aakeys):
    for j,atom2 in enumerate(aakeys):
        yukawa.pair_coeff.set(atom1,atom2,epsilon=aalist[atom1][1]*aalist[atom2][1]*1.73136, kappa=1.0, r_cut=3.5) 

## Group Particles
all = hoomd.group.all()

## Set up integrator
hoomd.md.integrate.mode_standard(dt=production_dt) # Time units in ps
temp = production_T*0.00831446
integrator = hoomd.md.integrate.langevin(group=all, kT=temp, seed=399991) # Temp is kT/0.00831446
for cnt,i in enumerate(aakeys):
    integrator.set_gamma(i,gamma=aamass[cnt]/1000.0)
## Outputs
hoomd.analyze.log(filename=fileroot+'.log', quantities=['potential_energy', 'pressure_xx', 'pressure_yy', 'pressure_zz', 'temperature','lx','ly','lz'], period=100000, overwrite=False, header_prefix='#')
hoomd.analyze.log(filename='stress.log', quantities=['pressure_xy', 'pressure_xz', 'pressure_yz'], period=100000, overwrite=False, header_prefix='#') # Output stress tensor?
hoomd.dump.gsd('restart_tmp1.gsd', period=1000000, group=all, truncate=True)
hoomd.dump.gsd('restart_tmp2.gsd', period=1000000, group=all, truncate=True, phase=500000)
hoomd.dump.dcd(fileroot+'_dump.dcd', period=100000, group=all, overwrite=False)

## Run simulation
hoomd.run_upto(production_steps, limit_hours=48)
########################################################################################################
