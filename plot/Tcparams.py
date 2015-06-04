f = open('LRSW_params', 'w')
f.write('''//LATTICE_LIBRARY="latticeLR.xml"
ALGORITHM = "LRSW"
ABSOLUTEPATH = "/home/k0010/k001007/LRI/LRSW/"
SWAPCONFIGURATION = "on"
BINARYWALKERFILE = "off"
NORMALIZATION = "on"
''')
f.write('SWEEPS = ')
f.write(str(2**10))
f.write('\nNUM_CLONES = ')
f.write(str(2**4))
f.write('''

INTERACTION = "LRI" // "LRI""nearest""meanfield"
LATTICE = "square_lattice" // "square_lattice""triangular_lattice"
sigma = 1.0
''')
f.write('dimension = ')
dim = 2
f.write(str(dim))
f.write('\n\n\n')

Nmin = 2**2
Nmax = 2**20
Tc = 0.88893
siteMag = 4.0


Tstr = 'T = ' + str(Tc) + '\n'
f.write(Tstr)
Lreal = Nmin**(1.0/float(dim))
L = int(round(Lreal))
if L%2==1:
  L = L-1
while L <= int(round(Nmax**(1.0/float(dim)))):
  Lstr = '{ L = ' + str(L) + ' }\n'
  f.writelines(Lstr)
  Lreal = Lreal*siteMag**(1.0/float(dim))
  L = int(round(Lreal))
  if L%2==1:
    L = L-1
f.write('\n')

f.close()


