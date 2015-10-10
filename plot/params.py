f = open('LRSW_params', 'w')
f.write('''//LATTICE_LIBRARY="latticeLR.xml"
ALGORITHM = "LRSW"
ABSOLUTEPATH = "/home/k0010/k001007/LRI/LRSW/"
SWAPCONFIGURATION = "on"
BINARYWALKERFILE = "off"
NORMALIZATION = "off"
''')
f.write('SWEEPS = ')
f.write(str(2**10))
f.write('\nNUM_CLONES = ')
f.write(str(2**2))
f.write('''

INTERACTION = "LRI" // "LRI""SR""MF"
LATTICE = "square_lattice" // "square_lattice""triangular_lattice"
sigma = 1.0
epsilon = 1.0
method = "metropolis" // "metropolis""FTloc""SW""FTclu"
''')
f.write('dimension = ')
dim = 2
f.write(str(dim))
f.write('\n\n\n')

Nmin = 2**4
Nmax = 2**16
Tc = 0.88893
Trange = 0.4 #T is obtained from Tc-Trange to Tc+Trange
Tdiv = 10
nu = 1.0
siteMag = 4.0
#TrangeMag = 2.0 
TrangeMag = siteMag**(1.0/(dim*nu))


Lreal = Nmin**(1.0/float(dim))
L = int(round(Lreal))
if L%2==1:
  L = L-1
while L <= int(round(Nmax**(1.0/float(dim)))):
  Lstr = 'L = ' + str(L) + '\n'
  f.writelines(Lstr)
  Lreal = Lreal*siteMag**(1.0/float(dim))
  L = int(round(Lreal))
  if L%2==1:
    L = L-1
  T = Tc - Trange
  i = 0
  while i < Tdiv + 1:
    Tstr = '{ T = ' + str(T) + ' }\n'
    f.write(Tstr)
    T += Trange*2.0/Tdiv
    i += 1
  Trange /= TrangeMag
  f.write('\n')

f.close()


