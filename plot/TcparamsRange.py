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
f.write(str(2**2))
f.write('''

INTERACTION = "LRI" // "LRI""nearest""meanfield"
LATTICE = "square_lattice" // "square_lattice""triangular_lattice"
sigma = 1.0
''')
f.write('dimension = ')
dim = 2
f.write(str(dim))
f.write('\n\n\n')

Nmin = 2**4
Nmax = 2**16
Tc = 0.88893
Tcsigma = 0.00000005 #T is obtained from Tc-Trange to Tc+Trange
Tcdiv = 5
siteMag = 4.0


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
  i = 0
  while i < Tcdiv:
    Tstr = '{ T = ' + str(Tc - Tcsigma + 2.0*Tcsigma*i/(Tcdiv-1) ) + ' }\n'
    f.write(Tstr)
    i += 1
  f.write('\n')

f.close()


