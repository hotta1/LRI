import subprocess
import math

#physicalQuantityList = ['binder','magsussca','absM','specheat']
physicalQuantityList = ['binder','magsussca','absM']

Tc = 8.56
nuinv = 2.0
omega = 0.0
gamma = 3.0
beta = -0.3
alpha = 0.1
zero = 0.0
theta0 = 1.0
theta1 = 1.0
theta2 = 1.0
theta3 = 1.0
theta4 = 1.0

numsizes = 3


filename = ''
for physicalQuantity in physicalQuantityList:
  filename += physicalQuantity + '_'
print filename
opFile = open(filename+'.dat','w')

numdata = 0
for physicalQuantity in physicalQuantityList:
  ipFile = open(physicalQuantity+'.dat','r')
  for line in ipFile:
    lineList = line[:-1].split()
    lineList.insert(0,str(numdata))
    opFile.write(" ".join(lineList))
    opFile.write("\n")
  ipFile.close()
  numdata +=1
opFile.close()

KHop = open(filename+'.op','w')
#subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870/CC2/new_bfss", "-n", str(i+1), "-f", "1", "-c", filename+'.dat', "1", str(Tc), "1", str(nuinv), "1", str(omega), "1", str(gamma)], stdout = KHop)
subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870/CC2/new_bfss", "-c", "-n", str(numdata), filename+'.dat', "1", str(Tc), "1", str(nuinv), "0", str(zero), "1", str(theta0), "1", str(theta1), "1", str(theta2), "1", str(gamma), "1", str(theta0), "1", str(theta1), "1", str(theta2), "1", str(beta), "1", str(theta0), "1", str(theta1), "1", str(theta2), "1", str(alpha), "1", str(theta0), "1", str(theta1), "1", str(theta2)], stdout = KHop)
KHop.close()
KHop2 = open(filename+"correction.op",'w')
subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870/CC2/new_bfss", "-c", "-f", "1", "-n", str(numdata), filename+'.dat', "1", str(Tc), "1", str(nuinv), "1", str(omega), "0", str(zero), "1", str(theta0), "1", str(theta1), "1", str(theta2), "1", str(theta3), "1", str(theta4), "1", str(gamma), "1", str(theta0), "1", str(theta1), "1", str(theta2), "1", str(theta3), "1", str(theta4), "1", str(beta), "1", str(theta0), "1", str(theta1), "1", str(theta2), "1", str(theta3), "1", str(theta4), "1", str(alpha), "1", str(theta0), "1", str(theta1), "1", str(theta2), "1", str(theta3), "1", str(theta4)], stdout = KHop2)
KHop2.close()


ipFile = open(filename+'.dat','r')
sizeList = []
ipList = []
for line in ipFile:
  lineList = line[:-1].split()
  ipList.append(lineList)
  if sizeList.count(lineList[1]) == 0:
    sizeList.append(lineList[1])
ipFile.close

print sizeList

pairList = []
for x in range(0,len(sizeList)-numsizes+1):
  ipPair = []
  for y in range(x,x+numsizes):
    ipPair.append(sizeList[y])
  pairList.append(ipPair)

print pairList

separatedFNList = []
for pair in pairList:
  opElement = []
  opElement.append(filename+pair[numsizes-numsizes/2-1])
  opElement.append(pair[numsizes-numsizes/2-1])
  separatedFNList.append(opElement)
  opFile = open(filename+pair[numsizes-numsizes/2-1]+'.dat','w')
  for ipLine in ipList:
    for pairelement in pair:
      if ipLine[1]==pairelement:
        opFile.write(" ".join(ipLine))
        opFile.write("\n")
  opFile.close()



print separatedFNList

for separatedFN in separatedFNList:
  KHop = open(separatedFN[0]+'.op','w')
  subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870/CC2/new_bfss", "-c", "-n", str(numdata), separatedFN[0]+'.dat', "1", str(Tc), "1", str(nuinv), "0", str(zero), "1", str(theta0), "1", str(theta1), "1", str(theta2), "1", str(gamma), "1", str(theta0), "1", str(theta1), "1", str(theta2), "1", str(beta), "1", str(theta0), "1", str(theta1), "1", str(theta2), "1", str(alpha), "1", str(theta0), "1", str(theta1), "1", str(theta2)], stdout = KHop)
  KHop.close()


TcList = []
nuinvList = []
C2_0List = []
C2_1List = []
C2_2List = []
C2_3List = []
for separatedFN in separatedFNList:
#  KHrun(separatedFN[0], physicalQuantity)
  KHip = open(separatedFN[0]+".op",'r')
  TcLineList = ['','','']
  nuinvLineList = ['','','']
  C2_0LineList = ['','','']
  C2_1LineList = ['','','']
  C2_2LineList = ['','','']
  C2_3LineList = ['','','']
  TcLineList[0] = separatedFN[1]
  nuinvLineList[0] = separatedFN[1]
  C2_0LineList[0] = separatedFN[1]
  C2_1LineList[0] = separatedFN[1]
  C2_2LineList[0] = separatedFN[1]
  C2_3LineList[0] = separatedFN[1]
  for line in KHip:
    ipLineList = line[:-1].split()
    if len(ipLineList)>=4:
      """
      if ipLineList[1]=='LMIN':
        TcLineList[0] = ipLineList[3]
        nuinvLineList[0] = ipLineList[3]
        C2LineList[0] = ipLineList[3]
      """
      if ipLineList[1]=='p[0]':
        if len(ipLineList)==4:
          TcLineList[1] = ipLineList[3]
        elif len(ipLineList)==5:
          TcLineList[1] = ipLineList[3]
          TcLineList[2] = ipLineList[4]
      if ipLineList[1]=='p[1]':
        if len(ipLineList)==4:
          nuinvLineList[1] = ipLineList[3]
        elif len(ipLineList)==5:
          nuinvLineList[1] = ipLineList[3]
          nuinvLineList[2] = ipLineList[4]
      if ipLineList[1]=='p[2]':
        if len(ipLineList)==4:
          C2_0LineList[1] = ipLineList[3]
        elif len(ipLineList)==5:
          C2_0LineList[1] = ipLineList[3]
          C2_0LineList[2] = ipLineList[4]
      if ipLineList[1]=='p[6]':
        if len(ipLineList)==4:
          C2_1LineList[1] = ipLineList[3]
        elif len(ipLineList)==5:
          C2_1LineList[1] = ipLineList[3]
          C2_1LineList[2] = ipLineList[4]
      if ipLineList[1]=='p[10]':
        if len(ipLineList)==4:
          C2_2LineList[1] = ipLineList[3]
        elif len(ipLineList)==5:
          C2_2LineList[1] = ipLineList[3]
          C2_2LineList[2] = ipLineList[4]
      if ipLineList[1]=='p[14]':
        if len(ipLineList)==4:
          C2_3LineList[1] = ipLineList[3]
        elif len(ipLineList)==5:
          C2_3LineList[1] = ipLineList[3]
          C2_3LineList[2] = ipLineList[4]
  TcList.append(TcLineList)
  nuinvList.append(nuinvLineList)
  C2_0List.append(C2_0LineList)
  C2_1List.append(C2_1LineList)
  C2_2List.append(C2_2LineList)
  C2_3List.append(C2_3LineList)
LvsTc = open(filename+"Tc-L.dat",'w')
Lvsnuinv = open(filename+"nuinv-L.dat",'w')
LvsC2_0 = open(filename+"C2_0-L.dat",'w')
LvsC2_1 = open(filename+"C2_1-L.dat",'w')
LvsC2_2 = open(filename+"C2_2-L.dat",'w')
LvsC2_3 = open(filename+"C2_3-L.dat",'w')
for op in TcList:
  LvsTc.write(" ".join(op))
  LvsTc.write("\n")
LvsTc.close()
for op in nuinvList:
  Lvsnuinv.write(" ".join(op))
  Lvsnuinv.write("\n")
Lvsnuinv.close()
for op in C2_0List:
  LvsC2_0.write(" ".join(op))
  LvsC2_0.write("\n")
LvsC2_0.close()
for op in C2_1List:
  LvsC2_1.write(" ".join(op))
  LvsC2_1.write("\n")
LvsC2_1.close()
for op in C2_2List:
  LvsC2_2.write(" ".join(op))
  LvsC2_2.write("\n")
LvsC2_2.close()
for op in C2_3List:
  LvsC2_3.write(" ".join(op))
  LvsC2_3.write("\n")
LvsC2_3.close()


