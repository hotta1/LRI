import subprocess
import math

def separate(fname):
  numsizes = 3
  ipFile = open(fname+'.dat','r')
  sizeList = []
  ipList = []
  for line in ipFile:
    lineList = line[:-1].split()
    ipList.append(lineList)
    if sizeList.count(lineList[0]) == 0:
      sizeList.append(lineList[0])
  ipFile.close
  
  print sizeList
  
  pairList = []
  for x in range(0,len(sizeList)-numsizes+1):
    ipPair = []
    for y in range(x,x+numsizes):
      ipPair.append(sizeList[y])
    pairList.append(ipPair)
  
  print pairList
  
  opList = []
  for pair in pairList:
    opElement = []
    opElement.append(fname+'_'+pair[numsizes-numsizes/2-1])
    opElement.append(pair[numsizes-numsizes/2-1])
    opList.append(opElement)
    opFile = open(fname+'_'+pair[numsizes-numsizes/2-1]+'.dat','w')
    for ipLine in ipList:
      for pairelement in pair:
        if ipLine[0]==pairelement:
          opFile.write(" ".join(ipLine))
          opFile.write("\n")
    opFile.close()
  return opList

def KHrun(fname, physicalQuantity):
  Tc = 8.78
  nuinv = 2.5
  gamma = 2.5
  alpha = 1.0
  beta = -1.0
  zero = 0.0
  omega = 1.0
  if fname == physicalQuantity:
    KHop = open(fname+".op",'w')
    KHop2 = open(fname+"_correction.op",'w')
    if physicalQuantity=='magsus' or physicalQuantity=='magsussca':
      subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870/CC2/new_bfss", "-c", fname+'.dat', "1", str(Tc), "1", str(nuinv), "1", str(gamma)], stdout = KHop)
      subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870//CC2/new_bfss", "-f", "1", "-c", fname+'.dat', "1", str(Tc), "1", str(nuinv), "1", str(omega), "1", str(gamma)], stdout = KHop2)
    elif physicalQuantity=='specheat':
      subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870/CC2/new_bfss", "-c", fname+'.dat', "1", str(Tc), "1", str(nuinv), "1", str(alpha)], stdout = KHop)
      subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870/CC2/new_bfss", "-f", "1", "-c", fname+'.dat', "1", str(Tc), "1", str(nuinv), "1", str(omega), "1", str(alpha)], stdout = KHop2)
    elif physicalQuantity=='absM':
      subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870/CC2/new_bfss", "-c", fname+'.dat', "1", str(Tc), "1", str(nuinv), "1", str(beta)], stdout = KHop)
      subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870/CC2/new_bfss", "-f", "1", "-c", fname+'.dat', "1", str(Tc), "1", str(nuinv), "1", str(omega), "1", str(beta)], stdout = KHop2)
    elif physicalQuantity=='binder':
      subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870/CC2/new_bfss", "-c", fname+'.dat', "1", str(Tc), "1", str(nuinv), "0", str(zero)], stdout = KHop)
      subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870/CC2/new_bfss", "-f", "1", "-c", fname+'.dat', "1", str(Tc), "1", str(nuinv), "1", str(omega), "0", str(zero)], stdout = KHop2)
    else:
      subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870/CC2/new_bfss", "-c", fname+'.dat', "1", str(Tc), "1", str(nuinv), "1", str(zero)], stdout = KHop)
      subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870/CC2/new_bfss", "-f", "1", "-c", fname+'.dat', "1", str(Tc), "1", str(nuinv), "1", str(omega), "1", str(zero)], stdout = KHop2)
    KHop.close()
    KHop2.close()
  else:
    KHop = open(fname+".op",'w')
    if physicalQuantity=='magsus' or physicalQuantity=='magsussca':
      subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870/CC2/new_bfss", "-c", fname+'.dat', "1", str(Tc), "1", str(nuinv), "1", str(gamma)], stdout = KHop)
    elif physicalQuantity=='specheat':
      subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870/CC2/new_bfss", "-c", fname+'.dat', "1", str(Tc), "1", str(nuinv), "1", str(alpha)], stdout = KHop)
    elif physicalQuantity=='absM':
      subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870/CC2/new_bfss", "-c", fname+'.dat', "1", str(Tc), "1", str(nuinv), "1", str(beta)], stdout = KHop)
    elif physicalQuantity=='binder':
      subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870/CC2/new_bfss", "-c", fname+'.dat', "1", str(Tc), "1", str(nuinv), "0", str(zero)], stdout = KHop)
    else:
      subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870/CC2/new_bfss", "-c", fname+'.dat', "1", str(Tc), "1", str(nuinv), "1", str(zero)], stdout = KHop)
    KHop.close()
#  subprocess.call(["/home/hotta/command/KenjiHarada-BSA-eb06870/CC2/new_bfss", "-c", fname+'.dat', "1", str(Tc), "1", str(nuinv), "1", str(zero), "1", "0", "1", "1", "1", "1"], stdout = KHop)
  return



physicalQuantityList = ['binder','magsussca','specheat','absM']
#physicalQuantityList = ['binder','magsus']
#physicalQuantityList = ['binder','magsus','absM','M^2']
#physicalQuantityList = ['specheatsca','specheat','absM']
#physicalQuantityList = ['binder','magsus','specheat','absM','M^2']
#physicalQuantityList = ['binder','combinderinv','combinder12','combinderopt125','comcombinder12inv','combinderinv2dNN']
#physicalQuantityList = ['binder']
#physicalQuantityList = ['magsus']
#physicalQuantityList = ['specheat']
#physicalQuantityList = ['specheat','absM']

for physicalQuantity in physicalQuantityList:
  KHrun(physicalQuantity,physicalQuantity)
  separatedFNList = separate(physicalQuantity)
  print separatedFNList
  LvsTc = open(physicalQuantity+"_Tc-L.dat",'w')
  Lvsnuinv = open(physicalQuantity+"_nuinv-L.dat",'w')
  LvsC2 = open(physicalQuantity+"_C2-L.dat",'w')
  TcList = []
  nuinvList = []
  C2List = []
  for separatedFN in separatedFNList:
    KHrun(separatedFN[0], physicalQuantity)
    KHip = open(separatedFN[0]+".op",'r')
    TcLineList = ['','','']
    nuinvLineList = ['','','']
    C2LineList = ['','','']
    TcLineList[0] = separatedFN[1]
    nuinvLineList[0] = separatedFN[1]
    C2LineList[0] = separatedFN[1]
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
            C2LineList[1] = ipLineList[3]
          elif len(ipLineList)==5:
            C2LineList[1] = ipLineList[3]
            C2LineList[2] = ipLineList[4]
    TcList.append(TcLineList)
    nuinvList.append(nuinvLineList)
    C2List.append(C2LineList)
    print TcLineList
  for op in TcList:
    LvsTc.write(" ".join(op))
    LvsTc.write("\n")
  LvsTc.close()
  for op in nuinvList:
    Lvsnuinv.write(" ".join(op))
    Lvsnuinv.write("\n")
  Lvsnuinv.close()
  for op in C2List:
    LvsC2.write(" ".join(op))
    LvsC2.write("\n")
  LvsC2.close()


