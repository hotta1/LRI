import subprocess

def separate(fname):
  File = open(fname+'.dat','r')
  sizeList = []
  ipList = []
  for line in File:
    lineList = line[:-1].split()
    ipList.append(lineList)
    if sizeList.count(lineList[0]) == 0:
      sizeList.append(lineList[0])
  File.close
  
  print sizeList
  
  pairList = []
  for x in range(0,len(sizeList)-1):
    pairList.append([sizeList[x],sizeList[x+1]])
  
  print pairList
  
  opList = []
  for pair in pairList:
    opList.append(fname+'_'+pair[0])
    opFile = open(fname+'_'+pair[0]+'.dat','a')
    for ipLine in ipList:
      if ipLine[0]==pair[0] or ipLine[0] == pair[1]:
        opFile.write(" ".join(ipLine))
        opFile.write("\n")
    opFile.close()
  return opList

def KHrun(fname):
  Tc = 1.0
  nuinv = 1.0
  zero = 0.0
  KHop = open(fname+".op",'w')
#  subprocess.call(["/home/hotta/KenjiHarada-BSA-efe2dbe/CC2/new_bfss", fname+'.dat', "1", str(Tc), "1", str(nuinv), "1", str(zero)], stdout = KHop)
  subprocess.call(["/home/hotta/KenjiHarada-BSA-efe2dbe/CC2/new_bfss", "-c", fname+'.dat', "1", str(Tc), "1", str(nuinv), "0", str(zero)], stdout = KHop)
#  subprocess.call(["/home/hotta/KenjiHarada-BSA-efe2dbe/CC2/new_bfss", "-c", fname+'.dat', "1", str(Tc), "1", str(nuinv), "1", str(zero), "1", "0", "1", "1", "1", "1"], stdout = KHop)
  KHop.close()
  return


physicalQuantityList = ['binder','combinderinv','combinder12','combinderopt125','comcombinder12inv','combinderinv2dNN']

for physicalQuantity in physicalQuantityList:
  separatedFNList = separate(physicalQuantity)
  print separatedFNList
  LvsTc = open(physicalQuantity+"_Tc-L.dat",'a')
  Lvsnuinv = open(physicalQuantity+"_nuinv-L.dat",'a')
  LvsC2 = open(physicalQuantity+"_C2-L.dat",'a')
  TcList = []
  nuinvList = []
  C2List = []
  for separatedFN in separatedFNList:
    KHrun(separatedFN)
    KHip = open(separatedFN+".op",'r')
    TcLineList = ['','','','']
    nuinvLineList = ['','','','']
    C2LineList = ['','','','']
    for line in KHip:
      ipLineList = line[:-1].split()
      if len(ipLineList)>=4:
        if ipLineList[1]=='LMIN':
          TcLineList[0] = ipLineList[3]
          nuinvLineList[0] = ipLineList[3]
          C2LineList[0] = ipLineList[3]
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


