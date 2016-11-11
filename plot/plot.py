import pyalps
import matplotlib.pyplot as plt
import pyalps.plot

name = [['|Magnetization|','absM'],['Magnetization^2','magnetization^2']]# name[x][0] is name of physical quantity. name[x][1] is fielname.
name.append(['Magnetization^4','magnetization^4'])
name.append(['Magnetic Susceptibility','magsus'])
name.append(['Magnetic Susceptibility for Scaling','magsussca'])
name.append(['Binder Ratio of Magnetization','binder'])
name.append(['Energy','energy'])
name.append(['Energy_normalized','energynormalized'])
name.append(['Specific Heat','specheat'])
name.append(['Magnetic Susceptibility disconnected','magsusdis'])
name.append(['Magnetic Susceptibility for Scaling disconnected','magsusscadis'])
name.append(['Binder Ratio of Magnetization connected','bindercon'])
name.append(['Binder Ratio of Magnetization disconnected','binderdis'])
name.append(['Specific Heat disconnected','specheatdis'])
name.append(['Specific Heat Conventional','specheatconv'])
name.append(['Specific Heat by FT','specheatft'])
#name.append(['Specific Heat for Scaling','specheatsca'])
#name.append(['Combined Binder 12','combinder12'])
#name.append(['Combined Binder Inverse','combinderinv'])
#name.append(['Combined Binder Inverse 2dNN','combinderinv2dNN'])
#name.append(['Combined Binder Inverse 3dNN','combinderinv3dNN'])
#only for cluster picture
#name.append(['Correlation Ratio 1','corratio1'])
#name.append(['Correlation Ratio 2','corratio2'])
#name.append(['Correlation Ratio 3','corratio3'])
#name.append(['Correlation_16','corfunc_16'])
#name.append(['Correlation_64','corfunc_64'])
#name.append(['Correlation_256','corfunc_256'])
#name.append(['Second Moment','secondmoment'])

for i in range(0,len(name)):
  data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='LRSW_params'),name[i][0])
  for item in pyalps.flatten(data):
    item.props['L'] = int(item.props['L'])
  graph = pyalps.collectXY(data,x='T',y=name[i][0],foreach=['L'])
  graph.sort(key=lambda item: item.props['L'])
  f1 = open(name[i][1]+'.plt','w')
  f1.write(pyalps.plot.makeGnuplotPlot(graph))
  f1.close()
  f2 = open(name[i][1]+'.dat','w')
  for j in graph:
    L=j.props['L']
    for k in range(0,len(j.x)):
      f2.write(str(L)+' '+str(j.x[k])+' '+str(j.y[k].mean)+' '+str(j.y[k].error)+'\n')
  f2.close()
  print 'finished to output ' + name[i][1] + '.plt and ' + name[i][1] + '.dat'


"""
for i in range(0,len(name)):
  data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='LRSW_params'),name[i][0])
  graph = pyalps.collectXY(data,x='T',y=name[i][0],foreach=['L'])
  f = open(name[i][1]+'.plt','w')
  f.write(pyalps.plot.makeGnuplotPlot(graph))
  f.close()

for i in range(0,len(name)):
  data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='LRSW_params'),name[i][0])
  graph = pyalps.collectXY(data,x='T',y=name[i][0],foreach=['L'])
  f = open(name[i][1]+'.dat','w')
  for j in graph:
    L=j.props['L']
    for k in range(0,len(j.x)):
      f.write(str(L)+' '+str(j.x[k])+' '+str(j.y[k].mean)+' '+str(j.y[k].error)+'\n')
  f.close()

#import sys
for i in range(0,len(name)):
  fw = open(name[i][1]+'.dat', 'w')
  j=0
  for line in open(name[i][1]+'.plt', 'r'):
    if j<5:
      pass
    elif line=='end \n':
      pass
    else:
      fw.write(line)
    j=j+1
  fw.close()
"""

