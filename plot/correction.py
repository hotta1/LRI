import pyalps
import matplotlib.pyplot as plt
import pyalps.plot

name = [['|Magnetization|','absM'],['Magnetization^2','magnetization^2']]# name[x][0] is name of physical quantity. name[x][1] is fielname.
name.append(['Magnetization^4','magnetization^4'])
name.append(['Magnetization^6','magnetization^6'])
name.append(['Magnetization^8','magnetization^8'])
name.append(['Energy','energy'])
name.append(['Energy_normalized','energynormalized'])
name.append(['Magnetic Susceptibility','magsus'])
name.append(['Specific Heat','specheat'])
name.append(['Binder Ratio of Magnetization','binder'])
name.append(['Binder Ratio of Magnetization 2','binder2'])
name.append(['Binder Ratio of Magnetization 3','binder3'])
name.append(['Binder Ratio of Magnetization 4','binder4'])
name.append(['Binder Ratio of Magnetization 5','binder5'])
name.append(['Binder Ratio of Magnetization 6','binder6'])
name.append(['Binder Ratio of Magnetization 7','binder7'])
name.append(['Binder Ratio of Magnetization 8','binder8'])
name.append(['Combined Binder 12','combinder12'])
name.append(['Combined Binder 13','combinder13'])
name.append(['Combined Binder 14','combinder14'])
name.append(['Combined Binder 15','combinder15'])
name.append(['Combined Binder 16','combinder16'])
name.append(['Combined Binder 17','combinder17'])
name.append(['Combined Binder 18','combinder18'])
name.append(['Combined Binder 23','combinder23'])
name.append(['Combined Binder 24','combinder24'])
name.append(['Combined Binder 25','combinder25'])
name.append(['Combined Binder 26','combinder26'])
name.append(['Combined Binder 27','combinder27'])
name.append(['Combined Binder 28','combinder28'])
name.append(['Combined Binder 34','combinder34'])
name.append(['Combined Binder 35','combinder35'])
name.append(['Combined Binder 36','combinder36'])
name.append(['Combined Binder 37','combinder37'])
name.append(['Combined Binder 38','combinder38'])
name.append(['Combined Binder 45','combinder45'])
name.append(['Combined Binder 46','combinder46'])
name.append(['Combined Binder 47','combinder47'])
name.append(['Combined Binder 48','combinder48'])
name.append(['Combined Binder 56','combinder56'])
name.append(['Combined Binder 57','combinder57'])
name.append(['Combined Binder 58','combinder58'])
name.append(['Combined Binder 67','combinder67'])
name.append(['Combined Binder 68','combinder68'])
name.append(['Combined Binder 78','combinder78'])
name.append(['Combined Binder Optimized 123','combinderopt123'])
name.append(['Combined Binder Optimized 456','combinderopt456'])
name.append(['Combined Binder Optimized 148','combinderopt148'])
name.append(['Combined Binder Optimized 125','combinderopt125'])
name.append(['Binder Ratio of Magnetization Inverse','binderinv'])
name.append(['Binder Ratio of Magnetization Inverse 2','binderinv2'])
name.append(['Combined Binder Inverse','combinderinv'])
name.append(['Combined Binder Inverse 2','combinderinv2'])
name.append(['Combined Binder Inverse 2dNN','combinderinv2dNN'])
name.append(['Combined Binder Inverse 3dNN','combinderinv3dNN'])
name.append(['Combined Combined Binder 12 Inverse','comcombinder12inv'])
name.append(['Derivative of Combined Binder by K','dercombinderK'])
name.append(['Derivative of Combined Binder by T','dercombinderT'])
name.append(['Derivative of Combined Binder inverse by K','dercombinderinvK'])
name.append(['Derivative of Combined Binder inverse by T','dercombinderinvT'])
#only for cluster picture
name.append(['Magnetization^2 by graph','magnetization^2graph'])
name.append(['Magnetization^4 by graph','magnetization^4graph'])
name.append(['Magnetization^6 by graph','magnetization^6graph'])
name.append(['Magnetization^8 by graph','magnetization^8graph'])
name.append(['Correlation Ratio 1','corratio1'])
name.append(['Correlation Ratio 2','corratio2'])
name.append(['Correlation Ratio 3','corratio3'])
name.append(['CorrelationFunction_16','corfunc_16'])
name.append(['CorrelationFunction_64','corfunc_64'])
name.append(['CorrelationFunction_256','corfunc_256'])
#name.append(['Second Moment','secondmoment'])

for i in range(0,len(name)):
  data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='LRSW_params'),name[i][0])
  for item in pyalps.flatten(data):
    item.props['L'] = int(item.props['L'])
  graph = pyalps.collectXY(data,x='L',y=name[i][0],foreach=['T'])
  f = open(name[i][1]+'.plt','w')
  f.write(pyalps.plot.makeGnuplotPlot(graph))
  f.close()
  print 'finished to output ' + name[i][1] + '.plt'

#import sys

for i in range(0,len(name)):
  fw = open(name[i][1]+'.dat', 'w')
  j=0
  for line in open(name[i][1]+'.plt', 'r'):
    if j<4:
      pass
    elif line=='end \n':
      pass
    else:
      fw.write(line)
    j=j+1
  fw.close()
  print 'finished to output ' + name[i][1] + '.dat'



