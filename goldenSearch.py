def myFunction(Youngs,Jobname,resultfile):
   from part import *
   from material import *
   from section import *
   from assembly import *
   from step import *
   from interaction import *
   from load import *
   from mesh import *
   from job import *
   from sketch import *
   from visualization import *
   from connectorBehavior import *
   from math import sqrt, pow
   mdb=openMdb('FemurModelV4.3.cae')
   mdb.models['GoldenSearchModelNL'].materials['NewMaterial'].elastic.setValues(table=((Youngs,0.3),))
   myJob = mdb.Job(name = str(Jobname), model='GoldenSearchModelNL', multiprocessingMode=DEFAULT, numCpus=4,    numDomains=4)
   myJob.submit()
   myJob.waitForCompletion()
   ##--- Odb business---- -------
   Outputfile = Jobname + '.odb'
   odb=openOdb( path=str(Outputfile))
   BottomIFM = odb.rootAssembly.instances['Cortical Inferior-1'].nodeSets['IFM']
   TopIFM = odb.rootAssembly.instances['Cortical Superior-1'].nodeSets['IFM']

   TopOriginalU = odb.steps['Step-5'].frames[-1].fieldOutputs['U'].getSubset(region=TopIFM)
   TopFinalU = odb.steps['Step-6'].frames[-1].fieldOutputs['U'].getSubset(region=TopIFM)
   TopOriginal = odb.steps['Step-6'].frames[-1].fieldOutputs['COORD'].getSubset(region=TopIFM)
   
   BottomOriginalU = odb.steps['Step-5'].frames[-1].fieldOutputs['U'].getSubset(region=BottomIFM)
   BottomFinalU = odb.steps['Step-6'].frames[-1].fieldOutputs['U'].getSubset(region=BottomIFM)
   BottomOriginal = odb.steps['Step-6'].frames[-1].fieldOutputs['COORD'].getSubset(region=BottomIFM)

   zprimeTop = (TopOriginal.values[0].data[0]+TopOriginalU.values[0].data[0])*sin(radians(6)) + (TopOriginal.values[0].data[2]+TopOriginalU.values[0].data[2])*cos(radians(6))
   zprimeBottom = (BottomOriginal.values[0].data[0]+BottomOriginalU.values[0].data[0])*sin(radians(6)) + (BottomOriginal.values[0].data[2]+BottomOriginalU.values[0].data[2])*cos(radians(6))
   zprime1deformed = (TopOriginal.values[0].data[0]+TopFinalU.values[0].data[0])*sin(radians(6)) + (TopOriginal.values[0].data[2]+TopFinalU.values[0].data[2])*cos(radians(6))
   zprime2deformed = (BottomOriginal.values[0].data[0]+BottomFinalU.values[0].data[0])*sin(radians(6)) + (BottomOriginal.values[0].data[2]+BottomFinalU.values[0].data[2])*cos(radians(6))

   Lastframe = odb.steps['Step-6'].frames[-1]
   print zprimeTop
   print zprime1deformed
   print zprimeBottom
   print zprime2deformed
   IFM = (zprimeTop-zprimeBottom)-(zprime1deformed-zprime2deformed)
   file1 = open(resultfile,'a')
   file1.write("%5.8f %5.8f\n" % (Youngs,IFM))
   file1.write("%s\n" %(Lastframe))
   file1.close()
   fx = abs(0.2-abs(abs(zprimeTop-zprimeBottom)-abs(zprime1deformed-zprime2deformed)))
   odb.close()  
   print fx
   return fx
   ## get ur stress corresponding to the specified point and return it

  

def search(f,a,b, resultfile,tol=0.1, jobname = 'NLSteelRight'):
   from math import log, sqrt
   nIter = -2.078087*log( tol/abs(b- a))
   #nIter = 5
   R = 0.618033989
   C = 1.0 - R
   # First telescoping
   x1 = R*a + C*b; x2 = C*a + R*b
   f1 = f(x1,jobname + (str(round(x1,3))).replace('.','-'),resultfile); f2 = f(x2,jobname + (str(round(x2,3))).replace('.','-'),resultfile)
   print a, b, x1, x2
   file1 = open(resultfile,'a')
   file1.write("%5.8f %5.8f %5.8f %5.8f\n" % (a, b, x1, x2,))
   file1.close()

   # Main loop
   for i in range(int(nIter) ):
       if f1 < f2:
		   b = x2;
		   x2 = x1;f2 = f1;
		   x1 = R*a +C*b;f1 = f(x1,jobname+(str(round(x1,3))).replace('.','-'),resultfile)
		   print a,b,x1,x2
		   file1 = open(resultfile,'a')
		   file1.write("%5.8f %5.8f %5.8f %5.8f\n" % (a, b, x1, x2,))
		   file1.close()
       else:
     	   a = x1;
           x1 = x2;f1 = f2;
           x2 = C*a + R*b;f2 = f(x2,jobname+(str(round(x2,3))).replace('.','-'),resultfile)
           print a,b,x1,x2
           file1 = open(resultfile,'a')
           file1.write("%5.8f %5.8f %5.8f %5.8f\n" % (a, b, x1, x2,))
           file1.close()
   else: return x2,f2
   

resultfile = 'NLSteelRight.txt'
file1 = open(resultfile,'a')
file1.write("START\n")
file1.close()
x,fMin = search(myFunction,300,800,resultfile)
file1 = open(resultfile,'a')
file1.write("DONE\n")
file1.close()

