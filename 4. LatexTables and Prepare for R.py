import os
import ast
from scipy import stats as stats
import numpy as np


Path = "d:/GISOUT/GC5/"
dataFiles = [Path+file for file in os.listdir(Path)]
pltPath = "D:/GISOUT/PLT/"


if not os.path.exists(pltPath):
    os.makedirs(pltPath)

import matplotlib.pyplot as plt
import matplotlib

replaces = [['END','\\\\\\hline']
            ,['-log','-LOG'],['-nb','-NB'],['-dt','-DT'],['-bn','-BN'],['-j48','-J48']
            ,['FIXED-','FX-'],['VAR-','VR-'],['LSHTune','LSH'],
            ['GEN-A','GIS'],['NNFILTER-A','NNF'],['.arff',''],['ALL-TOP-',''],
            ['ant-1.7','ANT'],['camel-1.6','CML'],['ivy-2.0','IVY'],['jedit-4.3','jED'],
            ['log4j-1.2','L4J'],['lucene-2.4','LUC'],['poi-3.0','POI'],['prop-6','PR6'],['synapse-1.2','SYN'],
            ['tomcat-6','TOM'],['velocity-1.6','VEL'],['xalan-2.7','XAL'],['xerces-1.4','XER'],
            ['SUPER-',''],
            ['CLFTU-TUNE-v20-vmx250--','LRNTune-'], ['TUNE-v20-vmx250--','-'], ['TUNEDCLF','LRNTune'],['LRNTune-','Tuned'],
            ['NNF-10','NNF'],['--','-'], 
            ['FX-VNN-GIS','GIS(FX-VNN)'],['VR-VNN-GIS','GIS(VR-VNN)'],['VR-VMUL-GIS','GIS(VR-VMUL)'],['FX-VMUL-GIS','GIS(FX-VMUL)']
            ,['TunedNNF-0','NNF'],['TuneNNF-0','NNF']
           ]


alp = 0.05
def cliffsDelta2(lst1,lst2,
                dull=[0.147, # small
                        0.33,  # medium
                        0.474 # large
                        ][0]): 
    n = gt = lt = 0.0
    for x in lst1:
        for y in lst2:
            n += 1
            if x > y:  gt += 1
            if x < y:  lt += 1
    if n == 0:
        return 0
    return abs(lt - gt) / n  #,abs(lt - gt)/n > dull

def cliffsDelta(lst1,lst2,
                dull=[0.147, # small
                        0.33,  # medium
                        0.474 # large
                        ][0],abst=False): 
    "Returns true if there are more than 'dull' differences"
    m, n = len(lst1), len(lst2)
    lst2 = sorted(lst2)
    j = more = less = 0
    for repeats,x in runs(sorted(lst1)):
        while j <= (n - 1) and lst2[j] < x: 
            j += 1
        more += j * repeats
        while j <= (n - 1) and lst2[j] == x: 
            j += 1
        less += (n - j) * repeats
    d = (more - less) / (m * n) 
    if abst:
        return abs(d) #  > dull
    return d


def runs(lst):
    "Iterator, chunks repeated values"
    for j,two in enumerate(lst):
        if j == 0:
            one,i = two,0
        if one!=two:
            yield j - i,one
            i = j
        one=two
    yield j - i + 1,two


eff = cliffsDelta
testfunc = stats.mannwhitneyu


mtdData = {}
times={}
topSizes = {}
#mes = ['prec', 'rmse','pr0/cnt0','pr1/cnt1','precp/cnt1','AAE','ARE']
mes = ['1rec','1prec','1acc','1F-m']

def doReplaces(data):
    for item in replaces:
        data = data.replace(item[0],item[1])        
    return data


for dataFile in dataFiles:

    file = open(dataFile,'r')
    alines = file.readlines()
    for lindex, line in enumerate(alines):
        
        line = line.replace('\n','')
        if line.startswith('GISCount Train Size:'):
            continue
        elif line.startswith('LSHTune') or line.startswith('TunedNNFILTER') or line.startswith('TuneNNFILTER') or line.startswith('FIXED-') or line.startswith('*FIXED-')  or line.startswith('VAR-') or line.startswith('*VAR-') :
            line = line.split(':') 
            try:
                
                mtd = line[0].strip()
                mtd = doReplaces(mtd)
                ds = line[1].strip()
                perf = line[2][line[2].find('['):line[2].find(']')+1].strip()
                if not mtd in mtdData.keys():
                    mtdData[mtd] = {}
                if not ds in mtdData[mtd].keys():
                    mtdData[mtd][ds] = []
                try:
                    mtdData[mtd][ds].append(ast.literal_eval(perf))
                except Exception as ex:
                    print(ex)
                    print (line, perf)
                    input()
            except Exception as ex2:
                print(ex2)
                print (line)
                input()
        else:

            if line.startswith('#TIME-FOR-IS:'):
                line = line[line.find(':')+1:]
                
                line = line.split(':')
                
                mtd = line[0].strip()
                mtd = doReplaces(mtd)
                ds = line[1].strip()
                t = line[2].strip()
                if not mtd in times.keys():
                    times[mtd] = {}
                if not ds in times[mtd].keys():
                    times[mtd][ds] = []

                times[mtd][ds].append(float(t)/1000.0)
            

print(len(mtdData.keys()))


methods = sorted(list(mtdData.keys()))
lrns = []
for m in methods:
    lrnname = m[m.rfind('-')+1:]
    
    lrns.append(lrnname)

lrns= list(set(lrns))
print (methods)
datasets = sorted(list(mtdData[methods[0]].keys()))


print (methods)
print (datasets)



cnt = len(mtdData[methods[0]][datasets[0]])
header = ','.join([m for m in methods])+'\n'

for item in replaces:
    header = header.replace(item[0],item[1])        



for field in [0,1,-2,-1]:
    for ds in datasets:
        o = open ('OutF/'+mes[field]+'--'+ds+".csv",'w')        
        
        o.write(header)
        for i in range(cnt):
            mvals = []
            for m in methods:
                try:                    
                    mvals.append(str(mtdData[m][ds][i][field]))                    
                except Exception as esx:
                    print (esx)
                    input()
            o.write(','.join(mvals)+'\n')
        o.close()

plt.close('all')
plt.cla()
plt.clf()

testfile = open('tests.txt','w')

plt.figure(figsize=(len(methods),12))
plt.subplots_adjust(wspace=0, hspace=0)

plt.cla()
plt.clf()

for field in [0,1,-2,-1]:



    print (mes[field])
    perfData = []
    mx = 0
    for i in range(len(methods)):
        md = []

        for key in mtdData[methods[i]].keys():
            for perfs in mtdData[methods[i]][key]:
                md.append(perfs[field])
         
        perfData.append(md)
        mx = max(mx,max(md))
    


    s = [np.median(pi) for pi in perfData]
    
    inds = sorted(range(len(s)), key=lambda k: s[k])
    perfData2 = [perfData[ind] for ind in inds]
    lbls = [methods[ind] for ind in inds]

    print ('For filed = ',mes[field])
    testfile.write('For filed = ' + mes[field]+'\n\n\n\n')
    sumefs = []
    sumefscnt = []
    for i in range(len(perfData)):
        sumf = 0
        sumfcnt = 0
        for j in range(len(perfData)):
            if i==j:
                continue
            lst1 = perfData[i]
            lst2 = perfData[j]

            s,p = testfunc(lst1,lst2)
            efs = eff(lst1,lst2)
            if efs>0:
                sumfcnt+=1
            sumf+=efs

            print (methods[i], ' vs. ' , methods[j])
            print (s, p)
            
            if mes[field] == 'ARE' or mes[field] == 'AAE' or mes[field] == 'rmse':
                print ('Effect (*) : ', -1*efs)
            else:
                print ('Effect : ', efs)

            testfile.write(methods[i]+ ' vs. '+ methods[j]+'\n')
            testfile.write (str(s) + '    p-val: '+str(p)+'\n')
            

            if mes[field] == 'ARE' or mes[field] == 'AAE' or mes[field] == 'rmse':
                testfile.write ('Effect *: '+str(-1 * efs)+'\n')
            else:
                testfile.write ('Effect  : '+str(efs)+'\n')



            testfile.write('-------------------------------------\n\n\n')
        sumefs.append(sumf)
        sumefscnt.append(sumfcnt)
    testfile.write('\n\n\n=============================================================\n\n\n')
    
    plt.cla()
    plt.clf()
    #plt.close('all')
    

    
    ax = plt.subplot(111)    
    
    if mx > 5: 
        mx = 5
    if mx < 1:
        mx = 1
    ax.set_ylim([0,mx])
    ax.set_xlim([0,len(methods)+1])
    


    plt.violinplot(perfData2, showmedians = True, showmeans =  True)
    
    
    

    plt.xticks([i + 1 for i in range(len(lbls))])

    pos = 1
    for lbl in lbls:        
        ax.text(pos - 0.5, float(mx) / 2.0, lbl, rotation=90)
        pos+=1

    pos = 1
    for i in range(len(perfData2)):        
        ax.text(pos - 0.2, np.median(perfData2[i]), "M:%.2f"%np.median(perfData2[i]))
        ax.text(pos - 0.2, np.mean(perfData2[i]), "A:%.2f"%np.mean(perfData2[i]))
        pos+=1

    ax.set_ylabel(mes[field])

    plt.savefig('PLT/'+mes[field]+'.png',format = 'png',dpi=200)
    #plt.show()
    
    tblfileM = open('tblMedian-For-'+mes[field]+'--median.txt','w')
    tblfileMt = open('tblMedianAndTime-For-'+mes[field]+'--median.txt','w')
    tblfileAS = open('tblMedian-For-'+mes[field]+'--meanStd.txt','w')
    tblfileALL = open('tblMedian-For-'+mes[field]+'--MedAvgStd.txt','w')
    tblfileM.write('Method & ' +' & '.join(datasets) + '& & END \n')
    tblfileMt.write('Method & ' +' & '.join(datasets) + '& & & END \n')
    tblfileAS.write('Method & ' +' & '.join(datasets) + 'END \n')
    tblfileALL.write('Method & ' +' & '.join(['multicolumn{2}{l|}{'+ds+'}' for ds in datasets])+'END \n')
    tblfileALL.write('Method & ' +' & '.join(['Md & Avg$\\pm$Std' for ds in datasets])+'END \n')
    dsVals = {}
    for ind in inds:
        mtdVals = []
        mtd = methods[ind]
        tblfileM.write(mtd+' & ')
        tblfileMt.write(mtd+' & ')
        tblfileAS.write(mtd+' & ')
        tblfileALL.write(mtd+' & ')
        dssMed = []
        dssMean = []
        dssStd = []
        mtdTimes = []
        for ds in datasets:
            if not ds in dsVals.keys():
                dsVals[ds] = []
            #tblfileM.write(ds+' & ')
            #tblfileAS.write(ds+' & ')
            #tblfileALL.write(ds+' & ')
            dspf = []
            for perfs in mtdData[mtd][ds]:
                dspf.append(perfs[field])
                mtdVals.append(perfs[field])
                dsVals[ds].append(perfs[field])
            dssMed.append(np.median(dspf))
            dssMean.append(np.mean(dspf))
            dssStd.append(np.std(dspf))
            mtdTimes+=times[mtd][ds]

        tblfileM.write(' & '.join([' %.2f '% dsm for dsm in dssMed])+' & %.2f & %.2f$\pm$%.2f END \n'%(np.median(mtdVals),np.mean(mtdVals),np.std(mtdVals)))
        tblfileMt.write(' & '.join([' %.2f '% dsm for dsm in dssMed])+' & %.2f & %.2f$\pm$%.2f & %.1f END \n'%(np.median(mtdVals),np.mean(mtdVals),np.std(mtdVals), np.mean(mtdTimes)))
        tblfileAS.write(' & '.join([' %.2f$\\pm$%.2f '% (dsa,dssStd[i]) for i,dsa in enumerate(dssMean)])+'END \n')
        tblfileALL.write(' & '.join([' %.2f & %.2f$\\pm$%.2f '% (dssMed[i],dsa,dssStd[i]) for i,dsa in enumerate(dssMean)])+'END \n')


    dsMeds = [np.median(dsVals[ds]) for ds in datasets]
    dsMeans = [np.mean(dsVals[ds]) for ds in datasets]
    dsStds = [np.std(dsVals[ds]) for ds in datasets]
    tblfileM.write('\\hline\n Median &' + ' & '.join([' %.2f '% dsm for dsm in dsMeds])+' &  &  END \n')
    tblfileM.write('Average &' + ' & '.join([' %.2f '% dsm for dsm in dsMeans])+' &  &  END \n')
    tblfileM.write('STD &' + ' & '.join([' %.2f '% dsm for dsm in dsStds])+' &  &  END \n')
    tblfileM.close()

    tblfileMt.write('\\hline\n Median &' + ' & '.join([' %.2f '% dsm for dsm in dsMeds])+' &  &  & END \n')
    tblfileMt.write('Average &' + ' & '.join([' %.2f '% dsm for dsm in dsMeans])+' &  &  & END \n')
    tblfileMt.write('STD &' + ' & '.join([' %.2f '% dsm for dsm in dsStds])+' &  &  & END \n')
    tblfileMt.close()

    tblfileAS.close()
    tblfileALL.close()            


    

    for file in os.listdir('.'):
        if file.startswith('tbl') and file.endswith('.txt'):
            o = open(file)
            fdata = o.read()
            o.close()
        
            for item in replaces:
                fdata = fdata.replace(item[0],item[1])        


            o = open(file,'w')
            o.write(fdata)
            o.close()





    
    
    inds = sorted(range(len(sumefs)), key=lambda k: sumefs[k])
    perfData2 = [perfData[ind] for ind in inds]
    lbls = [methods[ind] for ind in inds]


    plt.cla()
    plt.clf()
    #plt.close('all')

    
    ax = plt.subplot(111)    
    
    if mx > 5: 
        mx = 5
    if mx < 1:
        mx = 1
    ax.set_ylim([0,mx])
    ax.set_xlim([0,len(methods)+1])
    


    plt.violinplot(perfData2, showmedians = True, showmeans =  True)
    
    
    

    plt.xticks([i + 1 for i in range(len(lbls))])

    pos = 1
    for lbl in lbls:        
        ax.text(pos - 0.5, float(mx) / 2.0, lbl, rotation=90)
        pos+=1

    pos = 1
    for i in range(len(perfData2)):        
        ax.text(pos - 0.2, np.median(perfData2[i]), "M:%.2f"%np.median(perfData2[i]))
        ax.text(pos - 0.2, np.mean(perfData2[i]), "A:%.2f"%np.mean(perfData2[i]))
        ax.text(pos - 0.2, 0.05, "E:%.2f"%sumefs[inds[i]])
        pos+=1

    ax.set_ylabel(mes[field])


plt.clf()
plt.cla()
plt.close('all')


