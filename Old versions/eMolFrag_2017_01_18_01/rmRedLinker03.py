#This script is written by Tairan Liu. 

import sys
import shutil
import os
import os.path
import subprocess
from subprocess import Popen,PIPE
import time


#groupProp: ['T','1','C','1','N','0','O','0']
def RmLinkerRed(outputDir,inputL):
    pathList=[]
    with open('PathConfigure.log','r') as inf:
        tempList=inf.readlines()
        if len(tempList)==2:
            pathList.append(tempList[0].replace('\n',''))
            pathList.append(tempList[1].replace('\n',''))
        else:
            sys.exit()

    inputList=inputL[0]
    groupProp=inputL[1]
    outputPath=outputDir
    outputPath_log=outputDir+'output-log/'
    outputPath_chop_comb=outputDir+'output-chop-comb/'
    outputPath_linker=outputDir+'output-linker/'

    #get a list of input file path without '\n' end
    tempInputList=[]
    for tempinput in inputList:
        tempInputList.append(tempinput.replace('\n',''))


    #start to remove redundancy
    if (groupProp[1]=='1') and (groupProp[3]=='1'): #only one C
        appdList=[]
        pathList=[]
        for mol in tempInputList:
            molList=[]
            with open(mol,'r') as inf:
                molList=inf.readlines()
            appendIHead=list(filter(lambda x: '> <MAX-NUMBER-Of-CONTACTS ATOMTYPES>' in x, molList))
            indAppendIHead=molList.index(appendIHead[0])
            molEnd=list(filter(lambda x: '$$$$' in x, molList))
            indMolEnd=molList.index(molEnd[0])
            molAppd=molList[indAppendIHead+1:indMolEnd]
            molAppList=[]
            for templine in molAppd:
                if len(templine)>2:
                    tempStr=templine.replace('\n','')
                    tempList=tempStr.split()
                    molAppList.append(tempList)
            appdList.append(molAppList)
            pathList.append(mol)
        #find unique structure
        tempAppdList=[]
        for appd in appdList:
            if appd not in tempAppdList:
                tempAppdList.append(appd)
        uniqueFileSampleList=[]
        for appd in tempAppdList:
            ind=appdList.index(appd)
            uniqueFileSampleList.append(pathList[ind])
            #copy file to destination
            molBaseName=os.path.basename(pathList[ind])
            dest=outputPath_linker+molBaseName
            shutil.copyfile(pathList[ind],dest)

            with open(outputPath+'output-log/linkers-red-out.txt','at') as outf:
                outf.write(pathList[ind]+':'+pathList[ind]+':'+'1\n')
                outf.write(pathList[ind]+':'+pathList[ind]+':'+'1.000\n')
            #print(appd[0])
            #write log file
            with open(outputPath+'output-log/linker-log.txt','at') as outf:
                outf.write(time.asctime( time.localtime(time.time()) ))
                outf.write(' ')
                outf.write(pathList[ind])
                outf.write(' ')
                outf.write(str(appdList.count(appd)))
                outf.write('\n')
                outf.write('\tOne Carbon Case - '+' '.join(appd[0])+'\n')


    elif (groupProp[1]=='1') and (groupProp[5]=='1'): #only one N
        appdList=[]
        pathList=[]
        for mol in tempInputList:
            molList=[]
            with open(mol,'r') as inf:
                molList=inf.readlines()
            appendIHead=list(filter(lambda x: '> <MAX-NUMBER-Of-CONTACTS ATOMTYPES>' in x, molList))
            indAppendIHead=molList.index(appendIHead[0])
            molEnd=list(filter(lambda x: '$$$$' in x, molList))
            indMolEnd=molList.index(molEnd[0])
            molAppd=molList[indAppendIHead+1:indMolEnd]
            molAppList=[]
            for templine in molAppd:
                if len(templine)>2:
                    tempStr=templine.replace('\n','')
                    tempList=tempStr.split()
                    molAppList.append(tempList)
            appdList.append(molAppList)
            pathList.append(mol)
        #find unique structure
        tempAppdList=[]
        for appd in appdList:
            if appd not in tempAppdList:
                tempAppdList.append(appd)
        uniqueFileSampleList=[]
        for appd in tempAppdList:
            ind=appdList.index(appd)
            uniqueFileSampleList.append(pathList[ind])
            #copy file to destination
            molBaseName=os.path.basename(pathList[ind])
            dest=outputPath_linker+molBaseName
            shutil.copyfile(pathList[ind],dest)

            with open(outputPath+'output-log/linkers-red-out.txt','at') as outf:
                outf.write(pathList[ind]+':'+pathList[ind]+':'+'1\n')
                outf.write(pathList[ind]+':'+pathList[ind]+':'+'1.000\n')
            #print(appd[0])
            #write log file
            with open(outputPath+'output-log/linker-log.txt','at') as outf:
                outf.write(time.asctime( time.localtime(time.time()) ))
                outf.write(' ')
                outf.write(pathList[ind])
                outf.write(' ')
                outf.write(str(appdList.count(appd)))
                outf.write('\n')
                outf.write('\tOne Nitrogen Case - '+' '.join(appd[0])+'\n')



    elif (groupProp[1]=='1') and (groupProp[7]=='1'): #only one O
        appdList=[]
        pathList=[]
        for mol in tempInputList:
            molList=[]
            with open(mol,'r') as inf:
                molList=inf.readlines()
            appendIHead=list(filter(lambda x: '> <MAX-NUMBER-Of-CONTACTS ATOMTYPES>' in x, molList))
            indAppendIHead=molList.index(appendIHead[0])
            molEnd=list(filter(lambda x: '$$$$' in x, molList))
            indMolEnd=molList.index(molEnd[0])
            molAppd=molList[indAppendIHead+1:indMolEnd]
            molAppList=[]
            for templine in molAppd:
                if len(templine)>2:
                    tempStr=templine.replace('\n','')
                    tempList=tempStr.split()
                    molAppList.append(tempList)
            appdList.append(molAppList)
            pathList.append(mol)
        #find unique structure
        tempAppdList=[]
        for appd in appdList:
            if appd not in tempAppdList:
                tempAppdList.append(appd)
        uniqueFileSampleList=[]
        for appd in tempAppdList:
            ind=appdList.index(appd)
            uniqueFileSampleList.append(pathList[ind])
            #copy file to destination
            molBaseName=os.path.basename(pathList[ind])
            dest=outputPath_linker+molBaseName
            shutil.copyfile(pathList[ind],dest)

            with open(outputPath+'output-log/linkers-red-out.txt','at') as outf:
                outf.write(pathList[ind]+':'+pathList[ind]+':'+'1\n')
                outf.write(pathList[ind]+':'+pathList[ind]+':'+'1.000\n')
            #print(appd[0])
            #write log file
            with open(outputPath+'output-log/linker-log.txt','at') as outf:
                outf.write(time.asctime( time.localtime(time.time()) ))
                outf.write(' ')
                outf.write(pathList[ind])
                outf.write(' ')
                outf.write(str(appdList.count(appd)))
                outf.write('\n')
                outf.write('\tOne Oxygen Case - '+' '.join(appd[0])+'\n')


    else: #all other cases
        if len(tempInputList)>1:
            
            while len(tempInputList)>0:
                for mol1 in tempInputList:
                    molA=mol1
                    molAList=[]
                    with open(molA,'r') as infA:
                        molAList=infA.readlines()
                    appendAIHead=list(filter(lambda x: '> <MAX-NUMBER-Of-CONTACTS ATOMTYPES>' in x, molAList))
                    indAAppendIHead=molAList.index(appendAIHead[0])
                    molAEnd=list(filter(lambda x: '$$$$' in x, molAList))
                    indAMolEnd=molAList.index(molAEnd[0])
                    molAAppd=molAList[indAAppendIHead+1:indAMolEnd]
                    molAAppList=[]
                    for templine in molAAppd:
                        if len(templine)>2:
                            tempStr=templine.replace('\n','')
                            tempList=tempStr.split()
                            molAAppList.append(tempList)

                    similarList=[]
                    similarityList=[]
                    #similar fragments list
                    similarList.append(mol1)
                    similarityList.append("1.000")
                    #get a list of molecules except molA
                    restMolList=[]

                    for mol2 in tempInputList:
                        if mol2 != mol1:
                            restMolList.append(mol2)
                    
                    for molB in restMolList:        
                        try:
                            cmd1=Popen([pathList[1], '-A', molA, '-B', molB, '-oAm'],stdout=PIPE)
                            #cmd1=Popen(['/home/tliu7/apps/pkcombu', '-A', molA, '-B', molB, '-oAm'],stdout=PIPE)
                            #str1=cmd1.communicate()[0]
                            #cmd1.stdout.close()
                            tempstr=(cmd1.stdout.read())
                            str1=tempstr.decode('UTF-8')
                            
                            
                            
                            pos1=str1.index('#   Nmcs|tani|seldis:')
                            
                            str2=str1[pos1:]
                            strList=str2.split('\n')
                            infoLine=strList[1]
            
                            tnm=infoLine.split()[3]
                            ali=infoLine[infoLine.index('|')+1:]
                            
                        except:
                            tnm=str(0.00)
                        
                        if float(tnm)>0.99:
                                   
                            molBList=[]
                            with open(molB,'r') as infB:
                                molBList=infB.readlines()
                            appendBIHead=list(filter(lambda x: '> <MAX-NUMBER-Of-CONTACTS ATOMTYPES>' in x, molBList))
                            indBAppendIHead=molBList.index(appendBIHead[0])
                            molBEnd=list(filter(lambda x: '$$$$' in x, molBList))
                            indBMolEnd=molBList.index(molBEnd[0])
                            molBAppd=molBList[indBAppendIHead+1:indBMolEnd]
                            molBAppList=[]
                            for templine in molBAppd:
                                if len(templine)>2:
                                    tempStr=templine.replace('\n','')
                                    tempList=tempStr.split()
                                    molBAppList.append(tempList)
                
                            aliList=ali.split('|')
                            similarFlag=1
                            for alipair in aliList:
                                aliInd=alipair.split()
                                tempMolBApp=molBAppList[int(aliInd[1])-1]
                                tempMolAApp=molAAppList[int(aliInd[0])-1]
                    
                                if (tempMolAApp[0] == tempMolBApp[0]) and (tempMolAApp[1] == tempMolBApp[1]):
                                    pass
                                else:
                                    similarFlag=0
                            if similarFlag==1:
                                similarList.append(molB) 
                                similarityList.append(tnm)
                    #print(similarList)        
                    #print all the molecules and their similar molecules to the output
                    #for i in range(len(similarList)):
                    with open(outputPath+'output-log/linkers-red-out.txt','at') as outf:
                        outf.write(similarList[0]+':'+similarList[0]+':'+'1\n')
                    for j in range(len(similarList)):
                        with open(outputPath+'output-log/linkers-red-out.txt','at') as outf:
                            outf.write(similarList[0]+':'+similarList[j]+':'+similarityList[j]+'\n')

                    #copy file to destination
                        molBaseName=os.path.basename(mol1)
                        dest=outputPath_linker+molBaseName
                        shutil.copyfile(mol1,dest)

                    #remove all the similar molecules from the list, such that there are less molecules appear in the next loop
                    for i in range(len(similarList)):
                        tempInputList.remove(similarList[i])

                    #write log file
                    with open(outputPath+'output-log/linker-log.txt','at') as outf:
                        outf.write(time.asctime( time.localtime(time.time()) ))
                        outf.write(' ')
                        outf.write(mol1)
                        outf.write(' ')
                        outf.write(str(len(similarList)))
                        outf.write('\n')
                        outf.write('\t'+'\n\t'.join(similarList)+'\n')
                    #print(len(tempInputList))
   
        # only one molecule in the list, there will be no molecule similar to that molecule
        elif len(inputList)==1:
            for mol1 in inputList:
                with open(outputPath+'output-log/linkers-red-out.txt','at') as outf:
                    outf.write(mol1+':'+mol1+':'+'1\n')
                    outf.write(mol1+':'+mol1+':'+'1.000\n')

                #copy file to destination
                molBaseName=os.path.basename(mol1)
                dest=outputPath_linker+molBaseName
                shutil.copyfile(mol1,dest)


                with open(outputPath+'output-log/linker-log.txt','at') as outf:
                    outf.write(time.asctime( time.localtime(time.time()) ))
                    outf.write(' ')
                    outf.write(mol1)
                    outf.write(' ')
                    outf.write(str(1))
                    outf.write('\n')
                    outf.write('\t'+mol1+'\n')
        else:
            pass


