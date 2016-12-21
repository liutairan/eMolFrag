#This script is written by Tairan Liu. 

import sys
import shutil
import os
import os.path
import subprocess
from subprocess import Popen,PIPE
import time


def RmRigidRed(outputPath,inputList):
    pathList=[]
    with open('PathConfigure.log','r') as inf:
        tempList=inf.readlines()
        if len(tempList)==2:
            pathList.append(tempList[0].replace('\n',''))
            pathList.append(tempList[1].replace('\n',''))
        else:
            sys.exit()

    if len(inputList) >1:
        tempInputList=inputList
        while len(tempInputList)>0:
            for mol1 in tempInputList:
                #similar fragments list
                similarList=[]
                alignmentList=[]
                #similar fragments list
                similarList.append(mol1+'\n')
                #get a list of molecules except molA
                restMolList=[]
                for mol2 in inputList:
                    if mol2 != mol1:
                        restMolList.append(mol2)
                aliOutputName=os.path.basename(mol1)+'-alioutput.txt'
                #final result of molA and appendix
                finalMolA=[]
                for molB in restMolList:
            
                    molA=mol1
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

                    if float(tnm) > 0.97:
                    
                        similarList.append(molB+'\n')
                        alignmentList.append(ali)
                        #aliOutputName=os.path.basename(molA)+'-alioutput.txt'
                        #subprocess.call(['python', '/work/tliu7/run0406/script0406/mol-ali-04.py',outputPath, ali, molA, molB, aliOutputName])
                        subprocess.call(['python', pathList[0]+'mol-ali-04.py',outputPath, ali, molA, molB, aliOutputName])
                        #get file output, then add the info of output to a copy of molA
                        molAList=[]
                        with open(molA,'r') as inf:
                            molAList=inf.readlines()
                
                        appendIHead=list(filter(lambda x: '> <BRANCH @atom-number eligible-atmtype-to-connect>' in x, molAList))
                        indAAppendIHead=molAList.index(appendIHead[0])
                        molEnd=list(filter(lambda x: '$$$$' in x, molAList))
                        indAMolEnd=molAList.index(molEnd[0])
                
                        molABasic=molAList[:indAAppendIHead+1]
                
                        finalMolA=molABasic
                
                        

                #read aligned atom number and atom type able to connect, then add to the end of molecule.
                if len(similarList)>1:
                    appdList=[]
                    with open(outputPath+aliOutputName,'r') as inf:
                        appdList=inf.readlines()
                    tempList1=[] #remove empty and '\n' in each element
                    for app in appdList:
                        if len(app)>2:
                            tempPair=app.replace('\n','')
                            tempList1.append(tempPair)
                    tempList2=[] #remove repeat
                    for app in tempList1:
                        if app not in tempList2:
                            tempList2.append(app)
        
                    #group atom types to the same atom number
                    col1=[] #atom number list
                    col2=[] #atom type list,[[],[],[]] each sublist stores the atom types the corresponding atom can connect
                    for app in tempList2:
                        tempApp=app.split()   ###Problem comes from here.
                        if tempApp[0] not in col1: #atom number list
                            col1.append(tempApp[0])
                        tempInd=col1.index(tempApp[0]) #find out the index for each atom number
                        if tempInd+1>len(col2): #if find a new atom number which has not create corresponding atom type list
                            col2.append([])
                        col2[tempInd].append(tempApp[1]) #add new atom type to the sublist corresponding to the atom number
        
                    finalAppdLine=[]
                    for i in range(len(col1)):
                        finalAppdLine.append(' '.join([col1[i]]+col2[i])+'\n')
            
        
                    finalAppdLine.append('\n')
 
                    os.remove(outputPath+aliOutputName)
            
                    #generate final output and move to destination
                    finalMolA=finalMolA+finalAppdLine

                    finalMolA.append('\n> <fragments similar> \n')
                    finalMolA=finalMolA+similarList
                    finalMolA.append('\n')
                    finalMolA.append('$$$$\n')

                    inputFileName=os.path.basename(mol1)
                    outputFilePath=outputPath+'output-rigid/'+inputFileName
                    with open(outputFilePath,'w') as outf:
                        outf.writelines(finalMolA)
                    #finish process molecule molA
                
                    #start process the rest molecules in the similarList
                    for i in range(len(similarList)-1):
                        #the alignment info are stored in the list alignmentList, 1st mol with 2nd mol - alignment[0]; 1st with 3rd - alignment[1]; ...
                        tempMol=similarList[i+1].replace('\n','')
                        tempAli=alignmentList[i]
                        tempAliList=tempAli.split('|')
                        aliList=[]  #[[1,3],[2,4],[3,2],...]
                        for aliLine in tempAliList: 
                            aliList.append(aliLine.split())
                        tempAppd=[]
                        for aliSub in aliList:
                            if aliSub[0] in col1:
                                tempAliInd=col1.index(aliSub[0])
                                tempAliInfoList=[aliSub[1]]+col2[tempAliInd]
                                tempAliInfoLine=' '.join(tempAliInfoList)+'\n'
                                tempAppd.append(tempAliInfoLine)
                
                        #generate final output for this mol
                        tempMolList=[]
                        with open(tempMol,'r') as inf:
                            tempMolList=inf.readlines()

                        appendIHead=list(filter(lambda x: '> <BRANCH @atom-number eligible-atmtype-to-connect>' in x, tempMolList))
                        indAppendIHead=tempMolList.index(appendIHead[0])
                        molEnd=list(filter(lambda x: '$$$$' in x, tempMolList))
                        indMolEnd=tempMolList.index(molEnd[0])

                        tempMolBasic=tempMolList[:indAppendIHead+1]

                        tempFinalMol=tempMolBasic
                        tempFinalMol=tempFinalMol+tempAppd
                

                        tempFinalMol.append('\n> <fragments similar> \n')
                        tempFinalMol=tempFinalMol+similarList
                        tempFinalMol.append('\n')
                        tempFinalMol.append('$$$$\n')

                        inputFileName=os.path.basename(tempMol)
                        outputFilePath=outputPath+'output-rigid/'+inputFileName
                        with open(outputFilePath,'w') as outf:
                            outf.writelines(tempFinalMol)
                        #finish process molecule tempMol
                    #end process the rest molecules in the similarList


                else: # no molecule same to molA or cannot run pkcombu to get result, similar list only contains itself.
                    molA=similarList[0].replace('\n','')
                    molAList=[]
                    with open(molA,'r') as inf:
                        molAList=inf.readlines()

                    appendIHead=list(filter(lambda x: '> <BRANCH @atom-number eligible-atmtype-to-connect>' in x, molAList))
                    indAAppendIHead=molAList.index(appendIHead[0])
                    molEnd=list(filter(lambda x: '$$$$' in x, molAList))
                    indAMolEnd=molAList.index(molEnd[0])
                    finalMolA=molAList[:indAMolEnd]

                    finalMolA.append('\n> <fragments similar> \n')
                    finalMolA=finalMolA+similarList
                    finalMolA.append('\n')
                    finalMolA.append('$$$$\n')
        
                    inputFileName=os.path.basename(mol1)
                    outputFilePath=outputPath+'output-rigid/'+inputFileName
                    with open(outputFilePath,'w') as outf:
                        outf.writelines(finalMolA)
                    #finish process molecule molA

                #remove similarList from tempInputList
                for i in range(len(similarList)):
                    tempInputList.remove(similarList[i].replace('\n',''))

        
                #print rigid-red-out.txt        
                with open(outputPath+'output-log/rigids-red-out.txt','at') as outf:
                    tempSimilarList=list(map(lambda x: x.replace('\n',''), similarList)) #remove '\n' at the end of each molecule in the similarList
                    for i in range(len(tempSimilarList)):
                        outf.write(tempSimilarList[i]+':'+' '.join(tempSimilarList)+'\n')
        
                with open(outputPath+'output-log/rigid-log.txt','at') as outf:
                    outf.write(time.asctime( time.localtime(time.time()) ))
                    outf.write(' ')
                    outf.write(mol1)
                    outf.write(' ')
                    outf.write(str(len(similarList)))
                    outf.write('\n')
                    outf.write('\t'+'\t'.join(similarList))
            


    # only one molecule in the list, there will be no molecule similar to that molecule
    elif len(inputList)==1:
        for mol1 in inputList:
            molAList=[]
            with open(mol1,'r') as inf:
                molAList=inf.readlines()
        
            appendIHead=list(filter(lambda x: '> <BRANCH @atom-number eligible-atmtype-to-connect>' in x, molAList))
            indAAppendIHead=molAList.index(appendIHead[0])
            molEnd=list(filter(lambda x: '$$$$' in x, molAList))
            indAMolEnd=molAList.index(molEnd[0])    
            finalMolA=molAList[:indAMolEnd]
            finalMolA.append('\n> <fragments similar> \n')
            finalMolA.append(mol1+'\n')
            finalMolA.append('$$$$\n')

            inputFileName=os.path.basename(mol1)
            outputFilePath=outputPath+'output-rigid/'+inputFileName
            with open(outputFilePath,'w') as outf:
                outf.writelines(finalMolA)
            with open(outputPath+'output-log/rigids-red-out.txt','at') as outf:
                outf.write(mol1+':'+mol1+'\n')
            with open(outputPath+'output-log/rigid-log.txt','at') as outf:
                outf.write(time.asctime( time.localtime(time.time()) ))
                outf.write(' ')
                outf.write(mol1)
                outf.write(' ')
                outf.write(str(1))
                outf.write('\n')
                outf.write('\t'+mol1+'\n')
    else:
        pass


