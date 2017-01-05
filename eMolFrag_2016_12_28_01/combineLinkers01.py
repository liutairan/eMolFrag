#This script is written by Tairan Liu. 

import sys
import os
import time
import os.path
import shutil
import subprocess
from subprocess import Popen,PIPE

import rdkit
from rdkit import Chem

#after chop fragments, find out any linkers and their neighbor are chopped. If two linkers used to connect to each other, then connect them again to get larger linkers.
#input files should be a list of file paths of original molecule and rigids and linkers from that mol.
#input 1 is the output folder path: '/.../output/', the real output folder for the combined linkers is '/.../output/output-chop-comb'.
#input 2 is a list of files, with format: ['/.../CHEMBLxxxxx.mol2', '/.../r-CHEMBLxxxxx.mol2-000.sdf', '/.../l-CHEMBLxxxxx.mol2-000.sdf', ...].  

def parseMol2File(path):
    mol2AllList=[]
    with open(path,'r') as inf:
        mol2AllList=inf.readlines()

    atomHead=mol2AllList.index('@<TRIPOS>ATOM\n')
    bondHead=mol2AllList.index('@<TRIPOS>BOND\n')

    #remove empty line and '\n'
    tempMol2AtomInfo=mol2AllList[atomHead+1:bondHead]
    mol2AtomInfo=[]
    for templine in tempMol2AtomInfo:
        tempstr1=templine.replace('\n','')
        if len(tempstr1)>2:
            mol2AtomInfo.append(tempstr1)

    
    #remove empty line and '\n'
    tempMol2BondInfo=mol2AllList[bondHead+1:]
    mol2BondInfo=[]
    for templine in tempMol2BondInfo:
        tempstr1=templine.replace('\n','')
        if len(tempstr1)>2:
            mol2BondInfo.append(tempstr1)

    #remove mol2 Hs
    (mol2AtomInfo_n,mol2BondInfo_n)=removeMol2Hs(mol2AtomInfo,mol2BondInfo)


    #get coordinates and atom type
    mol2X=[]
    mol2Y=[]
    mol2Z=[]
    mol2A=[]
    mol2I=[]
    for i in range(len(mol2AtomInfo_n)):
        mol2Line=mol2AtomInfo_n[i].split()
        mol2I.append(i+1)
        mol2X.append(float(mol2Line[2]))
        mol2Y.append(float(mol2Line[3]))
        mol2Z.append(float(mol2Line[4]))
        mol2A.append(mol2Line[5])
    atomInfo=[mol2I,mol2X,mol2Y,mol2Z,mol2A]    

    #get bond info list
    bondInfo=[]
    for bond in mol2BondInfo_n:
        templist=bond.split()
        bondInfo.append([templist[1],templist[2]])

    return atomInfo,bondInfo

def parseSDFFile(path):
    #print(path+'\n')
# find V2000
    sdfInfoList=[]
    with open(path,'r') as inf:
        sdfInfoList=inf.readlines()

    fileHead=list(filter(lambda x: 'V2000' in x, sdfInfoList))
    fileHeadLineNum=sdfInfoList.index(fileHead[0])
    fileHeadList=fileHead[0].split()
    atomNum=int(fileHead[0][0:3])
    bondNum=int(fileHead[0][3:6])
    #print('SDF Bond')
    #print(bondNum)
    atomList=sdfInfoList[fileHeadLineNum+1:fileHeadLineNum+atomNum+1]
    bondList=sdfInfoList[fileHeadLineNum+atomNum+1:fileHeadLineNum+atomNum+bondNum+1]

    #get coordinates
    atomX=[]
    atomY=[]
    atomZ=[]
    atomI=[]
    atomOL=[]
    atomA=[]
    atomRemoveList=[]
    atomRemoveInd=[]
    for i in range(len(atomList)):
        atomLine=atomList[i].split()
        #print(atomList[i])
        if atomLine[3]=='H':
            atomRemoveList.append(atomList[i])
            atomRemoveInd.append(str(i+1))
        else:
            atomI.append(i+1)
            atomX.append(float(atomLine[0]))
            atomY.append(float(atomLine[1]))
            atomZ.append(float(atomLine[2]))
            atomOL.append(atomList[i])
            atomA.append(atomLine[3])

    atomInfo=[atomI,atomX,atomY,atomZ,atomOL,atomA]

    #get bond info list
    bondInfo=[]
    
    for bond in bondList:
        #templist=bond.split()
        templist=[bond[0:3],bond[3:6]]+bond[6:].split()
        if (str(int(templist[0])) not in atomRemoveInd) and (str(int(templist[1])) not in atomRemoveInd):
            bondInfo.append([str(int(templist[0])),str(int(templist[1])),bond])
    #print(bondInfo)
    return atomInfo,bondInfo


def removeMol2Hs(atomList,bondList):
    atomRemoveInfoList=[]
    atomRemoveIndList=[]
    for atomLine in atomList:
        templist=atomLine.split()
        prop=templist[5]
        ind=templist[0]
        if prop == 'H':
            atomRemoveInfoList.append(atomLine)
            atomRemoveIndList.append(ind)
    bondRemoveList=[]
    for bondLine in bondList:
        templist=bondLine.split()
        #templist=[bondLine[0:3],bondLine[3:6]]+bondLine[6:].split()
        if (str(int(templist[1])) in atomRemoveIndList) or (str(int(templist[2])) in atomRemoveIndList):
            bondRemoveList.append(bondLine)
    #remove 
    finalAtomList=[]
    for atom in atomRemoveInfoList:
        atomList.remove(atom)
    finalAtomList=atomList

    finalBondList=[]
    for bond in bondRemoveList:
        bondList.remove(bond)
    finalBondList=bondList

    return finalAtomList,finalBondList

def removeSDFHs():
    pass

def removeDummyAtoms():
    pass

def removeAtomAndBond():
    pass

def atomIndex(mol2AtomInfo,atomCoordinate):
    atomX=atomCoordinate[0]
    atomY=atomCoordinate[1]
    atomZ=atomCoordinate[2]
    infoI=mol2AtomInfo[0]
    infoX=mol2AtomInfo[1]
    infoY=mol2AtomInfo[2]
    infoZ=mol2AtomInfo[3]
    normList=[]
    for i in range(len(infoI)):
        norm=(atomX-infoX[i])*(atomX-infoX[i])+(atomY-infoY[i])*(atomY-infoY[i])+(atomZ-infoZ[i])*(atomZ-infoZ[i])
        normList.append(norm)
    minInd=normList.index(min(normList))
    index=minInd+1
    return index  #real atom index number, start from 1.


def GetAtomIndexList(mol2AtomInfo,sdfAtomInfo):

    atomIndexList=[]
    atomNum=len(sdfAtomInfo[0])
    for i in range(atomNum):
        templist=[sdfAtomInfo[1][i],sdfAtomInfo[2][i],sdfAtomInfo[3][i]]
        tempInd=atomIndex(mol2AtomInfo,templist)
        atomIndexList.append(tempInd)
        
    return atomIndexList

def findFragments(outputDir,mol2File,rigidList,linkerList):
    tempSDFName=os.path.basename(mol2File)+'.sdf'
    tempSDFPath=outputDir+'output-sdf/'+tempSDFName
    #print('SDF')
    sdfInfo=parseSDFFile(tempSDFPath)

    outputPath_chop_comb=outputDir+'output-chop-comb/'
    outputPath_log=outputDir+'output-log/'

    mol2Info=parseMol2File(mol2File)
    rigidAtomList=[]
    for rigidFile in rigidList:
        rigidBaseName=os.path.basename(rigidFile)
        destPath=outputPath_chop_comb+rigidBaseName
        shutil.copyfile(rigidFile,destPath)

        tempInfo=parseSDFFile(rigidFile)
        #print(rigidFile)
        tempAtomIndexList=GetAtomIndexList(mol2Info[0],tempInfo[0])
        rigidAtomList.append(tempAtomIndexList)
        atomCount=len(tempInfo[0][0])
        carbonCount=0
        nitrogCount=0
        oxygenCount=0
        for atom in tempInfo[0][5]:
            if atom == 'C':
                carbonCount=carbonCount+1
            elif atom == 'N':
                nitrogCount=nitrogCount+1
            elif atom == 'O':
                oxygenCount=oxygenCount+1
            else:
                pass
        countList=[destPath,atomCount,carbonCount,nitrogCount,oxygenCount]
        writeRLList(outputPath_log+'RigidListAll.txt',countList)


    linkerAtomList=[]
    for linkerFile in linkerList:
        tempInfo=parseSDFFile(linkerFile)
        tempAtomIndexList=GetAtomIndexList(mol2Info[0],tempInfo[0])
        linkerAtomList.append(tempAtomIndexList)

    rigidAtomAll=[]
    for i in range(len(rigidAtomList)):
        rigidAtomAll=rigidAtomAll+rigidAtomList[i]
    rigidAtomAll.sort()

    linkerAtomAll=[]
    for i in range(len(linkerAtomList)):
        linkerAtomAll=linkerAtomAll+linkerAtomList[i]
    linkerAtomAll.sort()

    newLinkerList=[]
    fragmentsList=[] #return value, may be empty
    fragmentsCountList=[] #[T,C,N,O]

    if len(linkerAtomList)>0:
        if len(linkerAtomList[0])>0:
            while len(linkerAtomAll)>0:
                templist1=list(filter(lambda x: linkerAtomAll[0] in x, linkerAtomList))
                templist2=templist1[0]  # atoms of the linker contain the first atom of linkerAtomAll

                conFlag=1
                templistlen=len(templist2)
                while conFlag:
                    for atom in templist2:
                        templist3=list(filter(lambda x: str(atom) in x, mol2Info[1])) # [[],[],[]]

                        for tempconnection in templist3: #['1','24'] in [['1','24']]
                            if (int(tempconnection[0]) not in templist2) and (int(tempconnection[0]) not in rigidAtomAll):
                                templist4=list(filter(lambda x: int(tempconnection[0]) in x, linkerAtomList))
                                templist5=templist4[0]
                                templist2=templist2+templist5
                            elif (int(tempconnection[1]) not in templist2) and (int(tempconnection[1]) not in rigidAtomAll):
                                templist4=list(filter(lambda x: int(tempconnection[1]) in x, linkerAtomList))
                                templist5=templist4[0]
                                templist2=templist2+templist5
                            else:
                                pass
    
                    templist6=[]
                    for tempatom in templist2:
                        if tempatom not in templist6:
                            templist6.append(tempatom)
                    templist2=templist6 #remove repeat atom index
                    if len(templist2)>templistlen: #connected atom number still increase, means still haven't find all the connected linkers
                        templistlen=len(templist2)
                    else: #connected atom number not change any more
                        conFlag=0
                #this new linker is found
                newLinkerList.append(sorted(templist2))
                for tempAtom in templist2:
                    linkerAtomAll.remove(tempAtom)
    
        else: # the first linker does not have atom 
            pass
    else: # no linker exist
        pass
    
    tempFileNameInd=0
    for linker in newLinkerList:
        tempStore1=[] #atom info
        tempStore2=[] #bond info
        tempStore3=[] #appendix info
        for atomIndex in linker:
            tempStore1.append(sdfInfo[0][4][atomIndex-1])

        #write atom count info to LinkerListAll.txt
        atomCount=len(tempStore1)
        carbonCount=0
        nitrogCount=0
        oxygenCount=0
        for atomline in tempStore1:
            templist=atomline.split()
            atomprop=templist[3]
            if atomprop =='C':
                carbonCount=carbonCount+1
            elif atomprop =='N':
                nitrogCount=nitrogCount+1
            elif atomprop =='O':
                oxygenCount=oxygenCount+1
            else:
                pass
        fragmentsCountList.append([atomCount,carbonCount,nitrogCount,oxygenCount])


        for bond in sdfInfo[1]:  #bond: ['1','3','  1  3  1  0']
            if (int(bond[0]) in linker) and (int(bond[1]) in linker):
                tempstr1=str(linker.index(int(bond[0]))+1)
                tempstr2=str(linker.index(int(bond[1]))+1)
                tempstr3=tempstr1.rjust(3)+tempstr2.rjust(3)+bond[2][6:]
                tempStore2.append(tempstr3)

        for atomIndex in linker:
            templist=list(filter(lambda x: str(atomIndex) in x,sdfInfo[1]))
            tempcount=len(templist)
            for templist2 in templist:
                if (int(templist2[0]) in linker) and (int(templist2[1]) in linker):
                    tempcount=tempcount-1
    
            tempStore3.append(str(tempcount)+' '+mol2Info[0][4][atomIndex-1]+'\n')


        newAtomNum=len(tempStore1)
        newBondNum=len(tempStore2)
        newHead=str(newAtomNum).rjust(3)+str(newBondNum).rjust(3)+'  0  0  0  0  0  0  0  0999 V2000\n'

        baseFileName=os.path.basename(mol2File)
        tempFileName='l-'+baseFileName+'-'+str(tempFileNameInd).zfill(3)+'.sdf'
        tempFileNameInd=tempFileNameInd+1        

        tempFileDataList=[]
        tempFileDataList.append(tempFileName+'\n')
        tempFileDataList.append('     RDKit          3D\n')
        tempFileDataList.append('\n')
        tempFileDataList.append(newHead)
        tempFileDataList=tempFileDataList+tempStore1
        tempFileDataList=tempFileDataList+tempStore2
        tempFileDataList.append('M  END\n')
        tempFileDataList.append('\n')
        tempFileDataList.append('> <MAX-NUMBER-Of-CONTACTS ATOMTYPES> \n')
        tempFileDataList=tempFileDataList+tempStore3
        tempFileDataList.append('\n')
        tempFileDataList.append('$$$$\n')

        fragmentsList.append(tempFileDataList)
    return fragmentsList,fragmentsCountList



def writeLinkers():
    pass

def writeSDFFile(path,content):
    with open(path,'w') as outf:
        #tempstr='\n'.join(content)
        outf.writelines(content)

def writeRLList(path,list):
    #RL List is a list of file path and atom count of total atom number and C, N, O numbers. It is used for the next step to do the group process and remove redundancy.
    with open(path,'at') as outf:
        tempStr=list[0]+' T '+str(list[1])+' C '+str(list[2])+' N '+str(list[3])+' O '+str(list[4])+'\n'
        outf.write(tempStr)
    pass 


def combineLinkers(outputDir,inputFileList):
    outputFolderPath_log=outputDir+'output-log/'
    outputFolderPath_chop=outputDir+'output-chop/'
    outputFolderPath_chop_comb=outputDir+'output-chop-comb/'
    outputFolderPath_sdf=outputDir+'output-sdf/'

    outputFolderPath_chop_comb=outputDir+'output-chop-comb/'
    if not os.path.exists(outputFolderPath_chop_comb):
        os.mkdir(outputFolderPath_chop_comb)

    outputFolderPath_sdf=outputDir+'output-sdf/'
    if not os.path.exists(outputFolderPath_sdf):
        os.mkdir(outputFolderPath_sdf)


    if len(inputFileList)>0: #not empty list
        originalFile=inputFileList[0]
        rigidList=[]
        linkerList=[]
        if len(inputFileList)>1: # rigids or linkers exist
            for i in range(len(inputFileList)-1):
                tempStr1=inputFileList[i+1]
                tempStr2=os.path.basename(tempStr1)
                if tempStr2[0]=='r':
                    rigidList.append(tempStr1)
                elif tempStr2[0]=='l':
                    linkerList.append(tempStr1)
                else:
                    pass
    
            #find fragments
            (fragmentsList,fragmentsCountList)=findFragments(outputDir,originalFile,rigidList,linkerList)
    
            #write linkers to file
            baseFileName=os.path.basename(originalFile) # base name, eg: xxx.mol2
            for i in range(len(fragmentsList)):
                fragment=fragmentsList[i]
                tempFileName='l-'+baseFileName+'-'+str(i).zfill(3)+'.sdf'
                writeSDFFile(outputFolderPath_chop_comb+tempFileName,fragment)
                fragmentCount=fragmentsCountList[i]

                #tempStr=tempFileName+' T '+str(fragmentCount[0])+' C '+str(fragmentCount[1])+' N '+str(fragmentCount[2])+' O '+str(fragmentCount[3])+'\n'
                tempList=[outputFolderPath_chop_comb+tempFileName]+fragmentCount
                writeRLList(outputFolderPath_log+'LinkerListAll.txt',tempList)




