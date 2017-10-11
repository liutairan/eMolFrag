from __future__ import print_function
import rdkit
from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Chem import rdmolops
suppl=Chem.MolFromMol2File('/work/tliu7/run0117/CHEMBL281957.mol2',sanitize=False)


def parseMolBlock(molBlock):
    # find V2000
    sdfInfoList=[]
    sdfInfoList = molBlock.split('\n')
    
    fileHead=list(filter(lambda x: 'V2000' in x, sdfInfoList))
    
    fileHeadLineNum=sdfInfoList.index(fileHead[0])
    fileHeadList=fileHead[0].split()
    atomNum=int(fileHead[0][0:3])
    bondNum=int(fileHead[0][3:6])
    
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
    
    return atomInfo,bondInfo

def FindDoubleBonds(inputMol):
    for i in range(inputMol.GetNumBonds()):
        typeValue = inputMol.GetBondWithIdx(i).GetBondTypeAsDouble()
        if abs(typeValue - 2.0) < 0.01:
            startAtom = inputMol.GetBondWithIdx(i).GetBeginAtomIdx()
            endAtom = inputMol.GetBondWithIdx(i).GetEndAtomIdx()
            startSymbol = inputMol.GetAtomWithIdx(startAtom).GetSymbol()
            endSymbol = inputMol.GetAtomWithIdx(endAtom).GetSymbol()
            
            if (startSymbol == 'C' and endSymbol == '*'):
                # meet L7
                return endAtom
            elif (startSymbol == '*' and endSymbol == 'C'):
                # meet L7
                return startAtom
            else:
                pass

    return -1


def atomIndex(parentAtomInfo,atomCoordinate):
    atomX=atomCoordinate[0]
    atomY=atomCoordinate[1]
    atomZ=atomCoordinate[2]
    infoI=parentAtomInfo[0]
    infoX=parentAtomInfo[1]
    infoY=parentAtomInfo[2]
    infoZ=parentAtomInfo[3]
    normList=[]
    for i in range(len(infoI)):
        norm=(atomX-infoX[i])*(atomX-infoX[i])+(atomY-infoY[i])*(atomY-infoY[i])+(atomZ-infoZ[i])*(atomZ-infoZ[i])
        normList.append(norm)
    minInd=normList.index(min(normList))
    index=minInd+1
    return index  #real atom index number, start from 1.


def GetAtomIndexList(parentAtomInfo,sdfAtomInfo):

    atomIndexList=[]
    atomNum=len(sdfAtomInfo[0])
    for i in range(atomNum):
        templist=[sdfAtomInfo[1][i],sdfAtomInfo[2][i],sdfAtomInfo[3][i]]
        tempInd=atomIndex(parentAtomInfo,templist)
        atomIndexList.append(tempInd)
        
    return atomIndexList


def ProcessDoubleBonds(parentMolblock, dbFragList):
    [parentAtomInfo, parentBondInfo] = parseMolBlock(parentMolblock)
    connectedList = []
    atomIndSetAll = []
    connectPointAll = []
    groupIndSetAll = []
    groupSymbolSetAll = []

    if len(dbFragList) >= 2:
        tempFragList = dbFragList
        
        while len(tempFragList) > 0:
            atomIndSet = []
            connectPoint = []
            groupIndSet = []
            groupSymbolSet = []
            tempFrag1 = tempFragList[0]
            [tempFrag1AtomInfo,tempFrag1BondInfo] = parseMolBlock(tempFrag1)
            fragInd1 = GetAtomIndexList(parentAtomInfo, tempFrag1AtomInfo)
            atomIndSet = fragInd1
            groupIndSet.append(fragInd1)
            groupSymbolSet.append(tempFrag1AtomInfo[5])
            tempFragList.remove(tempFrag1)

            restFrags = []
            for frag in tempFragList:
                if frag != tempFrag1:
                    restFrags.append(frag)
            for frag in restFrags:
                [tempFrag2AtomInfo,tempFrag2BondInfo] = parseMolBlock(frag)
                
                fragInd2 = GetAtomIndexList(parentAtomInfo, tempFrag2AtomInfo)
                interSet = list(set(atomIndSet).intersection(fragInd2))
                if len(interSet) >= 2:
                    atomIndSet = list(set(atomIndSet + fragInd2))
                    connectPoint = list(set(connectPoint + interSet))
                    groupIndSet.append(fragInd2)
                    groupSymbolSet.append(tempFrag2AtomInfo[5])
                    tempFragList.remove(frag)
                else:
                    pass
                    
            atomIndSetAll.append(atomIndSet)
            connectPointAll.append(connectPoint)
            groupIndSetAll.append(groupIndSet)
            groupSymbolSetAll.append(groupSymbolSet)
 
        # Get connected index sets print(atomIndSetAll)
        
        for i in range(len(atomIndSetAll)):
            atomIndSet = atomIndSetAll[i]
            atomX = []
            atomY = []
            atomZ = []
            atomA = []
            atomI = []
            atomOL = []
            for atomInd in atomIndSet:
                if atomInd in connectPointAll[i]:
                    # copy info from parent
                    atomI.append(atomInd)
                    atomA.append(parentAtomInfo[5][atomInd-1])
                    atomX.append(parentAtomInfo[1][atomInd-1])
                    atomY.append(parentAtomInfo[2][atomInd-1])
                    atomZ.append(parentAtomInfo[3][atomInd-1])
                    atomOL.append(parentAtomInfo[4][atomInd-1])                
                else:
                    # copy info from parent
                    tempInd = list(filter(lambda x: atomInd in x, groupIndSetAll[i]))
                    tempInd2 = groupIndSetAll[i].index(tempInd[0])
                    tempInd3 = groupIndSetAll[i][tempInd2].index(atomInd)
                    tempChar = groupSymbolSetAll[i][tempInd2][tempInd3]
                    if tempChar == 'R':
                        atomI.append(atomInd)
                        atomA.append('R')
                        atomX.append(parentAtomInfo[1][atomInd-1])
                        atomY.append(parentAtomInfo[2][atomInd-1])
                        atomZ.append(parentAtomInfo[3][atomInd-1])
                        atomOL.append(parentAtomInfo[4][atomInd-1])
                    else:
                        atomI.append(atomInd)
                        atomA.append(parentAtomInfo[5][atomInd-1])
                        atomX.append(parentAtomInfo[1][atomInd-1])
                        atomY.append(parentAtomInfo[2][atomInd-1])
                        atomZ.append(parentAtomInfo[3][atomInd-1])
                        atomOL.append(parentAtomInfo[4][atomInd-1])
                    
            atomInfo = [atomI,atomX,atomY,atomZ,atomOL,atomA]
            
            bondInfo = []
            for bond in parentBondInfo:
                if (int(bond[0]) in atomIndSetAll[i]) and (int(bond[1]) in atomIndSetAll[i]):
                    bondInfo.append(bond[2])
            tempMolblock = GenerateMolblock(atomInfo, bondInfo)
            connectedList.append(tempMolblock)
    else:
        pass

    return connectedList

def GenerateMolblock(atomInfo, bondInfo):
    
    tempMolblockList = []
    tempMolblockList.append('\n')
    tempMolblockList.append('     RDKit          3D\n')
    tempMolblockList.append('\n')
    
    newAtomNum = len(atomInfo[0])
    newBondNum = len(bondInfo)
    newHead=str(newAtomNum).rjust(3)+str(newBondNum).rjust(3)+'  0  0  0  0  0  0  0  0999 V2000\n'
    
    tempMolblockList.append(newHead)
    atomIndMapList = [] # [new] <-> [old]
    atomIndMapList.append(list(range(1,newAtomNum+1)))
    
    dummyIndList = []
    dummyAtomLineList = []
    normalIndList = []
    normalAtomLineList = []

    for i in range(newAtomNum):
        if len(atomInfo[4][i]) > 50:
            if (atomInfo[5][i] == 'R'):
                tempAtomLine = atomInfo[4][i][:31] + 'R ' + atomInfo[4][i][33:]
                dummyIndList.append(atomInfo[0][i])
                dummyAtomLineList.append(tempAtomLine)
            else:
                normalIndList.append(atomInfo[0][i])
                normalAtomLineList.append(atomInfo[4][i])
    
    atomIndMapList.append(normalIndList+dummyIndList)
    newAtomList = normalAtomLineList + dummyAtomLineList

    newBondList = []
    for bond in bondInfo:
        #print(bond)
        tempInd1 = int(bond[0:3])
        tempInd2 = int(bond[3:6])
        tempInfo = bond[6:]
        newInd1 = atomIndMapList[0][atomIndMapList[1].index(tempInd1)]
        newInd2 = atomIndMapList[0][atomIndMapList[1].index(tempInd2)]
        #print(newInd1,newInd2)
        newBond = str(newInd1).rjust(3) + str(newInd2).rjust(3) + tempInfo
        newBondList.append(newBond)
    
    tempMolblockList.append('\n'.join(newAtomList) + '\n')
    tempMolblockList.append('\n'.join(newBondList) + '\n')
    tempMolblockList.append('M  END\n')
    molblockReturn = ''.join(tempMolblockList)
    return molblockReturn


# Input: parent molecule and fragments in mol-object
# Output: fragments in mol-object
def ReconnectDoubleBond(parentMol, inputFrags):
    parentMolblock = Chem.MolToMolBlock(parentMol)
    fragmentMolblocks = []
    for i in range(len(inputFrags)):
        tempFragStr = Chem.MolToMolBlock(inputFrags[i])
        fragmentMolblocks.append(tempFragStr)

    newFragmentMolBlocks = []
    dbFragList = []
    for i in range(len(inputFrags)):
        tempValue = FindDoubleBonds(inputFrags[i])
        if tempValue >= 0:
            # Find C.2 = C.2 bond
            dbFragList.append(fragmentMolblocks[i])
        else:
            newFragmentMolBlocks.append(fragmentMolblocks[i])
    
    reconnectedDBFrags = ProcessDoubleBonds(parentMolblock, dbFragList)
    
    newFragmentMolBlocks = newFragmentMolBlocks + reconnectedDBFrags
    
    newFragmentMol = []
    for i in range(len(newFragmentMolBlocks)):
        tempFragMol = Chem.MolFromMolBlock(newFragmentMolBlocks[i],sanitize=False)
        newFragmentMol.append(tempFragMol)

    # return tuple mol-objects
    return tuple(newFragmentMol)

suppl2 = rdmolops.RemoveHs(suppl)

new2 = BRICS.BreakBRICSBonds(suppl2)
mfl2 = Chem.GetMolFrags(new2, asMols=True, sanitizeFrags=False)

mfl3 = ReconnectDoubleBond(suppl2,mfl2)

for i in range(len(mfl3)):
    print(Chem.MolToMolBlock(mfl3[i]))
    pass
    #FindDoubleBonds(mfl2[i])
