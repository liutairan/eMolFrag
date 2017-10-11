#This script is written by Tairan Liu. 

import sys

args=sys.argv

outputPath=args[1]
ali=args[2]
molA=args[3]
molB=args[4]
aliOutputName=args[5]

tempAliList=ali.split('|')
aliList=[]
for aliLine in tempAliList:
    aliList.append(aliLine.split())

molAInfo=[]
with open(molA,'r') as molAIn:
    molAInfo=molAIn.readlines()

appendAHead=list(filter(lambda x: '> <BRANCH @atom-number eligible-atmtype-to-connect>' in x, molAInfo))
indAAppendIHead=molAInfo.index(appendAHead[0])
molAEnd=list(filter(lambda x: '$$$$' in x, molAInfo))
indAMolEnd=molAInfo.index(molAEnd[0])
molAAppd=molAInfo[indAAppendIHead+1:indAMolEnd]

molBInfo=[]
with open(molB,'r') as molBIn:
    molBInfo=molBIn.readlines()

appendBHead=list(filter(lambda x: '> <BRANCH @atom-number eligible-atmtype-to-connect>' in x, molBInfo))
indBAppendIHead=molBInfo.index(appendBHead[0])
molBEnd=list(filter(lambda x: '$$$$' in x, molBInfo))
indBMolEnd=molBInfo.index(molBEnd[0])
molBAppd=molBInfo[indBAppendIHead+1:indBMolEnd]

molAAppdS=[]
for molApp in molAAppd:
    temp1=molApp.replace('\n','')
    if len(temp1) >2:
        temp2=temp1.split()
        molAAppdS.append(temp2)

molBAppdS=[]
for molApp in molBAppd:
    temp1=molApp.replace('\n','')
    if len(temp1) >2:
        temp2=temp1.split()
        molBAppdS.append(temp2)

newAppd=[]
for molApp in molBAppdS:
    tempInd1=molApp[0]
    tempInd2=''
    for aliSub in aliList:
        if aliSub[1]==tempInd1:
            tempInd2=aliSub[0]
    newAppd.append([tempInd2,molApp[1]])
newAppd=newAppd+molAAppdS

tempAppd=[]
for appd in newAppd:
    if appd not in tempAppd:
        tempAppd.append(appd)


#tempAppd=[['1','O.3'],['2','C.3'],['2','C.2'],['2','N.P']]

finalAppd=[]
finalAppdLine=[]
for i in range(len(tempAppd)):
    finalAppdLine.append(' '.join(tempAppd[i]))

with open(outputPath+aliOutputName,'at') as outf:
    tempStr='\n'.join(finalAppdLine)
    
    tempStr=tempStr+'\n'
    outf.writelines(tempStr)




