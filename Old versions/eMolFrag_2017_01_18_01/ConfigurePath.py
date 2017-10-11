import os
import sys


def LoadFolder(inputFolderPath):
    infilePathList=[]
    try:
        for root, dirs, files in os.walk(inputFolderPath):
            for file in files:
                infilePathList.append(inputFolderPath+file)
    except:
        print('Error Code: 0007')
    return infilePathList


def CheckExistence(filePath, prop): #folder prop=0, file prop=1
    try:
        if prop==0:
            try:
                if os.path.isdir(filePath):
                    fileList = LoadFolder(filePath)
                    if len(fileList)==0:
                        existFlag=1
                        return existFlag
                    elif filePath+'eMolFrag.py' in fileList:
                        existFlag=0
                        return existFlag
                    else:
                        existFlag=1
                        return existFlag
                else:
                    existFlag=1
                    return existFlag
            except:
                print('Error Code: 0005')
        elif prop==1:
            try:
                if os.path.isfile(filePath):
                    dirPath=os.path.dirname(filePath)
                    if dirPath[-1]=='/':
                        pass
                    else:
                        dirPath=dirPath+'/'
                    #print(dirPath)
                    fileList = LoadFolder(dirPath)
                    #print(len(fileList))
                    if len(fileList)==0:
                        existFlag=1
                        return existFlag
                    elif dirPath+'pkcombu' in fileList:
                        existFlag=0
                        return existFlag
                    else:
                        existFlag=1
                        return existFlag
                else:
                    existFlag=1
                    return existFlag
            except:
                print('Error Code: 0006')
        else:
            pass
        return existFlag  # 0: success, 1: fail
    except:
        print('Error Code: 0004')


def Configure():
    p1=''
    p2=''
    
    try:
        print('Please Type eMolFrag Scripts Directory Absolute Path:')
        flag1=1
        ver = sys.version[0]
        while flag1:
            if ver=='2':
                p1=raw_input()
            elif ver=='3':
                p1=input()
            else:
                print('Get python version failed.')

            if p1[-1] == '/':
                pass
            else:
                p1=p1+'/'
            fb1 = CheckExistence(p1,0)
            if fb1==0:
                flag1=0
            else:
                print('Configure path error, type again:')
    except:
        print('Error Code: 0001')

    try:
        print('Please Type pkcombu Path:')
        flag2=1
        while flag2:
            if ver=='2':
                p2=raw_input()
            elif ver=='3':
                p2=input()
            else:
                print('Get python version failed.')
        
            fb2 = CheckExistence(p2,1)
            if fb2==0:
                flag2=0
            else:
                print('Configure path error, type again:')
    except:
        print('Error Code: 0002')

    try:
        with open(p1+'PathConfigure.log','w') as outf:
            outf.write(p1+'\n'+p2+'\n')
        print('Path Correctly Configured.\nExit.')
    except:
        print('Error Code: 0003')

        
if __name__=="__main__":
    try:
        Configure()
    except:
        print('Error Code: 0000')
