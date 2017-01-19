import os
import sys
import pkgutil

def LoadFolder(inputFolderPath):
    infilePathList=[]
    try:
        for root, dirs, files in os.walk(inputFolderPath):
            for file in files:
                infilePathList.append(inputFolderPath+file)
    except:
        print('Error Code: 0015')
    return infilePathList

def Loader(mainPath0):
    flag=0
    
    try:
        mainEntryPath=mainPath0
        mainPath=os.path.dirname(mainEntryPath)
        if mainPath[-1]=='/':
            pass
        else:
            mainPath=mainPath+'/'

        fileList=LoadFolder(mainPath)
    except:
        print('Error Code: 0010')
    
    # Check path configure file
    try:
        if mainPath+'PathConfigure.log' in fileList:
            with open(mainPath+'PathConfigure.log','r') as inf:
                tempList=inf.readlines()
                if len(tempList)==2:
                    if os.path.exists(tempList[0].replace('\n','')):
                        if os.path.exists(tempList[1].replace('\n','')):
                            pass
                        else:
                            flag = 1
                            print('Pkcombu path is not configured correctly.\nExit.')
                            return 1
                    else:
                        flag = 1
                        print('eMolFrag path is not configured correctly.\nExit.')
                        return 1
                else:
                    flag = 1
                    print('Path configuration is not correctly.\nExit.')
                    return 1
        else:
            flag = 1
            print('Cannot find path configure file.\nExit.')
            return 1
    except:
        print('Error Code: 0011')

    # Check main entry eMolFrag.py
    try:
        if mainEntryPath in fileList:
            pass
        else:
            flag = 1
            print('Cannot find entrance: eMolFrag.py.\nExit.')
            return 1
    except:
        print('Error Code: 0012')

    # Check RDKit
    try:
        load1=pkgutil.find_loader('rdkit')
        #load2=pkgutil.find_loader('rdkit.Chem')
        if load1 == None:
            print('Cannot find RDKit.\nExit.')
            flag = 1
            return 1
        else:
            load2=pkgutil.find_loader('rdkit.Chem')
            if load2 == None:
                print('RDKit is not properly installed.\nExit.')
                flag = 1
                return 1
            else:
                pass
    except:
        print('Error Code: 0013')


    # Check pkcombu
    # Already checked previously.

    # Check function lib
    try:
        if os.path.exists(mainPath+'chopRDKit03.py'):
            pass
        else:
            flag = 1
            print('Cannot find part of script files.\nExit.')
            return 1

        if os.path.exists(mainPath+'combineLinkers01.py'):
            pass
        else:
            flag = 1
            print('Cannot find part of script files.\nExit.')
            return 1

        if os.path.exists(mainPath+'mol-ali-04.py'):
            pass
        else:
            flag = 1
            print('Cannot find part of script files.\nExit.')
            return 1

        if os.path.exists(mainPath+'rmRedLinker04.py'):
            pass
        else:
            flag = 1
            print('Cannot find part of script files.\nExit.')
            return 1

        if os.path.exists(mainPath+'rmRedBrick01.py'):
            pass
        else:
            flag = 1
            print('Cannot find part of script files.\nExit.')
            return 1
    except:
        print('Error Code: 0014')
    
    return flag
