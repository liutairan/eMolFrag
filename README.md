# eMolFrag
eMolFrag is a molecular fragmentation tool based on BRICS algorithm written in Python.

This README file is written by Tairan Liu. 
Last Update: 09/13/2016

========================

# Prerequisites:
1. Python 2.7.11 or Python 3.5.2
2. RDKit 2015.09.2 or newer (2016.03.3 has been tested). It is recommended to use Anaconda to install it and use the following command: "conda install -c rdkit rdkit=2015.09.2". Creating conda environment following the instructions on RDKit website may not work properly.
3. pkcombu
4. eMolFrag Scripts 2016_09_09_01
5. (Optional) Openbabel 2.3.1

========================

# Usage:
1. Use ConfigurePath.py to configure paths. Paths only need to be set before the first run if no changes to the paths. After script starts, instructions will be given for setting paths.
    - The first path is for eMolFrag scripts. The absolute path of the scripts folder is needed.
    - The second path is for pkcombu. The absolute path of pkcombu to be used is needed.  
2. Run scripts to process data: "/Path_to_Python/python /Path_to_scripts/eMolFrag.py /Path_to_input_directory/ /Path_to_output_directory/ Number-Of-Cores".
    - /Path_to_Python/python         | Can be simplified as python, ignore the path to it.
    - /Path_to_scripts/eMolFrag.py   | Main entrance to the scripts, relative path is also OK. 
    - /Path_to_input_directory/      | Path of the directory which contains input mol2 files, relative path is also OK.
    - /Path_to_output_directory/     | Path of the directory for output, relative path is also OK.
    - Number-Of-Cores                | Number of processes created in parallel step. It is better to set this parameter no larger than the number of cores of the system/node/cluster.

# Example:
1. mkdir TestEMolFrag/   # Any folder name you want
2. cp /.../test-set100.tar.gz /.../TestEMolFrag/
3. cd /.../TestEMolFrag/
4. tar -zxf test-set100.tar.gz   # Sample input data set
5. python ConfigurePath.py       # Configure path, use absolute path 
    - Step 1: Assume path to eMolFrag.py is /.../eMolFrag_2016_09_09_01/eMolFrag.py, then type: /.../eMolFrag_2016_09_09_01/ at this step.
    - Step 2: Assume path to pkcombu is /.../pkcombu, then type: /.../pkcombu at this step.
6. python /.../eMolFrag_2016_09_09_01/eMolFrag.py /.../TestEMolFrag/test-set100/ /.../TestEMolFrag/outputp-testset100-1/ 16    # Directory name can be changed to whatever you want. 
7. # Check output.

========================

# Output:
1. Find correct output folder, assume /.../Output-xxxx/
2. Then there should be 6 sub folders in this output directory:
   - output-log/        | Log files, some useful temporary files. 
   
      -- InputList               | File contains all the input *.mol2 file names.
      
      -- ListAll                 | File contains all the fragments before reconnect small linkers and total/carbon/nitrogen/oxygen atoms in each fragment.
      
      -- RigidListAll.txt        | File contains all the rigid fragments and total/carbon/nitrogen/oxygen atoms in each fragment.
      
      -- RigidGroupList.txt      | Rigid fragments are grouped by total/carbon/nitrogen/oxygen atoms in each fragment as property, and a number of how many fragments have the same property. 
      
      -- LinkerListAll.txt       | File contains all the linker fragments after reconnect and total/carbon/nitrogen/oxygen atoms in each fragment.
      
      -- LinkerGroupList.txt     | Linker fragments are grouped by total/carbon/nitrogen/oxygen atoms in each fragment as property, and a number of how many fragments have the same property.
      
      -- rigids-red-out.txt      | Rigid fragments and their similar fragments.
      
      -- rigid-log.txt           | Log file for remove redundancy of rigid fragments.
      
      -- linkers-red-out.txt     | Linker fragments and their similar fragments.
		
      -- linker-log.txt          | Log file for remove redundancy of linker fragments.
		
      -- Process.log             | Log file for the whole process.
      
   - output-chop/       | Fragments, rigids and small linkers.
   
   - output-chop-comb/  | Fragments, rigids and large linkers.
   
   - output-sdf/        | Temp files, *.sdf for each input molecule, can be deleted if not use.
   
   - output-rigid/      | Rigid fragments after remove redundancy, with similar fragments in the end of each file.
   
   - output-linker/     | Linker fragments after remove redundancy.

======================== 

# Update Log:
This script is written by Tairan Liu.

   Created       01/17/2016 - Chop

   Modification  01/17/2016 - Remove bug
   
   Modification  01/18/2016 - Reconnect linkers
   
   Modification  01/21/2016 - Remove redundancy
   
   Modification  02/29/2016 - Remove bug
   
   Modification  03/16/2016 - Remove bug
   
   Modification  03/17/2016 - Remove bug
   
   Modification  03/24/2016 - Remove bug
   
   Modification  03/25/2016 - Change each step to functions
   
   Modification  04/03/2016 - Remove bug
   
   Modification  04/06/2016 - Reduce temp output files
   
   Modification  04/06/2016 - Remove bug
   
   Modification  04/06/2016 - Start parallel with chop
   
   Modification  04/17/2016 - Improve efficiency
   
   Modification  04/18/2016 - Remove bug
   
   Modification  05/24/2016 - Start parallel with remove redundancy
   
   Modification  06/14/2016 - Add parallel option as arg
   
   Modification  06/14/2016 - Remove bug
   
   Modification  06/29/2016 - Remove bug
   
   Modification  07/08/2016 - Change similarity criteria of rigids from 0.95 to 0.97
   
   Modification  07/11/2016 - Improve efficiency
   
   Modification  07/18/2016 - Pack up, format.
   
   Modification  07/19/2016 - Solve python 2.x/3.x compatibility problem.
   
   Modification  07/20/2016 - Solve python 2.x/3.x compatibility problem.
   
   Modification  07/21/2016 - Solve python 2.x/3.x compatibility problem.
   
   Modification  07/22/2016 - Solve python 2.x/3.x compatibility problem.
   
   Modification  07/22/2016 - Modify README file

   Last revision 09/13/2016 - Solve output path conflict problem.

