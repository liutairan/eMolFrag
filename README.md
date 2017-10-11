# eMolFrag
eMolFrag is a molecular fragmentation tool based on BRICS algorithm written in Python.

This README file is written by Tairan Liu. 
Last Update: 06/19/2017

If you find this tool is useful to you, please cite this paper:

Liu, Tairan, Misagh Naderi, Chris Alvin, Supratik Mukhopadhyay, and Michal Brylinski. "Break Down in Order To Build Up: Decomposing Small Molecules for Fragment-Based Drug Design with e MolFrag." Journal of Chemical Information and Modeling 57, no. 4 (2017): 627-631.

# Prerequisites:
1. Python 2.7.11 or Python 3.5.2
2. RDKit 2015.09.2 or newer (2016.03.3 has been tested). It is recommended to use Anaconda to install it and use the following command: "conda install -c rdkit rdkit=2015.09.2". Creating conda environment following the instructions on RDKit website may not work properly.
3. pkcombu
4. eMolFrag Scripts 201x_xx_xx_xx
5. (Optional) Openbabel 2.3.1


# Usage:
1. Use ConfigurePath.py to configure paths. Paths only need to be set before the first run if no changes to the paths. After script starts, instructions will be given for setting paths.
    - The first path is for eMolFrag scripts. The absolute path of the scripts folder is needed.
    - The second path is for pkcombu. The absolute path of pkcombu to be used is needed.  
2. Run scripts to process data: `/Path_to_Python/python /Path_to_scripts/eMolFrag.py -i /Path_to_input_directory/ -o /Path_to_output_directory/ -p Number-Of-Cores -m Output-selection -c Output-format`.
    - `/Path_to_Python/python`         | Can be simplified as python, ignore the path to it.
    - `/Path_to_scripts/eMolFrag.py`   | Main entrance to the scripts, relative path is also OK. 
    - `/Path_to_input_directory/`      | Path of the directory which contains input mol2 files, relative path is also OK.
    - `/Path_to_output_directory/`     | Path of the directory for output, relative path is also OK.
    - `Number-Of-Cores`                | Number of processes created in parallel step. It is better to set this parameter no larger than the number of cores of the system/node/cluster.
    - `Output-selection`               | Different output select, remove redundancy or not.
    - `Output-format`                  | Keep or the files or put all the output bricks or linkers in 2 or 4 files.
    
|Parameter |  Optional |  Default argument |  Example of argument |  Description|
|:-------------:|:-------------:|:-----:|:-------------:|:-------------:|
|  -i      |      N    |      No default   |      /…/test-set100/    |    Input path | 
|  -o      |      N    |      No default   |      /…/output-100-1/   |    Output path |
|  -p      |      Y    |          1        |            16      |     Parallel cores to be used |
|  -m      |      Y    |          0        |             1      |     Output selection: 0: full process and output; 1: only chop (and reconnect); 2: chop and remove redundancy, but remove temp chop files, only output the rigids and linkers after remove redundancy | 
|  -c      |      Y    |          0        |             1      |     Output format: 1: all linkers in one file, all bricks in one file, all logs in one folder; 2: remove log files; 0: traditional format | 

# Example:
1. `mkdir TestEMolFrag/`   # Any folder name you want
2. `cp /.../test-set100.tar.gz /.../TestEMolFrag/`
3. `cd /.../TestEMolFrag/`
4. `tar -zxf test-set100.tar.gz`   # Sample input data set
5. `python ConfigurePath.py`       # Configure path, use absolute path 
    - Step 1: Assume path to eMolFrag.py is /.../eMolFrag_201x_xx_xx_xx/eMolFrag.py, then type: /.../eMolFrag_201x_xx_xx_xx/ at this step.
    - Step 2: Assume path to pkcombu is /.../pkcombu, then type: /.../pkcombu at this step.
6. `python /.../eMolFrag_201x_xx_xx_xx/eMolFrag.py -i /.../TestEMolFrag/test-set100/ -o /.../TestEMolFrag/outputp-testset100-1/ -p 16 -m 0 -c 0`   # Directory name can be changed to whatever you want. 
7. Check output.


# Output:
1. Find correct output folder, assume /.../Output-xxxx/
2. Then there should be 4 sub folders in this output directory if use "-m" 0 and "-c" 0:
   - `output-log/`        | Log files, some useful temporary files. 
   
      -- `InputList`               | File contains all the input *.mol2 file names.
      
      -- `ListAll`                 | File contains all the fragments before reconnect small linkers and total/carbon/nitrogen/oxygen atoms in each fragment.
      
      -- `BrickListAll.txt`        | File contains all the brick fragments and total/carbon/nitrogen/oxygen atoms in each fragment.
      
      -- `BrickGroupList.txt`      | Brick fragments are grouped by total/carbon/nitrogen/oxygen atoms in each fragment as property, and a number of how many fragments have the same property. 
      
      -- `LinkerListAll.txt`       | File contains all the linker fragments after reconnect and total/carbon/nitrogen/oxygen atoms in each fragment.
      
      -- `LinkerGroupList.txt`     | Linker fragments are grouped by total/carbon/nitrogen/oxygen atoms in each fragment as property, and a number of how many fragments have the same property.
      
      -- `bricks-red-out.txt`      | Brick fragments and their similar fragments.
      
      -- `brick-log.txt`           | Log file for remove redundancy of brick fragments.
      
      -- `linkers-red-out.txt`     | Linker fragments and their similar fragments.
		
      -- `linker-log.txt`          | Log file for remove redundancy of linker fragments.
		
      -- `Process.log`             | Log file for the whole process.
      
   
   - `output-chop-comb/`  | Fragments, bricks and large linkers.
   
   - `output-brick/`      | Brick fragments after remove redundancy.
   
   - `output-linker/`     | Linker fragments after remove redundancy.
 

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

   Modification  09/13/2016 - Solve output path conflict problem.
   
   Modification  12/01/2016 - Change to functions
   
   Modification  12/02/2016 - Error code.
   
   Modification  12/03/2016 - Error code.
   
   Modification  12/04/2016 - Error code.
   
   Modification  12/05/2016 - Error code.
   
   Modification  12/21/2016 - Python 2.x/3.x compatibility problem.
   
   Modification  12/22/2016 - Add new arguments.
   
   Modification  12/23/2016 - Add new arguments.
   
   Modification  12/25/2016 - Remove bug.
   
   Modification  12/26/2016 - Test.

   Modification  12/27/2016 - Minor changes.

   Modification  12/30/2016 - Test.

   Modification  01/12/2017 - Change linker redundancy definition.
   
   Modification  01/13/2017 - Change linker redundancy definition.
   
   Modification  01/18/2017 - Solve double bond problem.

   Last revision 01/18/2017 - Test.
