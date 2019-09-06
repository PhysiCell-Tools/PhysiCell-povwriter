# PhysiCell POV-writer (Version 1.0.0)
Copyright (c) Paul Macklin 2019, on behalf of PhysiCell & PhysiCell-Tools  
OSI License: BSD-3-Clause (see LICENSE.txt)  

This is an adaptation of Paul Macklin's 3-D cancer immune visualizations.

This processes a series of PhysiCell digital snapshots to create POV files for 3-D raytracing.

Note that users will then need to run POV-Ray on the .pov files, and potentially annotate them (e.g., with ImageMagick) and stitch them together into a movie (e.g., via mencocder). 

Please see http://mathcancer.org/blog/povwriter for a full tutorial. 

## Syntax 
================================================================================  
    povwriter		: run povwriter with config file ./config/settings.xml
    
    povwriter FILENAME.xml	: run povwriter with config file FILENAME.xml
    
    povwriter x:y:z		: run povwriter on data in FOLDER with indices from x 
                   		  to y in incremenets of z
    
                   		 Example: ./povwriter 0:2:10 processes files: 
                   		          ./FOLDER/FILEBASE00000000_physicell_cells.mat
                   		          ./FOLDER/FILEBASE00000002_physicell_cells.mat
                   		          ...
                   		          ./FOLDER/FILEBASE00000010_physicell_cells.mat
                   		 (See the config file to set FOLDER and FILEBASE)
    
    povwriter x1,...,xn	: run povwriter on data in FOLDER with indices x1,...,xn 
    
                   		 Example: ./povwriter 1,3,17 processes files: 
                   		          ./FOLDER/FILEBASE00000001_physicell_cells.mat
                   		          ./FOLDER/FILEBASE00000003_physicell_cells.mat
                   		          ./FOLDER/FILEBASE00000017_physicell_cells.mat
                   		 (Note that there are no spaces.)
                   		 (See the config file to set FOLDER and FILEBASE)
    [/code]                   


