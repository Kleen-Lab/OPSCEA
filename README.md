# OPSCEA for Matlab
================================

Omni-planar and surface casting of epileptiform activity (OPSCEA) (UC Case Number SF2020-281) jointly created by Dr. Jon Kleen, Ben Speidel, Dr. Robert Knowlton, and Dr. Edward Chang is licensed for non-commercial research use at no cost by the Regents of the University of California under CC BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/). Please contact innovation@ucsf.edu if you are interested in using OPSCEA for commercial purposes.  

The following copyright notice and citation is to be included in any publication, material or media wherein all or a part of Licensed Material is contained, “Certain materials incorporated herein are Copyright © 2016 The Regents of the University of California (REGENTS). All Rights Reserved.  

Please cite the following paper in your publications if you have used our software in your research, as well as any relevant toolboxes used herein as appropriate (img_pipe, FreeSurfer): 

Kleen JK, Speidel B, Baud MO, Rao VR, Ammanuel SG, Hamilton LS, Chang EF, Knowlton RC. Accuracy of omni-planar and surface casting of epileptiform activity for intracranial seizure localization. In press at Epilepsia.”  

The `OPSCEA` MATLAB package. All code included to generate a video of ictal activity projected on a reconstructed brain. Documentation with additional instructions and sample data are provided as well.

## Package setup.
### 1. Required Software

Must have MATLAB (r2019a?) or newer installed with Image Processing Toolbox.

### 2. Other dependencies

No other dependencies.

### 3. Installation
To install the package and have a copy of the code to edit locally, navigate to where you would like to store the package code in your terminal. Clone the package.
```
git clone https://github.com/Kleen-Lab/test-opscea3-for-paper.git
```
Note that you will have to configure paths in the code to wherever your data is stored. Paths to sample data (located in 'x' location) are included by default.

## References.
-Philipp Berens (2021). Circular Statistics Toolbox (Directional Statistics) (https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics), MATLAB Central File Exchange. Retrieved February 12, 2021.
