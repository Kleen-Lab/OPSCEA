# OPSCEA for Matlab
================================

Omni-planar and surface casting of epileptiform activity (OPSCEA) (UC Case Number SF2020-281) jointly created by Dr. Jon Kleen, Ben Speidel, Dr. Robert Knowlton, and Dr. Edward Chang is licensed for non-commercial research use at no cost by the Regents of the University of California under CC BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/). Please contact innovation@ucsf.edu if you are interested in using OPSCEA for commercial purposes.  

The following copyright notice and citation is to be included in any publication, material or media wherein all or a part of Licensed Material is contained, “Certain materials incorporated herein are Copyright © 2016 The Regents of the University of California (REGENTS). All Rights Reserved.  

Please cite the following paper in your publications if you have used our software in your research, as well as any relevant toolboxes used herein as appropriate (img_pipe, FreeSurfer): 

Kleen JK, Speidel B, Baud MO, Rao VR, Ammanuel SG, Hamilton LS, Chang EF, Knowlton RC. "Accuracy of omni-planar and surface casting of epileptiform activity for intracranial seizure localization." Epilepsia. 2021 Apr;62(4):947-959. https://doi.org/10.1111/epi.16841
(cover: https://doi.org/10.1111/epi.16080)

The `OPSCEA` MATLAB package. All code included to generate a video of ictal activity projected on a reconstructed brain. Documentation with additional instructions and sample data are provided as well.

## Package setup.
### 1. Required Software

Must have MATLAB (r2019a) or newer installed with Image Processing Toolbox. (note: has not been tested with previous MATLAB versions)

### 2. Other dependencies

No other dependencies.

### 3. Installation
To install the package and have a copy of the code to edit locally, navigate to where you would like to store the package code in your terminal. Clone the package.
```
git clone https://github.com/Kleen-Lab/OPSCEA.git
```
Note that you will have to configure paths in the code to wherever your data is stored. Recommended to store data folder in same location as this repository.

### 4. Download sample data (optional)
Sample data to make practice videos is available for download on OSF in a zip file called OPSCEADATA.zip (https://osf.io/49znp/). 

Once downloaded, unzip file and move it to the same folder as the code from this repo. The paths in the software default to using the sample data in this location.

### 5. Use your own data  
To use your own data with the OPSCEA package, you need to organize it in a specific way. Each patient should have a subdirectory within your data directory, as follows: `PatientName/`. Each patient subdirectory should contain : 
- The `Imaging/` subdirectory, which contains the following:
  - The `MRI/` subdirectory, which contains the MRI data as follows:
    - **brain.mgz.gz** - the t1 MRI scan
    - **aparc+aseg.mgz.gz** - the FreeSurfer full segmentation of the brain
  - The `Meshes/` subdirectory, containing processed surface mesh data as follows:
    - **PatientName_rh_pial.mat** - matlab data file containing the pial cortex mesh data for the right hemisphere
    - **PatientName_lh_pial.mat** - matlab data file containing the pial cortex mesh data for the left hemisphere
    - The `subcortical/` subdirectory, for meshes of the hippocampus and amygdala
      - **rHipp_subcort.mat** matlab data file containing the right hippocampus mesh data
      - **lHipp_subcort.mat** matlab data file containing the left hippocampus mesh data
      - **rAmgd_subcort.mat** matlab data file containing the right amygdala mesh data
      - **lAmgd_subcort.mat** matlab data file containing the left amygdala mesh data
  - The `Elecs/` subdirectory, containing electrode information as follows:
    - **Electrodefile.mat** - matlab data file containing:
      - *elecmatrix* - XYZ coordinates (units: mm; electrode rows by XYZ columns; numeric)
      - *eleclabels* - electrode labels (single column; cell array of strings)
- A subdirectory for each of the seizures to be imaged. These seizures subdirectories should be named `PatientName/PatientName_XX/`, where `XX` is the seizure ID (usually the seizures are numbered starting from 1), and they should contain the following: 
  - **PatientName_XX.mat** - a matlab data file containing:
    - *d* - matrix of clipped peri-ictal seizure ICEEG data (channels by samples), already should be pre-processed (e.g. notch filter, CAR, etc.)
    - *sfx* - sampling frequency (in Hz)
  - **PatientName_XX.badch.mat** - matlab data file containing logical index of bad (artifact-ridden) channels (bad=1, good=0)  
Additionally, the data should be preprocessed in a specific way. These requirements are explained in the [preprocessing](#preprocessing) section.
  
 ## Files of interest  
 The following files are the ones one should be familiar with in order to use OPSCEA:
 - **OPSCEA.m** - the main program
 - **OPSCEAparams.xlsx** - spreadsheet with parameter information for video (i.e. patient, seizure, color axis, size of linelength window, plotting parameters). More details in the [parameters](#parameters) section
 - **litebrain** - simple function that rotated brain into desired position and adds lighting
 - **OPSCEAsurfslice.m** - creates omni-planar slices with rest of brain behind them
 - **toaster.m** - adds heatmap to omni-planar slice
 - **cmOPSCEAcool.m**, **cmOPSCEAjet.m** - colormaps adapted from matlab's 

## Preprocessing  
ICEEG data should be:
- Downsampled to 512Hz
- A notch filter should be applied
- Bad channel indices should be recorded as explained in [Use your own data](#5-use-your-own-data)

## Parameters  
The **OPSCEAparams.xlsx** sheet contains the following parameters:
- *patient* - the `PatientName/`, used to name the patient directory
- *sz* - seizure ID, used to name the seizure subdirectory as `PatientName/PatientName_XX`, where `XX` is the seizure ID
  - If a patient has multiple seizures, each one will occupy a row in the sheet
- *VIDstart* - time in seconds from the start of the ICEEG data file from which the video should start (should be before seizure onset)
- *VIDstop* - time in seconds from the start of the ICEEG data file at which the video should stop
- *BLstart* - time in seconds from the start of the ICEEG data file from which the baseline period should start
- *BLstop* - time in seconds from the start of the ICEEG data file at which the baseline period should stop (should be before seizure onset, otherwise the "baseline" includes seizure data)
- *llw* - length of line length window for ECoG data
- *iceeg_scale* - scale of ECoG gain displayed in ECoG data
- *fps* - frames per seconds of the output video
  - use a very low fps value (e.g. `0.1`) to quickly generate a small video to verify everything is working
- *cax* - values of color axis to control saturation of color on heat map
- *gsp* - gaussian spreading parameter (how far to project the color from each electrode)
- *cm* - color map to determine the spectrum of colors used in the heatmap
  - **cmOPSCEAjet.m** and **cmOPSCEAcool.m** are included
- *iceegwin* - number pf secpmds pf ECoG data to display per frame
- *marg* - margin of data to show preceding the line length window (llw) in ECoG display
- *slicebright* - adjust brightess of the omniplanar slice  

Adjust these parameters for each seizure of each patient at your convenience.  
Additionally, each patient should have a separate sheet, named **PatientName**, in **OPSCEAparams.xlsx**. It should contain the following:
- Accurate depth electrode numbers (from **Electrodefile.m**) and their corresponding desired color in RGB format
- The surface views wanted (e.g. RAI would be a Right Anterior Inferior viewpoint)
- The subplots wanted. ECoG always goes on the right side. 

## Usage
After the data files have been organized as described in [Use your own data](#5-use-your-own-data), the data has been preprocessed as explained in [Preprocessing](#preprocessing), and the parameters have been adjusted as described in [Parameters](#parameters), one can generate a video by running the following command in matlab: 
```
OPSCEA(PatientName, SzId, ShowLabel, JumpTo)
```
Where `PatientName` is the patient ID, `SzId` is the seizure ID to be plotted, `ShowLabel` is a boolean value to decide whether to show the channel labels (1 shows the label; 0 hides the labels and randomizes the order; default is 1), and `JumpTo` is the number of seconds to offset the start of the video by (Default is 0).  
Example usage using the default data, as described in the [sample data](#4-download-sample-data) section:
```
OPSCEA('UCSF1', '01', 1, 0)
```
## References.
-Philipp Berens (2021). Circular Statistics Toolbox (Directional Statistics) (https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics), MATLAB Central File Exchange. Retrieved February 12, 2021.
