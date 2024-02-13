# BrugiaTracker
## Overview:
*Brugia malayi* are thread-like parasitic worms and one of the etiological agents of Lymphatic filariasis (LF). Existing anthelmintic drugs to treat LF are effective in reducing the larval microfilaria (mf) counts in human bloodstream but are less effective on adult parasites. To test potential drug candidates, we report a multi-parameter phenotypic assay based on tracking the motility of adult *B. malayi* and mf in vitro. For adult *B. malayi*, motility is characterized by the centroid velocity, path curvature, angular velocity, eccentricity, extent, and Euler Number. These parameters are evaluated in experiments with three anthelmintic drugs. For *B. malayi* mf, motility is extracted from the evolving body skeleton to yield positional data and bending angles at 74 key point. We achieved high-fidelity tracking of complex worm postures (self-occlusions, omega turns, body bending, and reversals) while providing a visual representation of pose estimates and behavioral attributes in both space and time scales.


## Run from the code
### Installation:
Download the MATLAB code files
### Tracking Worm
- Using MATLAB, run 'wellplate_tracking_script.m' in BrugiaAdultTracker or 'worm_tracking_script.m' in BrugiaMFTracker
- After selecting the Video file, the application will ask you to click the worm body, tail, and head.
- The Application automatically tracks the worm body movement and makes the result files, including tracked video and the coordinates xls.
### Analyze Worm Posture
- Using MATLAB, run 'analyze_data.m' in BrugiaAdultTracker\Analysis or 'calculate_all_parameters.m' in BrugiaMFTracker\Analysis
- After selecting the folder containing the results file from the tracker, the application automatically analyzes the parameters and makes the result files.


## Run without MATLAB
### Installation:
Download the Windows executable files from the releases
### Tracking Worm
- Run 'BrugiaAdultTracker.exe' or 'BrugiaMFTracker.exe'
- After selecting the Video file, the application will ask you to click the worm body, tail, and head.
- The Application automatically tracks the worm body movement and makes the result files, including tracked video and the coordinates xls.
### Analyze Worm Posture
- Run 'BrugiaAdultAnalyzer.exe' or'BrugiaMFAnalyzer.exe'
- After selecting the folder containing the results file from the tracker, the application automatically analyzes the parameters and makes the result files.


## Cite as:
Upender Kalwa, Yunsoo Park, Michael J. Kimber and Santosh Pandey, "An automated, high-resolution phenotypic assay for adult *Brugia malayi* and microfilaria," arXiv:2309.03235 [q-bio.QM], Sep. 2023.
https://doi.org/10.48550/arXiv.2309.03235
```
@misc{kalwa2023automated,
      title={An automated, high-resolution phenotypic assay for adult Brugia malayi and microfilaria}, 
      author={Upender Kalwa and Yunsoo Park and Michael J. Kimber and Santosh Pandey},
      year={2023},
      eprint={2309.03235},
      archivePrefix={arXiv},
      primaryClass={q-bio.QM}
}
```
