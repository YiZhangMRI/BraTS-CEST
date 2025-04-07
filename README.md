# BraTS-CEST
MATLAB code for BraTS-CEST data synthesis pipeline
![banner](https://github.com/8fanmao/BraTS-CEST/blob/main/banner.jpg)

## Get Started
### Use GUI (recommended) 
- Click `GUI_BraTS_CEST.mlapp` 
- Follow the tutorial `GUI_tutorial.pdf`
  
### Use script 
Run the following scripts in sequence to generate the data
- `step1_fastMRI_h5_to_mat.m`
- `step2_BraTS_r2c.m`
- `step3_BraTS_CEST.m`
  
Run the following script to view results or quality control
- `stepQC_BraTS_cpx.m`
- `stepQC_BraTS_CEST.m`

## Requirements
`Matlab R2018a` or above is required

Please make sure you have the following official toolboxs:
- `Image Processing Toolbox`
- `Curve Fitting Toolbox`
- `Partial Differential Equation Toolbox`
  
Please install `Parallel Computing Toolbox` for better performance

## Acknowledgments
[1] Menze BH, Jakab A, Bauer S, Kalpathy-Cramer J, Farahani K, Kirby J, Burren Y, Porz N, Slotboom J, Wiest R. The multimodal brain tumor image segmentation benchmark (BRATS). IEEE transactions on medical imaging 2014;34(10):1993-2024.

[2] Zbontar J, Knoll F, Sriram A, Murrell T, Huang Z, Muckley MJ, Defazio A, Stern R, Johnson P, Bruno M. fastMRI: An open dataset and benchmarks for accelerated MRI. arXiv preprint arXiv:181108839 2018.

   
