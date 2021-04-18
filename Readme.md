# Private Notes

- Use `Typora` software to edit this file
- List of available methods and their reference and citations should be added
- Please keep it as simple as possible and extend it as you want.
- 

# 3D Retinal OCT Denoising Package

This is a flexible and easy to use package for 3-D denoising of retinal optical coherence tomography (OCT) images. It is also equipped with some facilities to quantitatively assess the quality of denoising.



**Requirements:**
The toolbox is almost self-contained. We tested it on a 64bit PC with Windows 10 and MATLAB 2019. For other architectures, you may need to recompile MEX files in some packages. 



## Available Methods and Datasets

**Methods**

- Sparse K-SVD [1] and [2] (For K-SVD and wavelet based initialized version see [HERE](https://sites.google.com/site/rahelekafieh/research/state-of-the-art-method-for-oct-denoising/).)
- V-BM4D [3] and [4]
- WMF [5]
- BM4D [6] and [7] 
- Tensor Dictionary Learning (Tensor DL / TDL) [8] 
- Multiscale BM4D (MS BM4D) [9] implemented via `dwt3` and `dualtree3` 

**Datasets**

- Bioptigen Images [10]
- Topcon Images [11]

### References

[1]	Rubinstein, Ron, Michael Zibulevsky, and Michael Elad. "Double sparsity: Learning sparse dictionaries for sparse signal approximation." *IEEE Transactions on signal processing* 58.3 (2009): 1553-1564.

[2]    Kafieh, Raheleh, Hossein Rabbani, and Ivan Selesnick. "Three dimensional data-driven multi scale atomic representation of optical coherence tomography." *IEEE transactions on medical imaging* 34.5 (2014): 1042-1062. 

[3]	Maggioni, Matteo, et al. "Video denoising using separable 4D nonlocal spatiotemporal transforms." *Image Processing: Algorithms and Systems IX*. Vol. 7870. International Society for Optics and Photonics, 2011.

[4]	Maggioni, Matteo, et al. "Video denoising, deblocking, and enhancement through separable 4-D nonlocal spatiotemporal transforms." *IEEE Transactions on image processing* 21.9 (2012): 3952-3966.

[5]    Mayer, Markus A., et al. "Wavelet denoising of multiframe optical coherence tomography data." *Biomedical optics express* 3.3 (2012): 572-589.

[6]	Maggioni, Matteo, and Alessandro Foi. "Nonlocal transform-domain denoising of volumetric data with groupwise adaptive variance estimation." *Computational Imaging X*. Vol. 8296. International Society for Optics and Photonics, 2012.

[7]	Maggioni, Matteo, et al. "Nonlocal transform-domain filter for volumetric data denoising and reconstruction." *IEEE transactions on image processing* 22.1 (2012): 119-133.

[8]	Peng, Yi, et al. "Decomposable nonlocal tensor dictionary learning for multispectral image denoising." *Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition*. 2014.

[9]     ............................................................................. OUR WORK .....................

[10]   Fang, Leyuan, et al. "Fast acquisition and reconstruction of optical coherence tomography images via sparse representation." *IEEE transactions on medical imaging* 32.11 (2013): 2034-2049.

[11]   ........................ Reference to TOPCON dataset ....................



## Conventions of this package

1 - All datasets are saved into `./Datasets` folder. 

- The first dataset (`dt1`) is stored in `./Datasets/dt1_Bioptigen_SDOCT`.
- The second dataset (`dt2`) is stored in `./Datasets/dt2_topcon_oct1000_seg_normal`.

2 - Results are being saved into `./Results` folder in which a subfolder will be created for each experiment.

3 - All of the methods are stored in `./Methods`. The following functions are provided for running the methods in an easy and unified manner:

- `benchmark_X_on_dt1`: Runs the denoising method `X` on the 1st dataset.  A function handle is used (as a input argument) to specify the denoising method (`X`).
- `benchmark_X_on_dt2`: Runs the denoising method `X` on the 2nd dataset.  

4 - There are some utility functions for objective image quality assessment. These functions are stored in `./Metrics`. After running a denoising method, you can use the following functions to evaluate the metrics for each dataset in an easy and unified manner:

- `evaluate_metrics_dt1`
- `evaluate_metrics_dt2`



## Some Examples

There are some scripts to easily run the denoising methods on a dataset:

- 
- `Run_benchmarks_dt2` runs the selected methods on the second dataset



### Run a denoising method on the 1st dataset

**NOTE**: In `Run_benchmarks_dt1`, we provide a simple script to run every methods on the first dataset. 

#### Example 1:

Run the original BM4D on the 1st dataset to denoise two volumes which are specified by `test_indices`.

```matlab
addpath(genpath('./Methods'));
addpath('./Metrics/utils');

params.test_indices = [1,3] % denoise volume #1 and #3
params.save_mat = false; % Don't save the denoised volume as a separate MAT-file

X = @run_bm4d; % handle to the runner function for the denoising method

output_folder_name = 'benchmark_bm4d_dt1'; % it will be created in `./Results`
benchmark_X_on_dt1(output_folder_name,X,params)
evaluate_metrics_dt1(output_folder_name,params);
```

Outputs:

1. The denoised images are saved into <u>the output folder</u> (`./Results/benchmark_bm4d_dt1`)
2. The metric results are shown in a table (like below), and they are also stored in a separate Excel (`benchmark_bm4d_dt1.xlsx`) file inside the output folder. 

![](./Readme_files/benchmark_bm4d_dt1_sample.png)

3. All of the command line prompts are saved into a text file (`benchmark_bm4d_dt1.txt`) inside the output folder.



### Run a denoising method on the 2nd dataset

**NOTE**: `Run_benchmarks_dt2` contains full commands to run methods on the first dataset.



- Run Sparse K-SVD on the 2nd dataset to denoise frame number 10 and 60 in the 3rd volume .

```matlab
addpath(genpath('./Methods'));
addpath('./Metrics/utils');

params.test_indices = [3];
params.frame_numbers = ones(length(params.test_indices),1)*[10, 60];

params.n_frames = 5; % greater than or equal to 2
params.valid_rows = 150:512 - 1; % images in this dataset are cropped before processing

params.save_mat = false;
    
output_folder_name = 'benchmark_ksvds_dt2';

X = @run_ksvds;
params.get_params = @get_params_ksvds;

params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2 - 1;

output_folder_name = 'benchmark_ksvds_dt2';
benchmark_X_on_dt2(output_folder_name,X,params)
evaluate_metrics_dt2(output_folder_name,params);
```





# Extra notes

## Signatures of functions called in `benchmark_X_on_...`

Runner function for each denoising method:

`[denoised_imgs,run_time] = run_...(noisy_imgs,params)`

Noise estimation method used for a denoising method:

`sigma_value = estimate_noise_...(noisy_imgs)`

Getting specific parameters of a denoising method:

`params = get_params_....(noisy_imgs,sigma_value)`





