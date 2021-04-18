We use markdown to generate documentation for this toolbox. The best way to view/edit markdown files in your computer is to have a markdown editor.

In this file, we provide some notes for make it easier extending the toolbox.



# Notes

1) Create one script to benchmark a method over a specified dataset.

Every method has its own parameters and may need some extra considerations. To cope with these difficulties, a simple strategy of creating a separate script for each method and dataset can be used. For example, the script`benchmark_[Method Name]_[Dataset Name].m` takes care of how to prepare input images from the specified dataset, runs the method, and store the outputs. To run each method, we provide a function in `./Methods` subfolder. This function encapsulates all considarations which should be take before calling the method. For example, The function `./Methods/run_bm4d.m` is called in `benchmark_bm4d_dt1`. 



2) You may want to report the results of your method in a different way. 

For example, you may prefer to use a different set of metrics or you may want to define regions of interests (ROIs) which differ from the ones existed in the dataset. 

To cope with these varieties, report your results in a separate table.



For the dataset 1, we saved some ROIs in `.MAT` files in  `./Metrics/rois_for_dt1` . The first ROI in each one is background ROI and the others are foreground ROIs.



3) To insert your data into the readme file, you can use the following steps:

- Insert  your results into a table in an Excel sheet
  Create an image from the table (E.g., use print screen key and crop the image).
  Save  the image into `./Readme_files` folder. (E.g.,  `./Readme_files/results_table1.png`)
  Insert  it into the markdown using the following command in the source code of the markdown file.

`![name or title](./Readme_files/results_table1.png)`



4) You may find the following software useful:

- To edit or view markdown (`.md` files): Typora® 
- To take a snap shot from the screen: lightshot® 
- To rename multiple files at once: Rename Master ® 



# Demo for How to Define ROIs and Compute Metrics

`Demo_Compute_Metrics.m` shows how to compute the above mentioned metrics:


The other metrics do not need ground-truth image; however, they are computed based on previously selected local regions of interests (ROIs). This script (`Demo_Compute_Metrics`) shows how one can define ROIs and use them to compute these metrics.

The workflow of this script is as follows:

- Read 3 images: input, ground-truth, and output of a denoising method. 
- Compute PSNR and SSIM
- Force user to define  some number of ROIs
- Computes the other metrics (MSR, CNR, ENL, TP, and EP)