# designer_v2_dev
Temporary repo to test updates to designer. 

Jenny, Santi, please feel free to mess with this code however you want (in fact I encourage you to do so). Make sure everything works as expected. Ask me questions if you have them. 

Jenny, you might note that I now pull some data directly from the bids .json sidecar (pe direction, pf level). I also included the fast mp here, mostly just so I didnt have to wait for the longer version to run. If you want to switch to the nonlocal mp, just replace `import mpcomplex as mp` with `import mpdenoise_sj as mp` in designer. Again, please feel free to play with this and make whatever improvements you see fit.

## Dev Installation:
After development, designer-v2 will be uploaded to pip. Until then, here are instructions for installing locally.

- Note, designer still depends on FSL for topup and eddy, and mrtix3 for image-wise operations. 

- The current master version of mrtrix3 has a bug that stops external python modules from configuring properly. Therefore for now I am requiring users to install the dev version:\
` cd /my/path`\
` git clone https://github.com/MRtrix3/mrtrix3.git`\
` cd mrtrix3`\
` git checkout dev`\
` ./configure`\
` ./build`\
` export PYTHONPATH=/my/path/mrtrix3/lib`

- Clone this repository into a local directory:\
` cd /my/path`\
` git clone https://github.com/badesar1/designer_v2_dev.git`

- Install python dependancies:\
` cd designer_v2_dev`\
` python setup.py install`


## Example usage for meso data
An example script can be found in the examples folder. It is copied here as well. As you can see, preprocessing and fitting are now split into two separate functions: designer for preprocessing and tmi for fitting. 

```
#!/bin/bash

datapath=/mnt/labspace/Projects/MESO_v2.0/ALLSUBJS_2.0/M9734
meso1=M9734_074YM_DIFF_meso.nii
meso2=M9734_074YM__DIFF_meso_research.nii
pa=M9734_074YM_DIFF_meso_PA.nii

cd $datapath

python designer \
-denoise -algorithm jespersen -extent 7,7,7 \
-degibbs \
-mask \
-scratch designer2_processing_test -nocleanup \
$meso1,$meso2 designer2_test.mif

python tmi \
-DTIparams -DKIparams -WDKI -SMIparams \
-mask designer2_processing_test/brain_mask.nii \
-sigma designer2_processing_test/sigma.nii \
-nocleanup \
designer2_test.mif designer2_params_test
```

## TODO:
### Things that need to be done before this can be used in production:

- **Integration testing**. I created a directory called tests. Typically this has 2 types of scripts, integration test scripts and unit test scripts. We dont really need unit testing here, but we do need integration testing. 
We need a script that will take a few example cases (preferably different resolutions, b-shells, etc) and we run designer in as many permutations as possible. We select different input datatypes, output types, and use every available option in every combination. We test what breaks designer so that we can go in and do bug fixing in an organized manner.

- **TMI**. TMI is for the most part entirely untested. It does produce sensible parameter maps, but this is the most agressive change we are making to designer and we should meet and discuss the proper way to test its functionality. 

- **Documentation**. We need thorough documentation explaining all availabe options. What intput types can you use, what outputs, what different preprocessing steps are available, algorithms, etc.

- **Referencing**. The code for TMI is not yet properly referenced.

- **Comparison to Designer V1**. Now this this is 100% in python, we can expect some subtle changes in performance compared to matlab. We should run a few examples using both this version and the old, with the same processing under the hood and evaluate differences.