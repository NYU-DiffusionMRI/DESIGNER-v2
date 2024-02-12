This folder contains Docker files for the main designer dependencies.

- `Dockerfile_mrtrix`: currently builds de dev version of mrtrix. 
   The Container is pushed to [twom/mrtrix3:dev-latest](https://hub.docker.com/r/twom/mrtrix3/tags)

- `Dockerfile_fsl`: FSL installations. Currently [twom/fsl:6.0](https://hub.docker.com/r/twom/fsl/tags) 

- `Dockerfile_deps`: Combies FSL and Mrtrix3


## Updating the dependencies
These dependencies are used to create dependencie container, so designer containers can be easily and quickly build.


