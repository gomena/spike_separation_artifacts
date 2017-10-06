
# Large-scale separation of neural spikes from stimulation artifacts code

This repo contains (matlab) sample code for doing spike sorting with stimulation artifacts.
* See `Example.m` for the full execution of the code, along with details.
* `Art0.mat` contains the initial artifact that is used to find the model hyperparameteres. This Artifact is a three dimensional tensor (`Art0=Art0(j,e,t)`) with dimensions E=512 (electrodes), T=55 (time samples), J=34 (number of different electrical stimuli). This initial estimate was built using equation (7) in the paper. Also, hyperparameters based on this estimate are computed according to equation (8).
* `Data.mat` contains the actual data where spikes are to be found. Specifically
1. The fourth dimensional tensor `TracesAll(j,n_j,t,e)`  of recordings over array for all trials (`n_j=51` is the number of trials) at all amplitudes of stimulation.
2. The EI of 24 neurons (figure 1 here and in the paper), represented as a the cell array `templates`. Each of these EIs is represented as a matrix `V(e,t)` which states how a spike is recorded in each electrode.
3. The list of stimulation amplitudes, `listAmps`
4. Stimulating electrode indexes `stimElecs`(here, 404).
5. `breakpoints` contains the indexes of stimuli at which breakpoints (Sudden changes in the artifact measured in the stimulating electrode) occurred

![Figure 1: spatial arrangement of 24 EIs](/Utils/templates.png)

* After running `Example.m` two outputs should appear:
  * A sample of the spatial profile of the artifact, as figure 3 in the paper.
  ![](/Utils/Artifact.png)
  * The infered responses (activation curves) of each of the neurons of the array, and their EIs (another representation of the first figure shown here)
   ![](/Utils/NeuronSpikes.png)  



## Array support
This code is based on a 512-array, the one depicted below(each circle represents and electrode, and the number its index)
![figure 1](/Utils/array.png)

* The structure Array in `/Utils`, `Array.m` contains array-specific information (e.g. positions, spacing, etc). Specifics for a distinct, 519-array (see below) are also shown there, to illustrate that the algorithm depends on the array through those parameters, and therefore, extensions are straightforward as long as assumptions hold.
![](/Utils/array519.png)
