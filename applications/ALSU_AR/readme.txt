This set of scripts is used to perform the correction chain as explained in https://doi.org/10.1103/PhysRevAccelBeams.22.100702 
For a describtion of the correction chain please have a look at the comments in /Studies/main.m and the paper.

As an example we include the workflow which worked best for us:

- Symbolic links in the user's home directory named 'sc', 'at' and 'MML' point to the folders where the SC toolkit, AT and the MML are located, respectively. This allows for scripts working conveniently on different machines (e.g. laptop, cluster, different users) without modification.

- The files for different machines are located at e.g. ~/sc/applications/ALSU_AR/ 
In this folder all general functions regarding the ALS-U accumulator ring are stored, e.g. for loading the multipole tables ('set*Polynoms_ALSU_AR') or how the LOCO sequence ('performLOCO_ALSU_AR'), apertures ('defineApertures_ALSU_AR') are defined. 

- All multipole tables (in AT notation) are located in [...]/Multipoles/

- All lattices files are located in [...]/Lattices/

- All specific studies are carried out in [...]/Studies/
In this example there are no subfolders (e.g. different studies), thus the main script 'main.m' as well as the corresponding slurm submit script 'submit.sh' can be found here. The paths in 'main.m' have to be modified as well as the cluster properties in 'submit.sh'. A subsequent call of 'sbatch submit.sh main' in the corresponding directory on the cluster should evaluate 100 instances of 'main.m'. The results collected at different stages of the correction chain are saved in e.g. [...]/Studies/$RUNDIR/00001/results.mat together with the matlab printout in 'out' and the potential error message 'err'.

- We found it very beneficial to have the local code development and the cluster always in sync, thus in our case using the command
'rsync -Pavz --delete --max-size=10mb --exclude="/*/*/*/*/*/" ~/sc/ thellert@lrc-xfer.lbl.gov:scratch/sc'
synchronizes all scripts on the cluster with the local machine except all files and folders larger 10MB or with a depth of 5, thus e.g. [...]/Studies/$STUDYNAME/$RUNDIR to avoid unnecesarry file transfers. 

- The function 'crawlClusterJob' collects the data saved by the individual matlab instances on the cluster into one convenient data structure 