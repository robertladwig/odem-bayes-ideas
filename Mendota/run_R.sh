#!/bin/bash

# make .R directory
mkdir .R

# move Makevars to .R
mv Makevars .R/

# Command to enable modules, and then load an appropriate MP/MPI module
. /etc/profile.d/modules.sh
module load GCC/8.3.0

# untar your R installation. Make sure you are using the right version!
tar -xzf R361.tar.gz
# (optional) if you have a set of packages (created in Part 1), untar them also
tar -xzf rstan.tar.gz

# make sure the script will use your R installation, 
# and the working directory as its home location
export R_MAKEVARS_USER=$PWD/.R/Makevars
export PATH=$PWD/R/bin:$PATH
export RHOME=$PWD/R
export R_LIBS=$PWD/rstan

# run your script
Rscript run_odem_chtc_latest.R
#mpirun -np 8 ./my_script.R
