# Choose dune sw version
ls /cvmfs/dune.opensciencegrid.org/products/dune/dunesw/
export MY_DUNE_VERSION=v09_92_00

# Env setup
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup dunesw ${MY_DUNE_VERSION}d00 -q e26:prof
