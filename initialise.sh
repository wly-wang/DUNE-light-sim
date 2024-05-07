# Choose dune sw version
ls /cvmfs/dune.opensciencegrid.org/products/dune/dunesw/
export MY_DUNE_VERSION=v09_78_01

# Env setup
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup dunesw ${MY_DUNE_VERSION}d00 -q e20:prof
