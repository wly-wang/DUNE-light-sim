If you've set up your DUNE software environment, you can just run multiprocess.py with python.

Please note that this generation FCL relies on https://github.com/LArSoft/larana/pull/29
LARSOFT_SUITE_v09_78_01 or newer should do it

This takes about 30 hours on 30 CPUs. You can adjust the number of processes in the script to match your machine.
Note that this is extremely lazy multiprocessing for a single compute server, but it does avoid dealing with grid issues.

Take a look at Jiaoyang's repository if you want to do grid submission:
https://github.com/Li-Jiaoyang97/DUNE-light-sim
