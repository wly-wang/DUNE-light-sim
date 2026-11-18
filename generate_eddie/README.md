If you've set up your DUNE software environment, you can just run multiprocess.py with python.

Please note that this generation FCL relies on https://github.com/LArSoft/larana/pull/29
LARSOFT_SUITE_v09_78_01 or newer should do it

This takes about 30 hours on 30 CPUs. You can adjust the number of processes in the script to match your machine.
Note that this is extremely lazy multiprocessing for a single compute server, but it does avoid dealing with grid issues.

Take a look at Jiaoyang's repository if you want to do grid submission:
https://github.com/Li-Jiaoyang97/DUNE-light-sim

## Dense radial supplement on Eddie

`radial_campaign.py` prepares a scratch-based, full-optical campaign without
submitting anything by default. Its default 25 cm spacing proposes 12 symmetric
source-Z planes in each radial end region, uses 15 fitter-used X values and all
37 Y values, and automatically omits any Z plane already present in the
baseline sample.

Run `prepare` on the Eddie login host, outside Apptainer:

```bash
python3 radial_campaign.py prepare \
  --campaign-dir /exports/eddie/scratch/$USER/dune-light/radial-dense-25cm-v1
```

Inspect the generated manifest and jobs. Preview submission with:

```bash
python3 radial_campaign.py submit \
  --campaign-dir /exports/eddie/scratch/$USER/dune-light/radial-dense-25cm-v1
```

Only the same command with `--confirm` calls `qsub`. Check progress and create
a merge-ready list of successful raw ROOT files with:

```bash
python3 radial_campaign.py status \
  --campaign-dir /exports/eddie/scratch/$USER/dune-light/radial-dense-25cm-v1
```
