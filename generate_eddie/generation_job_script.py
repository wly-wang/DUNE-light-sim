import os

# === User Configuration ===
nProcesses = 500
# xVals = [ 0.,10., 20., 30., 40., 50., 100., 150., 200., 250., 280., 300., 310., 320., 330., 340., 350. ]

# yVals = [ 578.909, 570., 560., 550., 518.159, 510., 500., 490., 457.409, 450., 440., 430.,
#           396.659, 390., 380., 370., 335.909, 330., 320., 310., 275.159, 270., 260., 250.,
#           214.41, 210., 200., 190., 153.66, 150., 140., 130., 92.9099, 90., 80., 70., 32.16, 30., 20., 10.]

# zVals = [35.1188, 83.9188, 146.719, 195.519, 267.509, 316.309, 379.109, 427.909, 499.899, 548.699, 611.499, 660.299]

# For testing purposes, use a smaller set of points
xVals=[50.]
yVals=[31.1352]
zVals=[732.289]

totalPositions = len(xVals) * len(yVals) * len(zVals)

print(f"Total positions: {totalPositions}")
positionsPerJob = totalPositions // nProcesses

topDirectory = os.getcwd()
pointList = []
jobIndex = 0

for x in xVals:
    for y in yVals:
        for z in zVals:
            pointList.append([x, y, z])

            if len(pointList) >= positionsPerJob:
                # Create job directory
                jobDirectory = os.path.join(topDirectory, f"job{jobIndex}")
                os.makedirs(jobDirectory, exist_ok=True)

                # Write steering file
                steeringPath = os.path.join(jobDirectory, "myLightSourceSteering.txt")
                with open(steeringPath, "w") as f:
                    f.write("  x      y     z      t      dx      dy      dz      dt      p         dp        n\n")
                    for px, py, pz in pointList:
                        f.write(f"  {px}    {py}     {pz}     0.0     0.0     0.0     0.0     0.0     9.69     0.25     100000\n")

                # Write job script
                jobScriptPath = os.path.join(jobDirectory, "jobscript.sh")
                with open(jobScriptPath, "w") as jobScript:
                    jobScript.write("#!/bin/bash\n")
                    jobScript.write(f"#$ -N lightSimJob{jobIndex}\n")
                    jobScript.write("#$ -cwd\n")
                    jobScript.write("#$ -l h_rt=02:00:00\n")
                    jobScript.write("#$ -l h_vmem=8G\n")

                    jobScript.write("SECONDS=0\n")  # start timer

                    jobScript.write("/cvmfs/oasis.opensciencegrid.org/mis/apptainer/current/bin/apptainer exec "
                                    "--cleanenv "
                                    "-B /cvmfs,/home/s2106059,/opt,/run/user,/etc/hostname,/etc/hosts,/etc/krb5.conf "
                                    "/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:latest "
                                    "/bin/bash -lc "
                                    f"\"source /home/s2106059/setup/dune_light_sim_v10_03_01d02_setup.sh; "
                                    f"cd {jobDirectory}; "
                                    "echo Current directory: $(pwd); "
                                    "echo Starting job...; "
                                    f"lar -c ../protoduneHD_v6_refactored.fcl -n {len(pointList)*100} &> log\"\n")

                    jobScript.write("DURATION=$SECONDS\n")
                    jobScript.write("echo \"Job completed in $(($DURATION / 60)) min $(($DURATION % 60)) sec\" >> log\n")



                # Submit the job
                os.chdir(jobDirectory)
                os.system("qsub jobscript.sh")
                os.chdir(topDirectory)

                # Reset for next chunk
                pointList = []
                jobIndex += 1

# Catch remainder
if len(pointList) > 0:
    jobDirectory = os.path.join(topDirectory, f"job{jobIndex}")
    os.makedirs(jobDirectory, exist_ok=True)

    steeringPath = os.path.join(jobDirectory, "myLightSourceSteering.txt")
    with open(steeringPath, "w") as f:
        f.write("  x      y     z      t      dx      dy      dz      dt      p         dp        n\n")
        for px, py, pz in pointList:
            f.write(f"  {px}    {py}     {pz}     0.0     0.0     0.0     0.0     0.0     9.69     0.25     100000\n")

    jobScriptPath = os.path.join(jobDirectory, "jobscript.sh")
    with open(jobScriptPath, "w") as jobScript:
        jobScript.write("#!/bin/sh\n")
        jobScript.write(f"#$ -N lightSimJob{jobIndex}\n")
        jobScript.write("#$ -cwd\n")
        jobScript.write("#$ -l h_rt=02:00:00\n")
        jobScript.write("#$ -l h_vmem=8G\n")
        jobScript.write("source $HOME/setup/create_sl7_alias\n")
        jobScript.write("source $HOME/setup/activate_apptainer\n")
        jobScript.write("source $HOME/setup/dune_light_sim_v10_03_01d02_setup.sh\n")
        jobScript.write("echo Starting job...\n")
        jobScript.write("lar -c ../protoduneHD_v6_refactored.fcl -n " +
                        f"{len(pointList)*100} &> log\n")

    os.chdir(jobDirectory)
    os.system("qsub jobscript.sh")
    os.chdir(topDirectory)

print(f"All {jobIndex + 1} jobs submitted.")

