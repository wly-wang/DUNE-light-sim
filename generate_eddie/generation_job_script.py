import os

# === User Configuration ===
nProcesses = 1000
# xVals = [ 0.,10., 20., 30., 40., 50., 100., 150., 200., 250., 280., 300., 310., 320., 330., 340., 350. ]

# yVals = [ 578.909, 570., 560., 550., 518.159, 510., 500., 490., 457.409, 450., 440., 430.,
#           396.659, 390., 380., 370., 335.909, 330., 320., 310., 275.159, 270., 260., 250.,
#           214.41, 210., 200., 190., 153.66, 150., 140., 130., 92.9099, 90., 80., 70., 32.16, 30., 20., 10.]

# zVals = [35.1188, 83.9188, 146.719, 195.519, 267.509, 316.309, 379.109, 427.909, 499.899, 548.699, 611.499, 660.299]

# For testing purposes, use a smaller set of points
# xVals=[5.0]
# yVals=[3.11352]
# zVals=[73.2289]

# xVals = [ 1.25, 3.75, 6.25, 8.75, 11.25, 13.75, 16.25, 18.75, 21.25, 23.75, 26.25, 28.75, 31.25, 33.75, 36.25, 38.75, 41.25, 43.75, 46.25, 
#           48.75, 56.25, 61.25, 63.75, 68.75, 71.25, 78.75, 81.25, 86.25, 93.75, 96.25, 101.25, 106.25, 108.75, 113.75, 116.25, 118.75, 123.75, 
#           131.25, 138.75, 143.75, 146.25, 148.75, 151.25, 156.25, 166.25, 168.75, 178.75, 181.25, 183.75, 191.25, 193.75, 201.25, 206.25, 211.25, 
#           213.75, 218.75, 231.25, 233.75, 236.25, 243.75, 253.75, 258.75, 261.25, 271.25, 276.25, 281.25, 288.75, 303.75, 306.25, 308.75, 316.25, 
#           318.75, 323.75, 326.25, 341.25, 343.75, 348.75]

# yVals = [44.0, 120.0, 197.0, 273.0, 341.0, 409.0, 483.0, 560.0, 636.0, 711.0, 783.0, 856.0]

# zVals = [ 35.1188, 195.519, 379.109, 548.699, 732.289, 892.689, 1076.28, 1245.87, 1429.46, 1589.86, 
#           1773.45, 1943.04, 2126.63, 2287.03, 2470.62, 2640.21, 2823.8, 2984.2, 3167.79, 3337.38, 3520.97, 
#           3681.37, 3864.96, 4034.55, 4218.14, 4378.54, 4562.13, 4731.72, 4915.31, 5075.71, 5259.3, 5428.89, 
#           5612.48, 5772.88]


xVals = [ 0.,10., 20., 30., 40., 50., 100., 150., 200., 250., 280., 300., 310., 320., 330., 340., 350. ]

yVals = [ 31.135, 46.703, 62.27, 77.838, 93.406, 108.973, 124.541, 140.108, 155.676, 171.243, 186.811, 202.379, 
          217.946, 233.514, 249.082, 264.649, 280.217, 295.784, 311.352, 326.92, 342.487, 358.055, 373.623, 389.19, 
          404.758, 420.325, 435.893, 451.461, 467.028, 482.596, 498.163, 513.731, 529.298, 544.866, 560.433, 576.001, 591.569 ]

zVals= [ 35.119, 146.719, 267.509, 379.109, 499.899, 611.499, 732.289, 843.889, 964.679, 1076.28, 1197.07, 1308.67, 1429.46, 
         1541.06, 1661.85, 1773.45, 1894.24, 2005.84, 2126.63, 2238.23, 2359.02, 2470.62, 2591.41, 2703.01, 2823.8, 2935.4, 
         3056.19, 3167.79, 3288.58, 3400.18, 3520.97, 3632.57, 3753.36, 3864.96, 3985.75, 4097.35, 4218.14, 4329.74, 4450.53, 
         4562.13, 4682.92, 4794.52, 4915.31, 5026.91, 5147.7, 5259.3, 5380.09, 5491.69, 5612.48, 5724.08 ]



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
                    jobScript.write("#$ -l h_rt=10:00:00\n")
                    jobScript.write("#$ -l h_vmem=8G\n")

                    jobScript.write("SECONDS=0\n")  # start timer

                    jobScript.write("/cvmfs/oasis.opensciencegrid.org/mis/apptainer/current/bin/apptainer exec "
                                    "--cleanenv "
                                    "-B /cvmfs,/home/s2106059,/exports/eddie/scratch/wwang,/opt,/run/user,/etc/hostname,/etc/hosts,/etc/krb5.conf "
                                    "/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:latest "
                                    "/bin/bash -lc "
                                    f"\"source /home/s2106059/setup/dune_light_sim_v10_16_00d00_setup.sh; "
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
        jobScript.write("#$ -l h_rt=10:00:00\n")
        jobScript.write("#$ -l h_vmem=8G\n")
        jobScript.write("source $HOME/setup/create_sl7_alias\n")
        jobScript.write("source $HOME/setup/activate_apptainer\n")
        jobScript.write("source $HOME/setup/dune_light_sim_v10_16_00d00_setup.sh\n")
        jobScript.write("echo Starting job...\n")
        jobScript.write("lar -c ../protoduneHD_v6_refactored.fcl -n " +
                        f"{len(pointList)*100} &> log\n")

    os.chdir(jobDirectory)
    os.system("qsub jobscript.sh")
    os.chdir(topDirectory)

print(f"All {jobIndex + 1} jobs submitted.")

