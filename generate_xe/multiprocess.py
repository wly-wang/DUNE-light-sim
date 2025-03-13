import os, subprocess

def startJob( pointList, processCount ):

  print( "Starting process with " + str( len( pointList ) ) + " positions" )

  # Make subdirectory for process
  workingDirectory = os.getcwd()
  processDirectory = os.path.join( workingDirectory, "process"+str( processCount ) )
  os.makedirs( processDirectory )
  os.chdir( processDirectory )

  # Make steering file
  steeringFile = open( "myLightSourceSteering.txt", "w" )
  steeringFile.write( "  x      y     z      t      dx      dy      dz      dt      p         dp        n\n" )
  for x, y, z in pointList:
    steeringFile.write( "  "+str(x)+"    "+str(y)+"     "+str(z)+"     0.0     0.0     0.0     0.0     0.0     7.13     0.25     100000\n" )

  # Execute the LArSoft command (or a test)
  command = "lar -c ../protoduneHD_v6_refactored.fcl -n " + str( len( pointList )*100 ) + " &> log" # full
  process = subprocess.Popen( command, shell=True )

  os.chdir( workingDirectory )
  return process

# Photon origins 
# a quadrant of coverage
xVals = [ -350., -340., -330., -320., -310., -300., -280., -250., -200., -150., -100., -50.,
          0., 50., 100., 150., 200., 250., 280., 300., 310., 320., 330., 340., 350. ]

yVals = [ 578.909, 570., 560., 550., 518.159, 510., 500., 490., 457.409, 450., 440., 430.,
          396.659, 390., 380., 370., 335.909, 330., 320., 310., 275.159, 270., 260., 250.,
          214.41, 210., 200., 190., 153.66, 150., 140., 130., 92.9099, 90., 80., 70., 32.16, 30., 20., 10.]

zVals = [ 660.299, 548.699, 427.909, 316.309, 195.519, 35.1188]

# Full Coverage 
# xVals = [ -350., -340., -330., -320., -310., -300., -280., -250., -200., -150., -100., -50.,
#           0., 50., 100., 150., 200., 250., 280., 300., 310., 320., 330., 340., 350. ]

# yVals = [ 578.909, 570., 560., 550., 518.159, 510., 500., 490., 457.409, 450., 440., 430.,
#           396.659, 390., 380., 370., 335.909, 330., 320., 310., 275.159, 270., 260., 250.,
#           214.41, 210., 200., 190., 153.66, 150., 140., 130., 92.9099, 90., 80., 70., 32.16, 30., 20., 10., 
#           -10., -20., -30., -32.16, -70., -80., -90., -92.9099, -130., -140., -150., -153.66, -190., -200.,
#           -210., -214.41, -250., -260., -270., -275.159, -310., -320., -330., -335.909, -370., -380., -390., 
#           -396.659, -430., -440., -450., -457.409, -490., -500., -510., -518.159, -550., -560., -570., -578.909]

# zVals = [ 1357.47, 1308.67, 1245.87, 1197.07, 1125.08, 1076.28, 1013.48, 964.679, 892.689, 843.889, 781.089, 
#           732.289, 660.299, 611.499, 548.699, 499.899, 427.909, 379.109, 316.309, 267.509, 195.519, 146.719, 
#           83.9188, 35.1188 ]

# Define number of processes to use
nProcesses = 25
totalPositions = len( xVals ) * len( yVals ) * len( zVals )
print( "Total positions: " + str( totalPositions ) )

# Divide up all the position combinations into a set of processes
pointList = []
processes = []
for x in xVals:
  for y in yVals:
    for z in zVals:
      pointList.append( [ x, y, z ] )

      # Launch a set of positions as a process
      if len( pointList ) >= totalPositions / nProcesses:
        processes.append( startJob( pointList, len( processes ) ) )
        pointList = []

# Catch remainder
if len( pointList ) > 0:
  processes.append( startJob( pointList, len( processes ) ) )

# Wait for completion
print( "Waiting for " + str( len( processes ) ) + " processes" )
for process in processes:
  process.wait()
print( "All processes complete" )
