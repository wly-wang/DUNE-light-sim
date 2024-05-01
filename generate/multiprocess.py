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
    steeringFile.write( "  "+str(x)+"    "+str(y)+"     "+str(z)+"     0.0     0.0     0.0     0.0     0.0     9.69     0.25     100000\n" )

  # Execute the LArSoft command (or a test)
  command = "lar -c ../protoduneHD_v6_refactored.fcl -n " + str( len( pointList )*100 ) + " &> log" # full
  process = subprocess.Popen( command, shell=True )

  os.chdir( workingDirectory )
  return process

# Photon origins
xVals = [ -350., -340., -330., -320., -310., -300., -280., -250., -200., -150., -100., -50.,
          0., 50., 100., 150., 200., 250., 280., 300., 310., 320., 330., 340., 350. ]
yVals = [ 578.909, 570., 560., 550., 518.159, 510., 500., 490., 457.409, 450., 440., 430.,
          396.659, 390., 380., 370., 335.909, 330., 320., 310., 275.159, 270., 260., 250.,
          214.41, 210., 200., 190., 153.66, 150., 140., 130., 92.9099, 90., 80., 70., 32.16, 30., 20., 10. ]
zVals = [ 427.071, 377.921, 316.671, 267.521 ]

# Define number of processes to use
nProcesses = 30
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
