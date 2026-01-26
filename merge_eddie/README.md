You can just source run.sh, or examine the individual commands inside.

Assumes relative position of the generate folder with results inside. Also note there is a hard-coded parameter for the number of generated photons, assuming you use the default. Change if necessary.


source makeFileList.sh
source makeChannelMap.sh

make

bin/FastOpticalMerge List_of_files.txt
