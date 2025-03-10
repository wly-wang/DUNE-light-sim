for FILE in `ls --color=never ../generate/process*/*.root`
do
  echo $FILE >> List_of_files.txt
done

for FILE in `ls --color=never ../generate_extra_photon/process*/*.root`
do
  echo $FILE >> List_of_files.txt
done

