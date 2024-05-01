for FILE in `ls --color=never ../generate/process*/*.root`
do
  echo $FILE >> List_of_files.txt
done
