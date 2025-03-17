for FILE in `ls --color=never ../generate_xe/process*/*.root`
do
  echo $FILE >> List_of_files.txt
done