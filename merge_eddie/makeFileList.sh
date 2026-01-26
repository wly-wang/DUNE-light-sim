for FILE in `ls --color=never ../generate_eddie/job*/*.root`
do
  echo $FILE >> List_of_files.txt
done

