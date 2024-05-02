SEARCH_STRING="Optical Channels positions:  "
INPUT_FILE=../generate/process0/log
INPUT_FILE=../../runWithWiresv6/run1/process0/log
OUTPUT_FILE=protodune_optical_mapping.txt

MAP_HEADER=`grep "${SEARCH_STRING}" ${INPUT_FILE}`
export CHANNEL_NUMBER=${MAP_HEADER#${SEARCH_STRING}}
echo CHANNEL_NUMBER=${CHANNEL_NUMBER}

grep "${SEARCH_STRING}" -A ${CHANNEL_NUMBER} ${INPUT_FILE} | grep -v "${SEARCH_STRING}" > ${OUTPUT_FILE}
echo Created map file ${OUTPUT_FILE}
