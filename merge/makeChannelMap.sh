SEARCH_STRING="Optical Channels positions:  "
INPUT_FILE=../generate/process0/log
OUTPUT_FILE=protodune_optical_mapping.txt

HEADER=`grep "${SEARCH_STRING}" ${INPUT_FILE}`
CHANNEL_NUMBER=${HEADER#${SEARCH_STRING}}
echo CHANNEL_NUMBER=${CHANNEL_NUMBER}: update src/main.cpp
grep "${SEARCH_STRING}" -A ${CHANNEL_NUMBER} ${INPUT_FILE} | grep -v "${SEARCH_STRING}" > ${OUTPUT_FILE}
