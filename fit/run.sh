RUN=run3
mkdir plots
cp ../mergeWithWiresv6/CompiledMerger/${RUN}/output.${RUN}.root output.root
root Semi_Mode_Gen_v2.C
mv plots ${RUN}defaultPlusRange
