#!/usr/bin/env bash

export OMP_PROC_BIND=spread
export OMP_PLACES=threads

#Event number
ievt=$1

date > ${OUTPUT}/event${ievt}.log


#copy event to temp dir
rsync -a ${WORKDIR}/configs ${TMP}/event${ievt}

cd ${TMP}/event${ievt}

${WORKDIR}/models/trento/build/src/trento -c ${WORKDIR}/configs/trento_config.dat

GRID_MAX=$(grep "grid-max" "${WORKDIR}/configs/trento_config.dat" | awk -F'=' '{print $2}' | awk '{print $1}')
GRID_STEP=$(grep "grid-step" "${WORKDIR}/configs/trento_config.dat" | awk -F'=' '{print $2}' | awk '{print $1}')

# Print the values
echo "Grid Max: $GRID_MAX"
echo "Grid Step: $GRID_STEP"

python3 ${WORKDIR}/models/CCAKE/utilities/trento2ccake.py ${TMP}/event${ievt}/trento_ic/ic0.dat ${GRID_MAX} ${GRID_STEP} ${TMP}/event${ievt}/trento_ccake.dat 


sed -i "s|file: DUMMY|file: ${TMP}/event${ievt}/trento_ccake.dat|g" "${TMP}/event${ievt}/configs/input_parameters_ccake.yaml"

${WORKDIR}/models/CCAKE/build/ccake ${TMP}/event${ievt}/configs/input_parameters_ccake.yaml ${TMP}/event${ievt}


echo "Dir" >> ${OUTPUT}/event${ievt}.log
echo pwd >> ${OUTPUT}/event${ievt}.log
echo "Dir contents:" >> ${OUTPUT}/event${ievt}.log

#run trento





date >> ${OUTPUT}/event${ievt}.log
