#!/usr/bin/env bash

export OMP_PROC_BIND=spread
export OMP_PLACES=threads

#Event number
ievt=$1

date > ${OUTPUT}/event${ievt}.log


#copy event to temp dir
rsync -a ${WORKDIR}/configs ${TMP}/event${ievt}
rsync -a ${WORKDIR}/models/CCAKE/EoS/Houston ${TMP}/event${ievt}

sed -i "s|entropy-dict-dir = DUMMY|entropy-dict-dir = ${WORKDIR}/misc/PbPb_5020_sigma0d30_3M.dat|g" "${WORKDIR}/configs/trento_config.dat"

cd ${TMP}/event${ievt}

${WORKDIR}/models/trento/build/src/trento -c ${WORKDIR}/configs/trento_config.dat

GRID_MAX=$(grep "grid-max" "${WORKDIR}/configs/trento_config.dat" | awk -F'=' '{print $2}' | awk '{print $1}')
GRID_STEP=$(grep "grid-step" "${WORKDIR}/configs/trento_config.dat" | awk -F'=' '{print $2}' | awk '{print $1}')

# Print the values
echo "Grid Max: $GRID_MAX"
echo "Grid Step: $GRID_STEP"

python3 ${WORKDIR}/models/CCAKE/utilities/trento2ccake.py ${TMP}/event${ievt}/trento_ic/ic0.dat ${GRID_STEP} ${GRID_MAX} ${TMP}/event${ievt}/trento_ccake.dat 


sed -i "s|file: DUMMY|file: ${TMP}/event${ievt}/trento_ccake.dat|g" "${TMP}/event${ievt}/configs/input_parameters_ccake.yaml"
sed -i "s|path: DUMMY|path: ${TMP}/event${ievt}/Houston/|g" "${TMP}/event${ievt}/configs/input_parameters_ccake.yaml"


${WORKDIR}/models/CCAKE/build/ccake ${TMP}/event${ievt}/configs/input_parameters_ccake.yaml ${TMP}/event${ievt}

cd ${WORKDIR}

mkdir -p ${TMP}/event${ievt}/input
mkdir -p ${TMP}/event${ievt}/results

mv ${TMP}/event${ievt}/freeze_out.dat ${TMP}/event${ievt}/input/surface.dat

cp ${WORKDIR}/models/iS3D/iS3D ${TMP}/event${ievt}
cp ${TMP}/event${ievt}/configs/iS3D_parameters.dat ${TMP}/event${ievt}
cp -rf ${WORKDIR}/models/iS3D/PDG ${TMP}/event${ievt}
cp -rf ${WORKDIR}/models/iS3D/deltaf_coefficients ${TMP}/event${ievt}
cp -rf ${WORKDIR}/models/iS3D/generate_delta_f_coefficients ${TMP}/event${ievt}
cp -rf ${WORKDIR}/models/iS3D/tables ${TMP}/event${ievt}

cd ${TMP}/event${ievt}
./iS3D

rm -rf ${TMP}/event${ievt}/PDG
rm -rf ${TMP}/event${ievt}/deltaf_coefficients
rm -rf ${TMP}/event${ievt}/generate_delta_f_coefficients
rm -rf ${TMP}/event${ievt}/tables
rm -rf ${TMP}/event${ievt}/iS3D
rm -rf ${TMP}/event${ievt}/input

mv average_thermodynamic_quantities.dat ${TMP}/event${ievt}/results
mv iS3D_parameters.dat ${TMP}/event${ievt}/results


smash -i ${TMP}/event${ievt}/configs/config.yaml \
          -c "Modi: { List: { File_Directory: ${TMP}/event${ievt}/results } }" \
          -o ${OUTPUT}/event${ievt}/

rm -rf ${OUTPUT}/event${ievt}/tabulations  
#copy all configs 
cp -rf ${WORKDIR}/configs ${OUTPUT}/event${ievt}

echo "Dir" >> ${OUTPUT}/event${ievt}.log
echo pwd >> ${OUTPUT}/event${ievt}.log
echo "Dir contents:" >> ${OUTPUT}/event${ievt}.log

#run trento





date >> ${OUTPUT}/event${ievt}.log
