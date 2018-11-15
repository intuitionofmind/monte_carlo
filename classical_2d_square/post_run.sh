#!/bin/sh

SIZE=3232
START_BETA=0.01
STEP=0.01
NUM=100

mv mag.dat ising_mag${numEval}_size${SIZE}_beta${START_BETA}_${STEP}_${NUM}.dat
mv *.dat ${HOME}/Downloads/data/

# while [ ${folder} != ${end} ]; do
    # cd ${DIR}/${folder}
    # echo ${folder}
    # if [ $((${folder} % 2)) = 0 ]; then
        # mv eigenvalues.dat ladderJ03_eigenvalues${numEval}_size${size}_alpha${startAlpha}_${step}_${num}_bcPO_sigmaFalse_anisoBoth.dat
        # mv eigenvectors.dat ladderJ03_eigenvectors${numEval}_size${size}_alpha${startAlpha}_${step}_${num}_bcPO_sigmaFalse_anisoBoth.dat
    # elif [ $((${folder} % 2)) = 1 ]; then
        # # mv eigenvalues.dat ladder_eigenvalues${numEval}_size${size}_alpha${startAlpha}_${step}_${num}_bcPO_sigmaFalse_anisoJ.dat
        # # mv eigenvectors.dat ladder_eigenvectors${numEval}_size${size}_alpha${startAlpha}_${step}_${num}_bcPO_sigmaFalse_anisoJ.dat
        # mv eigenvalues.dat ladderJ03_eigenvalues${numEval}_size${size}_alpha${startAlpha}_${step}_${num}_bcPO_sigmaTrue_anisoBoth.dat
        # mv eigenvectors.dat ladderJ03_eigenvectors${numEval}_size${size}_alpha${startAlpha}_${step}_${num}_bcPO_sigmaTrue_anisoBoth.dat
        # mv *.dat ${HOME}/data
    # fi
    # cd ..
    # folder=$((${folder} + 1))
# done

exit 0
