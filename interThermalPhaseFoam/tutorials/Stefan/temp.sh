SOLVER=interThermalPhaseChangeFoam
stefan=CFD-PC/interThermalPhaseFoam/tutorials/Stefan
./Allclean
./InitStep &> log.txt &
sleep 15s
echo "killed"
killall $SOLVER
echo "resumed screen"
cd ../
