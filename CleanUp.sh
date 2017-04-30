today=`date '+%d_%m_%Y_%H%M'`;
mkdir data/DipoleDipole/run_${today}
mv InitialState.dat Output.dat Lattice_* Equili* Run* data/DipoleDipole/run_${today}/

