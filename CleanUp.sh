today=`date '+%d_%m_%Y_%H%M'`;
mkdir data/run_${today}
mv InitialState.dat Output.dat Lattice_* Equili* data/run_${today}

