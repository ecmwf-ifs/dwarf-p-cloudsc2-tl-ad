export LC_ALL=C 
module purge
module load pi ecbuild ninja hdf5 python3 gcc 
python3 -m venv ~/cstest
source ~/cstest/bin/activate
pip install --upgrade pip
pip install f90wrap h5py
cd build
rm -rf *
cmake -G Ninja .. && ninja
ninja
cd bin
cd pythonexec
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:.
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.
ipython ./dwarfdriver.py

