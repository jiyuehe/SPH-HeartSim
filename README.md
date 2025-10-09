# Introduction
------------------------------  
This work is based on the open-source framework SPHinXsys (Smoothed Particle Hydrodynamics for Industrial Complex Systems, https://www.sphinxsys.org). The framework provides numerous demonstration codes for phenomena such as dam breaks, elastic water gates, and heat transfer, as well as a bi-ventricle electromechanical simulation. I have extended the codebase to enable cardiac electromechanical simulations on arbitrary geometries, including patient-specific atria and ventricles.

# How to run simulation
------------------------------  
Open terminal in folder /SPH-HeartSim  
mkdir build && cd build  
cmake .. -DCMAKE_BUILD_TYPE=Release  
(NOTE: if no -DCMAKE_BUILD_TYPE=Release at the end, the compiled code runs much slower because it will be in debug mode)  
make -j$(nproc)  
Then open terminal in folder /SPH-HeartSim/build/sim/bin/  
./sim  

Use software ParaView (https://www.paraview.org/) to view the simulation result  
Results will be in folder /SPH-HeartSim/build/sim/bin/output  
Install ParaView: brew install paraview  
ParaView can open PhysiologyHeart_0000000000.vtp and play the simulation movie  

There are useful codes in the folder "tools".  

# Install dependencies
------------------------------  
## Manjaro Linux
sudo pacman -Syu  
sudo pacman -S base-devel cmake git python python-pip eigen  
sudo pacman -S cmake git gcc-fortran  
sudo pacman -S tbb  
sudo pacman -S boost-libs  

Install Simbody:  
cd Desktop  
git clone https://github.com/simbody/simbody.git  
cd simbody  
mkdir build && cd build  
cmake .. -DCMAKE_BUILD_TYPE=Release  
make -j$(nproc)  
sudo make install  
export CMAKE_PREFIX_PATH=/usr/local  

------------------------------  
## Ubuntu
sudo apt install libsimbody-dev  
sudo apt install libeigen3-dev  
sudo apt install libtbb-dev  
sudo apt install libboost-all-dev  

------------------------------  
## MacOS
brew install eigen  
brew install tbb  
brew install boost  
brew install googletest  

Install somebody:  
git clone https://github.com/simbody/simbody.git  
cd simbody  
mkdir build && cd build  
cmake .. -DCMAKE_BUILD_TYPE=Release  
make -j$(nproc)  
sudo make install  
Now can delete the simbody folder  

Most likely it will take forever to build on a Macbook Air, can try the below, but did not work on my 2025 Macbook Air.
To make the c++ build faster:  
brew install ccache  
export PATH="/opt/homebrew/opt/ccache/libexec:$PATH"  

brew install llvm
then choose the compiler in Visual Studio Code: 
Cmd+Shift+P --> CMake: Select a Kit --> /opt/homebrew/opt/llvm/bin/clang++
