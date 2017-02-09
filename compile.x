
#g++ -I/home/niepin/code/6.software/lammps-30Jul16/src/ -I/home/niepin/code/6.software/lammps-30Jul16/src/STUBS/ -I/home/niepin/code/6.software/lammps-30Jul16/examples/COUPLE/library/ -c -Wall wrap.cpp
g++ -I/home/niepin/code/6.software/lammps-30Jul16/src/  -I/home/niepin/code/6.software/lammps-30Jul16/src/STUBS/ -c -Wall wrap.cpp
g++ -L/home/niepin/code/6.software/lammps-30Jul16/src/  -L/home/niepin/code/6.software/lammps-30Jul16/src/STUBS/ -L/home/niepin/code/6.software/lammps-30Jul16/examples/COUPLE/library/ wrap.o -llammps_serial -lmpi_stubs -lpthread -o wrap.x
