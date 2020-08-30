gfortran -c -O3 -Wall -Wno-unused-dummy-argument -Wno-unused-function -Wno-uninitialized -Wno-conversion -Wno-unused-variable -static -ffree-line-length-0 -std=gnu -o dof_kriging.o dof_kriging.f08 -Igalib\include
gfortran -c -O3 -Wall -Wno-unused-dummy-argument -Wno-unused-function -Wno-uninitialized -Wno-conversion -Wno-unused-variable -static -ffree-line-length-0 -std=gnu -o check_krig.o check_krig.f08 -Igalib\include
gfortran -O3 -Wall -Wno-unused-dummy-argument -Wno-unused-function -Wno-uninitialized -Wno-conversion -Wno-unused-variable -static -ffree-line-length-0 -std=gnu -o check_krig.exe dof_kriging.o check_krig.o -Lgalib\lib -lfec -llapack -lrefblas
copy check_krig.exe ..
