rmdir build /s /q
mkdir build
cd build
cmake -G "MinGW Makefiles" .. ^
  -DCMAKE_MAKE_PROGRAM="c:/gcc/bin/make.exe" ^
  -DCMAKE_Fortran_COMPILER="c:/gcc/bin/gfortran.exe" ^
  -DCMAKE_C_COMPILER="c:/gcc/bin/gcc.exe" ^
  -DTPL_BLAS_LIBRARIES="../Binaries/lib/libopenblas.dll.a" ^
  -DMYSTRAN_BLAS_LAPACK=SYSTEM ^
  -DTPL_ENABLE_METISLIB=OFF ^
  -DMYSTRAN_USE_EXTERNAL_SUPERLU=ON ^
  -DUSE_SUPERLU_MT=OFF ^
  -DCMAKE_BUILD_TYPE=RELEASE -DMYSTRAN_EXTERNAL_SUPERLU_INCLUDE_DIR="../Binaries/include" 

