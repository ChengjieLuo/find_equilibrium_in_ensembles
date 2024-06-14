#### in the cluster compile using the fftw installed in vcpkg. Note that -D_GLIBCXX_USE_CXX11_ABI=0 is needed to be consistent.


# LD_LIBRARY_PATH=/usr/ds/gcc-12.1/lib64:${LD_LIBRARY_PATH}

g++ -std=c++11 solve_gibbs_with_input.cpp -I/home/cluo/softwares/vcpkg/installed/x64-linux/include -L/home/cluo/softwares/vcpkg/installed/x64-linux/lib -lfftw3 -O3 -o solve_gibbs_with_input.out

#### locally compile using the fftw installed in homebrew

# g++ -std=c++11 solve_gibbs_with_input.cpp -I/opt/homebrew/include/ -L/opt/homebrew/lib -lfftw3 -O3 -o solve_gibbs_with_input.out