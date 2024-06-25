g++ -c -fPIC wirecell_fiducial_volume.cxx -o wirecell_fiducial_volume.o
g++ -shared -Wl,-soname,lib_wirecell_fiducial_volume.so -o lib_wirecell_fiducial_volume.so  wirecell_fiducial_volume.o
