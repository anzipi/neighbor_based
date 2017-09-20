su root
echo 1 > /proc/sys/vm/drop_caches

编译：

g++ -fopenmp -o var getVar.cpp -I /soft/share/include/ -L /soft/share/lib -lgdal
g++ -fopenmp -o predP_L predictPoint_line.cpp -I /soft/share/include/ -L /soft/share/lib -lgdal
g++ -fopenmp -o predP_C predictPoint_circle.cpp -I /soft/share/include/ -L /soft/share/lib -lgdal
g++ -fopenmp -o predP_S predictPoint_sfd.cpp -I /soft/share/include/ -L /soft/share/lib -lgdal
g++ -fopenmp -o predP_M predictPoint_mfd.cpp -I /soft/share/include/ -L /soft/share/lib -lgdal
g++ -fopenmp -o predA_C predictArea_circle.cpp -I /soft/share/include/ -L /soft/share/lib -lgdal
g++ -fopenmp -o predA_S predictArea_sfd.cpp -I /soft/share/include/ -L /soft/share/lib -lgdal
g++ -fopenmp -o predA_M predictArea_mfd.cpp -I /soft/share/include/ -L /soft/share/lib -lgdal


