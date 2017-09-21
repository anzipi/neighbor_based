su root
echo 1 > /proc/sys/vm/drop_caches

编译：


g++ -fopenmp -o predP_C predictPoint_circle.cpp -I /soft/share/include/ -L /soft/share/lib -lgdal
g++ -fopenmp -o predP_S predictPoint_sfd.cpp -I /soft/share/include/ -L /soft/share/lib -lgdal
g++ -fopenmp -o predP_M predictPoint_mfd.cpp -I /soft/share/include/ -L /soft/share/lib -lgdal
g++ -fopenmp -o predA_C predictArea_circle.cpp -I /soft/share/include/ -L /soft/share/lib -lgdal
g++ -fopenmp -o predA_S predictArea_sfd.cpp -I /soft/share/include/ -L /soft/share/lib -lgdal
g++ -fopenmp -o predA_M predictArea_mfd.cpp -I /soft/share/include/ -L /soft/share/lib -lgdal


运行：

./predP_M ./heshan/data/demFil.tif ./heshan/data/slope.tif#./heshan/data/plan.tif#./heshan/data/prof.tif#./heshan/data/twi.tif ./heshan/data/samples_yl_38p.csv ./heshan/data/samples_regular.csv X Y org ./heshan/result/pred_mfd_env4_10.txt 10
./predP_S ./heshan/data/demFil.tif ./heshan/data/flowDir.tif ./heshan/data/slope.tif#./heshan/data/plan.tif#./heshan/data/prof.tif#./heshan/data/twi.tif ./heshan/data/samples_yl_38p.csv ./heshan/data/samples_regular.csv X Y org ./heshan/result/pred_env4_20.txt 20
