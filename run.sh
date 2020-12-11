#!bin/bash
schemes="upw cnt lxw tow mcc"
for i in $schemes;do
    mkdir $i
    sed -i "s/scheme=.*/scheme=$i/" input.ini
    mpirun -n 3 ./a.out > "$i".out
    mv sgn* $i/
    mv "$i".out $i/
done
