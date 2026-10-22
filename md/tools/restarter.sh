

a='restart_a'
b='restart_b'
for i in */; do
cd $i
for T in T*; do
cd $T
 for k in run*; do
cd $k
for res in $(ls -lt | awk '{print $9}'); do
    if [ $res == $a ] ||  [ $res == $b ]; then
    mkdir continue
    cp $res ./continue/restart.data
    cd continue
    ln -s ../../../../md.lmp_$T ./md.lmp
    cd ../
    break
fi;
done;
cd ../; done;
cd ../; done;
cd ../; done

