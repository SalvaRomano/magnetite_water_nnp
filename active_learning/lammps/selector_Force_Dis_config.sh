
N=10  #one over N survive (opposite add !)

index=1
for i in  For*; do
  if   (( !  $index % $N == 0)); then
     echo $i
        rm -r  $i
fi
  ((index++))
done
