for fol in */; do cd $fol;
echo $fol
 for i in run*; do cd $i; echo $i;  tail -n 50 log* | grep NNP -B1 | tail -n 2| grep -v NNP | awk '{print $1}' | tee -a  ../tmp; cd ../; done
echo "Total =  $(awk '{ sum += $1 } END { print sum }'  tmp)"
rm tmp
cd ../; done                                                                                                                                      
