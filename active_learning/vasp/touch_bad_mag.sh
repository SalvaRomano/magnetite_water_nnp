
#echo "rm totally bad mag outcar"
#(for i in $(find . -mindepth 3 -type d | sort )  ; do cd  $i; if [ -f "OUTCAR" ]; then echo $i;  tail -n 1 OSZICAR; fi; cd ../../../; done) | grep -v -e  'mag=     0' -e  'mag=    -0' -e 'mag=     1' -e  'mag=    -1' -e  'mag=     2' -e  'mag=     -2' -e 'mag=    56'  -e 'mag=    55' | grep mag -B1 | grep _ > bad_mag_total.txt
#
#for i in $(cat bad_mag_total.txt); do echo $i; rm $i/O*; done




echo "create  'bad_mag' where iron oxide has slightly wrong mag 56"

(for i in $(find . -mindepth 3 -type d | sort )  ; do cd  $i; if [ -f "OUTCAR" ]; then echo $i;  tail -n 1 OSZICAR; fi; cd ../../../; done) | grep -v -e '0\.00' -e '56\.0' -e '55\.9' | grep -e 'mag=    56' -e 'mag=    55' -B1 | grep _ > bad_mag_iron_folders.txt

for i in $(cat bad_mag_iron_folders.txt); do echo $i ; touch $i/bad_mag;  done

echo "create  'bad_mag' where  water has slightly wrong mag 0"

(for i in $(find . -mindepth 3 -type d | sort )  ; do cd  $i; if [ -f "OUTCAR" ]; then echo $i;  tail -n 1 OSZICAR; fi; cd ../../../; done) | grep -v -e '0\.00' -e '56\.0' -e '55\.9'  | grep -e 'mag=     0' -e 'mag=     1' -e  'mag=     2' -e 'mag=    -' -B1 | grep _ > bad_mag_water_folders.txt

for i in $(cat bad_mag_water_folders.txt); do echo $i ; touch $i/bad_mag; done
