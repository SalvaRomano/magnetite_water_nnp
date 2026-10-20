echo 'Total'
ls */*/* | grep F | wc -l

echo 'Done'
ls */*/*/OUTCAR | wc -l

echo "to go"
(for i in $(find . -mindepth 3 -type d  | sort)  ; do cd  $i; if [ ! -f "OUTCAR" ]; then  echo $i; fi; cd ../../../; done)  | wc -l

echo 'bad magnetization'
(for i in $(find . -mindepth 3 -type d | sort )  ; do cd  $i; if [ -f "OUTCAR" ]; then echo $i;  tail -n 1 OSZICAR; fi; cd ../../../; done) | grep -v -e  'mag=     0' -e  'mag=    -0' -e 'mag=     1' -e  'mag=    -1' -e  'mag=     2' -e  'mag=     -2' -e 'mag=    56'  -e 'mag=    55' | grep mag | wc -l 

echo 'DAV error'

(for i in $(find . -mindepth 3 -type d | sort )  ; do cd  $i; if [ -f "OUTCAR" ]; then echo $i;  tail -n 1 OSZICAR; fi; cd ../../../; done) | grep DAV | wc -l   #remove |wc -l and add -B1 to see folders 
echo 'CGA 250'
(for i in $(find . -mindepth 3 -type d   | sort ) ; do cd  $i; if [ -f "OSZICAR" ]; then echo $i; tail -n 2 OSZICAR; fi; cd ../../../; done) | grep 'CGA: 250' |wc -l  #remove |wc -l and add -B1 to see folders


