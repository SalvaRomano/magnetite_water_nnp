#!/bin/bash
rm bad_folders.txt 2>/dev/null

echo "get totally bad mag folders"
(for i in $(find . -mindepth 3 -type d | sort )  ; do cd  $i; if [ -f "OUTCAR" ]; then echo $i;  tail -n 1 OSZICAR; fi; cd ../../../; done) | grep -v -e  'mag=     0' -e  'mag=    -0' -e 'mag=     1' -e  'mag=    -1' -e  'mag=     2' -e  'mag=     -2' -e 'mag=    56'  -e 'mag=    55' | grep mag -B1 | grep For | tee -a  bad_folders.txt


echo "get DAV folders"
(for i in $(find . -mindepth 3 -type d | sort )  ; do cd  $i; if [ -f "OUTCAR" ]; then echo $i;  tail -n 1 OSZICAR; fi; cd ../../../; done) | grep DAV -B1 | grep _ | tee -a bad_folders.txt



echo "get CGA 250 folders"
# GET CGA 250 folders
(for i in $(find . -mindepth 3 -type d | sort )  ; do cd  $i; if [ -f "OUTCAR" ]; then echo $i;  tail -n 2 OSZICAR; fi; cd ../../../; done) | grep 'CGA: 250' -B1 | grep _ | tee -a bad_folders.txt


echo "remove tot_bad_mag/DAV/CGA folders"
for i in $(cat bad_folders.txt); do echo "removing $i"; rm -r $i;   done

rm bad_folders.txt
