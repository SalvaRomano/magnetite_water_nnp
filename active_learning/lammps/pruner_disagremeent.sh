#!/bin/bash
#launch it in the folder containing disag. in a loop if multiple folders has disagreement.txt
output_file="pruned_disagreement.txt"
interval=1000   #"least time interval between configs"

prev_number=""

if [ -f 'disagreement.txt' ]; then
echo 'disagreement'
for number in $(cat disagreement.txt | awk '{ print $1}' | uniq   ); do
    # Check if the line is empty or contains non-numeric characters
    if [[ -z "$number" || ! "$number" =~ ^[0-9.-]+$ ]]; then
        continue
    fi

    if [ -z "$prev_number" ]; then
        echo "$number"  #first number taken
        prev_number="$number"
    elif (( number - prev_number >= interval )); then
        echo "$number"
        prev_number=$number
    fi
done > "$output_file"
else 
    echo 'no disagreement.txt in the folder'
fi

echo "pruned disagrements: $(cat pruned_disagreement.txt |wc -l)"
