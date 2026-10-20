#!/bin/bash

# Function to generate a random number between 0 and 1000000 with spacing of 100
generate_random_number() {
  local min=0
  local max=1000000
  local spacing=200

  # Calculate the number of possible values within the range and spacing
  local num_possible_values=$(( (max - min) / spacing + 1 ))

  # Generate a random number within the possible values range
  local random_index=$(( RANDOM % num_possible_values ))

  # Calculate the actual random number based on the index and spacing
  local random_number=$(( min + random_index * spacing ))

  echo $random_number
}

# Number of random numbers to generate (N)
N=80

# File name to write the random numbers
file_name="disagreement.txt"

# Generate N random numbers and write them to the file
for (( i=0; i<N; i++ )); do
  random_number=$(generate_random_number)
  echo $random_number >> "$file_name"
done

echo "Random numbers written to $file_name"

