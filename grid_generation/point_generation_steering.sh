#!/bin/bash

# ------------ CONFIGURABLE VARIABLES ------------
number_of_process=1  # <-- Change this as needed
# ------------------------------------------------

#------------------------------------------------------------
# dune10k_1x2x6_v6_full

# x_=(10 25 50 75 100 125 150 175 200 225 250 275 300 325 350 360)
# y_=(31.1352 35 40 45 50 93.4056 95 100 105 110 155.676 160 165 170 175 \
#     217.946 220 225 230 235 280.217 285 290 295 300 342.487 345 350 355 360 \
#     404.758 410 415 420 425 467.028 470 475 480 485 529.298 535 540 545 550 \
#     585 590 591.569 595)
# z_=(732.289 843.889 964.679 1076.28 1197.07 1308.67 1357.47)

# testing
x_=( 50. 150. 300. )
y_=( 31.1352 )
z_=( 732.289 )

# Compute total number of points
total_points=$(( ${#x_[@]} * ${#y_[@]} * ${#z_[@]} ))
points_per_file=$(( (total_points + number_of_process - 1) / number_of_process ))  # ceil division

# Initialize counters
count=0
file_index=0
file="myLightSourceSteering_${file_index}.txt"
# Create and initialize the first file
> "$file"
echo "x      y      z      t      dx      dy      dz      dt      p      dp      n" >> "$file"

# Generate files
for x in "${x_[@]}"; do
  for y in "${y_[@]}"; do
    for z in "${z_[@]}"; do
      if (( count % points_per_file == 0 && count != 0 )); then
        ((file_index++))
        file="myLightSourceSteering_${file_index}.txt"
        > "$file"
        echo "x      y      z      t      dx      dy      dz      dt      p      dp      n" >> "$file"
      fi

      echo "$x  $y  $z  0.0  0.0  0.0  0.0  0.0  9.69  0.25  100000" >> "$file"

      ((count++))
    done
  done
done

echo "Generated $((file_index + 1)) files with approximately $points_per_file points each."
