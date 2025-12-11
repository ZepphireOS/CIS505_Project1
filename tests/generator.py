#!/usr/bin/env python3
import csv
import random
import sys

# Usage:
#   python points.py 100 out.csv
#   python points.py 100 out.csv square

def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <num_points> <output.csv> [square]")
        sys.exit(1)

    num_points = int(sys.argv[1])
    filename = sys.argv[2]
    use_square = len(sys.argv) >= 4 and sys.argv[3].lower() == "square"

    if num_points <= 0:
        print("num_points must be > 0")
        sys.exit(1)

    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)

        if not use_square:
            # Uniform random points in [0,100] x [0,100]
            for _ in range(num_points):
                x = random.uniform(0.0, 100.0)  # random float in range. [web:33][web:36]
                y = random.uniform(0.0, 100.0)
                writer.writerow([x, y])
        else:
            if num_points < 4:
                print("Need at least 4 points when using 'square'")
                sys.exit(1)

            # Square corners from (0,0) to (100,100)
            x_min, y_min = 0.0, 0.0
            x_max, y_max = 100.0, 100.0

            # First 4 points: corners
            writer.writerow([x_min, y_min])  # bottom-left
            writer.writerow([x_max, y_min])  # bottom-right
            writer.writerow([x_max, y_max])  # top-right
            writer.writerow([x_min, y_max])  # top-left

            # Remaining points: random inside square
            for _ in range(num_points - 4):
                x = random.uniform(x_min, x_max)
                y = random.uniform(y_min, y_max)
                writer.writerow([x, y])

    print(f"Generated {num_points} points into {filename}")

if __name__ == "__main__":
    main()
