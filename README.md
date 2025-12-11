# Convex Hull Algorithms

## Overview
This project implements two convex hull algorithms in C++:
1. **Jarvis March (Gift Wrapping)** - O(n·h) time complexity
2. **Graham Scan** - O(n log n) time complexity

Both algorithms compute the smallest convex polygon that encloses a set of 2D points and include comprehensive performance metrics.

## System Requirements

- **Operating System**: Windows, macOS, or Linux
- **Compiler**: C++ compiler with C++11 support
  - Windows: MinGW-g++, Visual Studio, or MSYS2
  - macOS: Clang
  - Linux: g++
- **RAM**: 512 MB minimum
- **Disk Space**: 10 MB minimum

## Compilation Instructions

After cloning the repository, compile the code files by running the following.

```bash
g++ -std=c++11 -o jarvis.exe jarvis.cpp
g++ -std=c++11 -o graham.exe graham.cpp
```

g++ can be replaced by any alternatives based on the hardware used (eg. clang++ for macOS).

## Creating Input Files

### Input File Format
Create a CSV (Comma-Separated Values) file with 2D coordinates. Each line should contain one point with x and y coordinates separated by a comma.

### Input File Rules:
1. **Format**: Each line must have exactly 2 values separated by a comma
2. **Values**: Must be numeric (integers or floating-point numbers)
3. **Whitespace**: Spaces around values are allowed
4. **Empty Lines**: Empty lines are ignored
5. **Scientific Notation**: Supported (e.g., `1.5e10`)

### Input File Generation
For random point generation within the range of 0.0 to 100.0, run the python file generator.py within the 'tests' folder. This was run on Python 3.12.3 and should not need extra packages. The term 'square' can be optionally added to the end to enclose all the points within a large square.

```bash
python3 tests/generator.py \<num_points\> \<output.csv\> \[square\]
```

## Running the Programs

### Prepare Your Input File
Create a CSV file with your 2D points.

### Run the Program

Run the program by running the executable. Rename it as required. Keep in mind that adding the input csv file and output file name to the program call is also a viable means of running the program. If not added, they will be asked by the program directly.

```bash
./jarvis.exe
./graham.exe
```

or

```bash
./jarvis.exe \<inputFilename\> \<outputFilename\>
./graham.exe \<inputFilename\> \<outputFilename\>
```

## Accessing Tests

Tests for the files have been placed in the 'tests' folder. The output of each test can be found in the algorithm's respective '*_out' folder. 