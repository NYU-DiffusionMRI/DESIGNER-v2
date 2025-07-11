#!/bin/bash

# Check if correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <directory1> <directory2>"
    echo "Example: $0 tests/benchmark tests/tmp_benchmark"
    exit 1
fi

dir1="$1"
dir2="$2"

# Check if both directories exist
if [ ! -d "$dir1" ]; then
    echo "Error: Directory '$dir1' does not exist"
    exit 1
fi

if [ ! -d "$dir2" ]; then
    echo "Error: Directory '$dir2' does not exist"
    exit 1
fi

# Get the directory of this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Initialize counters
total_files=0
different_files=0

# Find all JSON files in the first directory and compare with corresponding files in second directory
find "$dir1" -type f -name "*.json" | while read -r file1; do
    # Get the relative path from dir1
    rel_path="${file1#$dir1/}"
    file2="$dir2/$rel_path"
    
    # Check if corresponding file exists in dir2
    if [ -f "$file2" ]; then
        echo "Comparing: $file1 with $file2"
        total_files=$((total_files + 1))
        
        # Run the Python comparison script
        if ! python3 "$SCRIPT_DIR/compare_json.py" "$file1" "$file2"; then
            different_files=$((different_files + 1))
            echo "----------------------------------------"
        fi
    else
        echo "Warning: No matching file found for: $file1"
    fi
done

echo "========================================="
echo "Comparison Summary:"
echo "Total files compared: $total_files"
echo "Files with differences: $different_files" 