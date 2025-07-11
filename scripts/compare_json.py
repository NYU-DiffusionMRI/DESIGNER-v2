import json
import numpy as np
import argparse
from typing import Any, List

def load_and_compare_json(file1_path: str, file2_path: str, tolerance: float = 1e-10) -> bool:
    # Load JSON files
    with open(file1_path, 'r') as f1, open(file2_path, 'r') as f2:
        data1 = json.load(f1)
        data2 = json.load(f2)
    
    differences = []
    
    def compare_values(v1: Any, v2: Any, path: List[str]) -> bool:
        current_path = '.'.join(path)
        
        # Check for different types
        if type(v1) != type(v2):
            differences.append(f"Type mismatch at {current_path}: {type(v1).__name__} vs {type(v2).__name__}")
            return False
            
        # Compare based on type
        if isinstance(v1, (int, str)):
            if v1 != v2:
                differences.append(f"Value mismatch at {current_path}: {v1} vs {v2}")
                return False
        elif isinstance(v1, float) or (hasattr(v1, 'dtype') and np.issubdtype(v1.dtype, np.floating)):
            if abs(v1 - v2) >= tolerance:
                differences.append(f"Value mismatch at {current_path}: {v1} vs {v2} (diff: {abs(v1 - v2)})")
                return False
        elif isinstance(v1, list):
            if len(v1) != len(v2):
                differences.append(f"List length mismatch at {current_path}: {len(v1)} vs {len(v2)}")
                return False
            return all(compare_values(x, y, path + [f"[{i}]"]) for i, (x, y) in enumerate(zip(v1, v2)))
        elif isinstance(v1, dict):
            # Check for key differences
            keys1, keys2 = set(v1.keys()), set(v2.keys())
            if keys1 != keys2:
                missing_in_2 = keys1 - keys2
                missing_in_1 = keys2 - keys1
                if missing_in_2:
                    differences.append(f"Keys missing in second JSON at {current_path}: {missing_in_2}")
                if missing_in_1:
                    differences.append(f"Keys missing in first JSON at {current_path}: {missing_in_1}")
                return False
            return all(compare_values(v1[k], v2[k], path + [k]) for k in v1)
        return True

    # Compare the JSONs
    are_equal = compare_values(data1, data2, [])
    
    if are_equal:
        print("JSONs are equal")
    else:
        print("JSONs are different. Differences found:")
        for diff in differences:
            print(f"- {diff}")
    
    return are_equal

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare two JSON files for equality')
    parser.add_argument('file1', type=str, help='Path to the first JSON file')
    parser.add_argument('file2', type=str, help='Path to the second JSON file')
    parser.add_argument('--tolerance', type=float, default=1e-10, 
                      help='Tolerance for floating point comparisons (default: 1e-10)')
    
    args = parser.parse_args()
    load_and_compare_json(args.file1, args.file2, args.tolerance)
