import json
import numpy as np
import argparse
from typing import Any, List, Dict
from dataclasses import dataclass, field

@dataclass
class DiffStats:
    max_rel_diff: float = 0.0
    max_abs_diff: float = 0.0
    path: str = ""
    value1: float = 0.0  # Store first value for printing when equal
    value2: float = 0.0  # Store second value for printing when equal

@dataclass
class StatsSummary:
    b0_mean: DiffStats = field(default_factory=DiffStats)
    b0_std: DiffStats = field(default_factory=DiffStats)
    fa_mean: DiffStats = field(default_factory=DiffStats)
    fa_std: DiffStats = field(default_factory=DiffStats)
    wm_ratio: DiffStats = field(default_factory=DiffStats)

def load_and_compare_json(file1_path: str, file2_path: str, tolerance: float = 1e-10) -> bool:
    # Load JSON files while preserving order
    with open(file1_path, 'r') as f1, open(file2_path, 'r') as f2:
        data1 = json.load(f1)
        data2 = json.load(f2)
    
    differences = []
    stats_summary = StatsSummary()
    
    def get_list_index_description(index: int) -> str:
        """Return descriptive name for special list indices."""
        if index == 0:
            return "mean"
        elif index == 1:
            return "std"
        return str(index)
    
    def format_list_path(path: List[str]) -> str:
        """Format path with descriptive names for special indices."""
        result = []
        for p in path:
            if p.startswith('[') and p.endswith(']'):
                # Extract index from "[N]" format
                try:
                    index = int(p[1:-1])
                    desc = get_list_index_description(index)
                    if desc.isdigit():
                        result.append(p)  # Keep original format for non-special indices
                    else:
                        result.append(f"[{index}:{desc}]")
                except ValueError:
                    result.append(p)
            else:
                result.append(p)
        return '.'.join(result)

    def update_stats_summary(path: str, abs_diff: float, rel_diff_v1: float, rel_diff_v2: float, v1: float, v2: float) -> None:
        """Update statistics summary based on the path and differences."""
        max_rel_diff = max(rel_diff_v1, rel_diff_v2)
        
        # Helper function to update DiffStats if the new difference is larger
        def update_diff_stats(stats: DiffStats, abs_diff: float, rel_diff: float, path: str, v1: float, v2: float) -> None:
            if rel_diff > stats.max_rel_diff:
                stats.max_rel_diff = rel_diff
                stats.max_abs_diff = abs_diff
                stats.path = path
            # Always store the values for WM ratio
            if "wm_ratio" in path.lower():
                stats.value1 = v1
                stats.value2 = v2

        # Check if path contains b0, FA, or WM ratio statistics
        path_lower = path.lower()
        if "b0" in path_lower:
            if "mean" in path_lower:
                update_diff_stats(stats_summary.b0_mean, abs_diff, max_rel_diff, path, v1, v2)
            elif "std" in path_lower:
                update_diff_stats(stats_summary.b0_std, abs_diff, max_rel_diff, path, v1, v2)
        elif "fa" in path_lower:
            if "mean" in path_lower:
                update_diff_stats(stats_summary.fa_mean, abs_diff, max_rel_diff, path, v1, v2)
            elif "std" in path_lower:
                update_diff_stats(stats_summary.fa_std, abs_diff, max_rel_diff, path, v1, v2)
        elif "wm_ratio" in path_lower:
            update_diff_stats(stats_summary.wm_ratio, abs_diff, max_rel_diff, path, v1, v2)
    
    def compare_values(v1: Any, v2: Any, path: List[str]) -> bool:
        current_path = format_list_path(path)
        is_equal = True
        
        # Check for different types
        if type(v1) != type(v2):
            differences.append(f"Type mismatch at {current_path}: {type(v1).__name__} vs {type(v2).__name__}")
            return False
            
        # Compare based on type
        if isinstance(v1, (int, str)):
            if v1 != v2:
                differences.append(f"Value mismatch at {current_path}: {v1} vs {v2}")
                is_equal = False
        elif isinstance(v1, float) or (hasattr(v1, 'dtype') and np.issubdtype(v1.dtype, np.floating)):
            if abs(v1 - v2) >= tolerance:
                abs_diff = abs(v1 - v2)
                # Calculate relative differences, handling division by zero
                rel_diff_v1 = abs_diff / abs(v1) if v1 != 0 else float('inf')
                rel_diff_v2 = abs_diff / abs(v2) if v2 != 0 else float('inf')
                
                # Update statistics summary
                update_stats_summary(current_path, abs_diff, rel_diff_v1, rel_diff_v2, v1, v2)
                
                differences.append(
                    f"Value mismatch at {current_path}:\n"
                    f"  Values: {v1} vs {v2}\n"
                    f"  Absolute difference: {abs_diff}\n"
                    f"  Relative difference (to first value): {rel_diff_v1}\n"
                    f"  Relative difference (to second value): {rel_diff_v2}"
                )
                is_equal = False
            elif "wm_ratio" in current_path.lower():
                # For WM ratio, store the values even if they're equal
                update_stats_summary(current_path, 0.0, 0.0, 0.0, v1, v2)
        elif isinstance(v1, list):
            if len(v1) != len(v2):
                differences.append(f"List length mismatch at {current_path}: {len(v1)} vs {len(v2)}")
                is_equal = False
            # Continue comparing elements even if lengths differ
            for i, (x, y) in enumerate(zip(v1, v2)):
                is_equal &= compare_values(x, y, path + [f"[{i}]"])
            # If first list is longer, compare remaining elements with None
            for i in range(len(v2), len(v1)):
                desc = get_list_index_description(i)
                differences.append(f"Extra element in first list at {current_path}[{i}:{desc}]: {v1[i]}")
                is_equal = False
            # If second list is longer, compare remaining elements with None
            for i in range(len(v1), len(v2)):
                desc = get_list_index_description(i)
                differences.append(f"Extra element in second list at {current_path}[{i}:{desc}]: {v2[i]}")
                is_equal = False
        elif isinstance(v1, dict):
            # Get all keys while preserving order
            keys1 = list(v1.keys())
            keys2 = list(v2.keys())
            
            # Compare keys in order
            if keys1 != keys2:
                # Check for key order differences first
                if sorted(keys1) == sorted(keys2):
                    differences.append(f"Key order mismatch at {current_path}:\n"
                                    f"  First file order: {keys1}\n"
                                    f"  Second file order: {keys2}")
                    is_equal = False
                else:
                    # If keys are actually different (not just order)
                    missing_in_2 = set(keys1) - set(keys2)
                    missing_in_1 = set(keys2) - set(keys1)
                    if missing_in_2:
                        differences.append(f"Keys missing in second JSON at {current_path}: {missing_in_2}")
                        is_equal = False
                    if missing_in_1:
                        differences.append(f"Keys missing in first JSON at {current_path}: {missing_in_1}")
                        is_equal = False
            
            # Compare values for all keys that exist in both dictionaries
            # Use the order from the first file for comparison
            for k in keys1:
                if k in v2:
                    is_equal &= compare_values(v1[k], v2[k], path + [k])
            
            # Check for any keys in second file that aren't in first file
            # (we've already reported them as missing, but need to compare their values)
            for k in keys2:
                if k not in v1:
                    differences.append(f"Additional key in second JSON at {current_path}: {k} = {v2[k]}")
        
        return is_equal

    # Compare the JSONs
    are_equal = compare_values(data1, data2, [])
    
    if are_equal:
        print("JSONs are equal")
    else:
        print("JSONs are different. Differences found:")
        for diff in differences:
            print(f"- {diff}", end="\n\n")
        
        # Print summary of maximum differences
        print("\nSummary of Maximum Differences:")
        print("\nB0 Statistics:")
        if stats_summary.b0_mean.max_rel_diff > 0:
            print(f"- Mean: max relative diff = {stats_summary.b0_mean.max_rel_diff}, "
                  f"absolute diff = {stats_summary.b0_mean.max_abs_diff} "
                  f"(at {stats_summary.b0_mean.path})")
        if stats_summary.b0_std.max_rel_diff > 0:
            print(f"- Std: max relative diff = {stats_summary.b0_std.max_rel_diff}, "
                  f"absolute diff = {stats_summary.b0_std.max_abs_diff} "
                  f"(at {stats_summary.b0_std.path})")
        
        print("\nFA Statistics:")
        if stats_summary.fa_mean.max_rel_diff > 0:
            print(f"- Mean: max relative diff = {stats_summary.fa_mean.max_rel_diff}, "
                  f"absolute diff = {stats_summary.fa_mean.max_abs_diff} "
                  f"(at {stats_summary.fa_mean.path})")
        if stats_summary.fa_std.max_rel_diff > 0:
            print(f"- Std: max relative diff = {stats_summary.fa_std.max_rel_diff}, "
                  f"absolute diff = {stats_summary.fa_std.max_abs_diff} "
                  f"(at {stats_summary.fa_std.path})")

        print("\nWM Ratio:")
        if stats_summary.wm_ratio.max_rel_diff > 0:
            print(f"- Max relative diff = {stats_summary.wm_ratio.max_rel_diff}, "
                    f"absolute diff = {stats_summary.wm_ratio.max_abs_diff} "
                    f"(at {stats_summary.wm_ratio.path})")
        else:
            print(f"- Values are identical: {stats_summary.wm_ratio.value1}")

    return are_equal

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare two JSON files for equality')
    parser.add_argument('file1', type=str, help='Path to the first JSON file')
    parser.add_argument('file2', type=str, help='Path to the second JSON file')
    parser.add_argument('--tolerance', type=float, default=1e-10, 
                      help='Absolute difference tolerance for floating point comparisons (default: 1e-10)')
    
    args = parser.parse_args()
    load_and_compare_json(args.file1, args.file2, args.tolerance)
