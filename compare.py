import sys
import uproot
import awkward as ak
import numpy as np

def main():
    if len(sys.argv) != 3:
        print("Usage: python compare_trees.py <cpp_output.root> <macro_output.root>")
        sys.exit(1)

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]
    tree_name = "events"

    print(f"Comparing '{file1_path}' against '{file2_path}'...")

    with uproot.open(file1_path) as f1, uproot.open(file2_path) as f2:
        if tree_name not in f1 or tree_name not in f2:
            print(f"Error: '{tree_name}' tree not found in one or both files.")
            sys.exit(1)

        t1 = f1[tree_name]
        t2 = f2[tree_name]

        keys1 = set(t1.keys())
        keys2 = set(t2.keys())
        common_keys = sorted(list(keys1.intersection(keys2)))

        print(f"Found {len(common_keys)} common branches. Running diff...")

        all_match = True
        for key in common_keys:
            arr1 = t1[key].array()
            arr2 = t2[key].array()

            # 1. Check Top-Level Event Counts
            if len(arr1) != len(arr2):
                print(f"❌ [FAIL] {key}: Event count mismatch! ({len(arr1)} vs {len(arr2)})")
                all_match = False
                continue

            try:
                # 2. Handle Nested/Jagged Arrays (1-to-N)
                if arr1.ndim > 1:
                    # First, ensure the number of inner items matches exactly
                    lengths_match = ak.all(ak.num(arr1, axis=1) == ak.num(arr2, axis=1))
                    if not lengths_match:
                        print(f"❌ [FAIL] {key}: Jagged array shapes (inner list lengths) do not match.")
                        all_match = False
                        continue

                    # Sort the inner lists so [11, 22] matches [22, 11]
                    arr1 = ak.sort(arr1, axis=-1)
                    arr2 = ak.sort(arr2, axis=-1)

                # 3. Check Exact Values (with float tolerance)
                is_equal = ak.all(ak.isclose(arr1, arr2, equal_nan=True))

                if is_equal:
                    print(f"✅ [PASS] {key}")
                else:
                    print(f"❌ [FAIL] {key}: Value mismatch detected after sorting!")
                    all_match = False

            except Exception as e:
                # Fallback for complex structural differences
                print(f"⚠️ [WARN] {key}: Could not perform direct comparison ({e})")
                all_match = False

        print("-" * 40)
        if all_match:
            print("✅ SUCCESS: All common branches match perfectly!")
        else:
            print("❌ FAILURE: Differences were found.")

if __name__ == "__main__":
    main()
