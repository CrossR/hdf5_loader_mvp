import sys

import awkward as ak
import numpy as np
import uproot


def main():
    if len(sys.argv) != 3:
        print("Usage: python compare_trees.py <cpp_output.root> <macro_output.root>")
        sys.exit(1)

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]
    tree_name = "events"

    print(f"Comparing '{file1_path}' against '{file2_path}'...")

    with uproot.open(file1_path) as f1, uproot.open(file2_path) as f2:
        t1 = f1[tree_name]
        t2 = f2[tree_name]
        common_keys = sorted(list(set(t1.keys()).intersection(set(t2.keys()))))

        all_match = True
        for key in common_keys:
            arr1 = t1[key].array()
            arr2 = t2[key].array()

            if len(arr1) != len(arr2):
                print(
                    f"❌ [FAIL] {key}: Event count mismatch! ({len(arr1)} vs {len(arr2)})"
                )
                all_match = False
                continue

            try:
                if arr1.ndim > 1:
                    # Check outer dimension (Number of Hits)
                    hits_match = ak.all(ak.num(arr1, axis=1) == ak.num(arr2, axis=1))
                    if not hits_match:
                        print(
                            f"❌ [FAIL] {key}: Number of hits in array differs! (C++ has {ak.sum(ak.num(arr1, axis=1))}, Python has {ak.sum(ak.num(arr2, axis=1))})"
                        )
                        all_match = False
                        continue

                    if arr1.ndim > 2:
                        # Check inner dimension (Number of Matches per Hit)
                        matches_match = ak.all(
                            ak.num(arr1, axis=2) == ak.num(arr2, axis=2)
                        )
                        if not matches_match:
                            print(
                                f"❌ [FAIL] {key}: Number of truth matches per hit differs!"
                            )
                            all_match = False
                            continue

                        # Sort to guarantee order-invariance
                        arr1 = ak.sort(arr1, axis=-1)
                        arr2 = ak.sort(arr2, axis=-1)

                is_equal = ak.all(ak.isclose(arr1, arr2, equal_nan=True))
                if is_equal:
                    num_events = len(arr1)
                    num_hits = ak.sum(ak.num(arr1, axis=1)) if arr1.ndim > 1 else 0
                    print(f"✅ [PASS] {key}: Both equal with {num_events} events and {num_hits} hits.")
                else:
                    print(f"❌ [FAIL] {key}: Value mismatch detected after sorting!")
                    all_match = False

            except Exception as e:
                print(f"⚠️ [WARN] {key}: Could not perform direct comparison ({e})")
                all_match = False

        print("-" * 40)
        if all_match:
            print("✅ SUCCESS: All common branches match perfectly!")
        else:
            print(
                "❌ FAILURE: Differences were found (Likely proving the Python truncation bug)."
            )


if __name__ == "__main__":
    main()
