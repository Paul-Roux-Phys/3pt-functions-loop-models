BASE_DIR=$(git rev-parse --show-toplevel)/results || exit 1

for subdir in "$BASE_DIR"/*/; do
    [ -d "$subdir" ] || continue  # skip if not a directory
    echo "Processing directory: $subdir"
    cd "$subdir" || { echo "Failed to enter $subdir"; continue; }
    for L in {4..12}; do
	echo "Merging CSVs for L=$L"
	tmpfile=$(mktemp)

	( head -n 1 "L=${L}_lambda=0.43.csv" \
	      && tail -n +2 -q L=${L}_lambda=*.csv ) > "$tmpfile" \
	    && mv "$tmpfile" "L=${L}.csv" \
		|| { echo "Failed for L=$L"; rm -f "$tmpfile"; }
    done
done
