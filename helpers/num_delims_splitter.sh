#!/bin/bash

# Check arguments
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <file>"
    echo "Work on tab-delimited files.  Takes <file> and splits into subfiles"
    echo "depending on number of tabs per line."
    echo " data_file.ext --> data_file_tab5.ext, data_file_tab_8.ext, etc."

    # echo "Delimiter should be ',' or 'TAB'"
    # echo "Set remove_quotes to FALSE if you don't want to count delimiters within double quotes"
    # echo "e.g. for Ted,2018-06-54,\"loves cars, dogs, spoons\",ABC"
    # echo "  remove_quotes=TRUE gives a delimiter count of 5"
    # echo "  remove_quotes=FALSE gives a delimiter count of 3"
    exit 1
fi

file="$1"
directory=$(dirname "$file")
base_name=$(basename "$file")
ext="${base_name##*.}"
base_name="${base_name%.*}"
# delimiter="$2"

# Check if the file has an extension
if [[ "$ext" == "$base_name" ]]; then
    ext=""  # No extension found
else
    ext=".$ext"  # Ensure the dot is included if an extension exists
fi

# Process file with awk
awk -F'\t' -v base="$base_name" -v ext="$ext" -v path="$directory" '
NF {
    filename = path "/" base "_tab" (NF-1) ext;
    print >> filename;
}' "$file"

echo "Splitting complete. Files created: ${base_name}_tab*${ext}"
