#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Define paths
json_file="version_switcher.json"
docs_path="docs"
base_url="https://martinpdes.github.io/SuPyMode"

# Check if docs_path exists
if [ ! -d "$docs_path" ]; then
  echo "Error: $docs_path does not exist."
  exit 1
fi

# Get sorted version directories in descending order
version_dirs=$(ls -1v "$docs_path" | tac)

# Generate the JSON array using jq
jq -n --arg base_url "$base_url" --arg docs_path "$docs_path" '
  [inputs | {name: ., version: ., url: ($base_url + "/" + $docs_path + "/" + .)}]
' <<< "$version_dirs" > "$json_file"

# Add and commit changes
git add "$json_file"
git commit --allow-empty -m "Update version_switcher.json"
git push origin HEAD

