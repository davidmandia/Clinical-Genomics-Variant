#!/bin/bash

# Define the database paths
DB_PATHS=("output/database/GRCh38_indels_variant.db" "output/database/GRCh37_indels_variant.db")

# Define the directory containing the analysis scripts
ANALYSIS_DIR="analysis"

# Iterate over each database path
for db_path in "${DB_PATHS[@]}"; do
    echo "Running analysis scripts for database: $db_path"
    # Iterate over each Python script in the analysis directory
    for script in "$ANALYSIS_DIR"/*.py; do
        echo "Executing script: $script with database: $db_path"
        python3 "$script" "$db_path"
    done
done
