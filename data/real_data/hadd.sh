
#!/bin/bash

# Configuration
BASE_DIR="/sps/nemo/scratch/granjon/double_beta/cimrman_curve/ptd_pol2/root_with_max_threshold_working"
MIN_RUN=3470
MAX_RUN=5448
CHUNK_SIZE=100

# Find existing ROOT files
FILES=()
MAX_FOUND=0

for ((run=MIN_RUN; run<=MAX_RUN; run++)); do
    FILE="${BASE_DIR}/run_${run}/extracted_data.root"

    if [[ -f "$FILE" ]]; then
        FILES+=("$FILE")
        MAX_FOUND=$run
        echo "Found: run_${run}"
    fi
done

# Check if any files were found
if [[ ${#FILES[@]} -eq 0 ]]; then
    echo "ERROR: No ROOT files found!"
    exit 1
fi

OUTPUT="run_${MIN_RUN}-${MAX_FOUND}.root"
TEMP_DIR="hadd_temp"

mkdir -p "$TEMP_DIR"

echo ""
echo "Found ${#FILES[@]} ROOT files."
echo "Merging in chunks of ${CHUNK_SIZE}..."

# Merge files in chunks
TEMP_FILES=()
CHUNK=0

for ((i=0; i<${#FILES[@]}; i+=CHUNK_SIZE)); do
    CHUNK_FILES=("${FILES[@]:i:CHUNK_SIZE}")
    TEMP_FILE="${TEMP_DIR}/chunk_${CHUNK}.root"

    echo "Merging chunk ${CHUNK} (${#CHUNK_FILES[@]} files)..."

    hadd -k -f "$TEMP_FILE" "${CHUNK_FILES[@]}"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Failed merging chunk ${CHUNK}!"
        exit 1
    fi

    TEMP_FILES+=("$TEMP_FILE")
    ((CHUNK++))
done

# Final merge
echo ""
echo "Merging ${#TEMP_FILES[@]} temporary files into ${OUTPUT}..."

hadd -k -f "$OUTPUT" "${TEMP_FILES[@]}"

if [[ $? -eq 0 ]]; then
    echo "Successfully created: ${OUTPUT}"
    rm -rf "$TEMP_DIR"
    echo "Temporary files removed."
else
    echo "ERROR: Final merge failed! Temporary files preserved in ${TEMP_DIR}/"
    exit 1
fi
