#!/bin/bash

# --- CONFIGURATION ---
# WIN_BASE is where saved files go this is my local file path as an example
WIN_BASE="/mnt/c/Users/Sam/Phylogenetic-Tree-Building-Analysis/ROSE"

# The four templates [cite: 14]
MODELS="JC_low JC_high HKY_low HKY_high"

# Lists of taxa and sequence length vareity desired
TAXA_LIST="10 20 35 50 75 100 250 500"
LENGTH_LIST="10 50 100 250 500 1000 2500 5000"

# Change to generate multiple files per run, must update name to accept runs
RUNS=1

for MODEL in $MODELS
do
    # Define the subfolder paths for this specific model
    MODEL_DIR="${WIN_BASE}/to_calc"
    PHY_DIR="${MODEL_DIR}/phylip"
    FAS_DIR="${MODEL_DIR}/fasta"
    TR_DIR="${MODEL_DIR}/true_trees"

    # Create the structure if it doesn't exist
    mkdir -p "$PHY_DIR"
    mkdir -p "$TR_DIR"
    mkdir -p "$FAS_DIR"

    for LEN in $LENGTH_LIST
    do
        for TAXA in $TAXA_LIST
        do
            for RUN in $(seq 1 $RUNS)
            do
                # Unique name identifying the simulation
                NAME="${MODEL}_${TAXA}_${LEN}"
                #
                echo "Processing $MODEL: Taxa=$TAXA, Length=$LEN"

                # Calculate Relatedness based on Taxa count, these numbers worked
		# There are probably better options
                if [[ $MODEL == *"high"* ]]; then
                    REL=$((50 + TAXA * 5 / 10))
                else
                    REL=$((TAXA * 4 / 10))
                fi

                #
                (cat "$MODEL";
                 echo "SequenceNum = $TAXA;";
                 echo "SequenceLen = $LEN;";
                 echo "Relatedness = $REL;";
                 echo "SequenceOutputLen = $((LEN * 3));";
                 echo "OutputFilebase = \"$NAME\";") | rose -

                # Move files to their specific sub-folders
                [ -f "${NAME}.phy" ] && mv "${NAME}.phy" "$PHY_DIR/"
                [ -f "${NAME}.tree" ] && mv "${NAME}.tree" "$TR_DIR/"
		[ -f "${NAME}.fas" ] && mv "${NAME}.fas" "$FAS_DIR/"
            done
        done
    done
done

echo "Workflow complete. Files are sorted by Model -> File Type."
