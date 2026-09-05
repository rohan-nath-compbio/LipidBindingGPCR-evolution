#!/usr/bin/bash
set -euo pipefail

# Configuration
# Set by script.sh or provide it when calling this
# script directly: GENE_NEIGHBORHOOD_DB=/path/to/gene_neighborhood.tsv.
: "${GENE_NEIGHBORHOOD_DB:?Set GENE_NEIGHBORHOOD_DB to the gene-neighborhood database path.}"
DB="$GENE_NEIGHBORHOOD_DB"

if [[ ! -s "$DB" ]]; then
    echo "Error: gene-neighborhood database is missing or empty: $DB" >&2
    exit 1
fi

TAXA=("Mammalia" "Archosauria" "Lepidosauria" "Testudines" "Amphibia" "Dipnomorpha" "Coelacanthiformes" "Actinopterygii" "Chondrichthyes" "Cyclostomata")

# Cleanup function
cleanup() {
    rm -f *.gene.neighborhood
}

# Remove empty files
remove_empty_files() {
    local pattern="$1"
    local count=0
    for file in $pattern; do
        [[ ! -f "$file" ]] && continue
        if [[ ! -s "$file" ]]; then
            rm -f "$file"
            count=$((count + 1))
        fi
    done
    if [[ $count -gt 0 ]]; then
        echo "    Removed $count empty files"
    fi
}

# Process 1: Extract accessions by taxa
extract_accessions() {
    echo "Step 1: Extracting accessions by taxa..."
    
    local found_files=0
    local created_count=0
    local skipped_count=0
    
    for lineage_file in *.acc.taxid.lineage; do
        # Check if glob matched any files
        [[ ! -f "$lineage_file" ]] && continue
        
        found_files=1
        echo "  Processing lineage file: $lineage_file"
        
        name="${lineage_file%%.*}"
        
        # Check if ALL acc files already exist for this gene
        local all_exist=1
        for taxa in "${TAXA[@]}"; do
            lower_taxa=$(echo "$taxa" | tr '[:upper:]' '[:lower:]')
            output_file="${name}.${lower_taxa}.acc"
            if [[ ! -f "$output_file" ]] || [[ ! -s "$output_file" ]]; then
                all_exist=0
                break
            fi
        done
        
        # Also check invertebrate file
        if [[ ! -f "${name}.invertebrate.acc" ]] || [[ ! -s "${name}.invertebrate.acc" ]]; then
            all_exist=0
        fi
        
        # If all files exist, skip this entire gene
        if [[ $all_exist -eq 1 ]]; then
            echo "    All acc files exist for $name, skipping..."
            continue
        fi
        
        # Extract by taxa
        for taxa in "${TAXA[@]}"; do
            lower_taxa=$(echo "$taxa" | tr '[:upper:]' '[:lower:]')
            output_file="${name}.${lower_taxa}.acc"
            
            # Skip if file already exists and is not empty
            if [[ -f "$output_file" ]] && [[ -s "$output_file" ]]; then
                skipped_count=$((skipped_count + 1))
                continue
            fi
            
            # Extract matching lines
            grep -E "$taxa" "$lineage_file" 2>/dev/null | awk '{print $1}' > "$output_file" || true
            
            if [[ -s "$output_file" ]]; then
                echo "    Created: $output_file ($(wc -l < "$output_file") entries)"
                created_count=$((created_count + 1))
            else
                rm -f "$output_file"
            fi
        done
        
        # Extract invertebrates
        output_file="${name}.invertebrate.acc"
        
        # Skip if file already exists and is not empty
        if [[ -f "$output_file" ]] && [[ -s "$output_file" ]]; then
            skipped_count=$((skipped_count + 1))
        else
            # Extract non-vertebrate lines
            grep -v "Vertebrata" "$lineage_file" 2>/dev/null | awk '{print $1}' > "$output_file" || true
            
            if [[ -s "$output_file" ]]; then
                echo "    Created: $output_file ($(wc -l < "$output_file") entries)"
                created_count=$((created_count + 1))
            else
                rm -f "$output_file"
            fi
        fi
    done
    
    # Clean up any empty .acc files that might have been created
    echo "  Cleaning up empty .acc files..."
    remove_empty_files "*.acc"
    
    if [[ $found_files -eq 0 ]]; then
        echo "  WARNING: No *.acc.taxid.lineage files found!"
    else
        echo "  Summary: Created $created_count files, skipped $skipped_count existing files"
    fi
}

# Process 2: Generate gene neighborhoods (with skip if exists)
generate_neighborhoods() {
    local acc_file="$1"
    local neighbor_file="${acc_file%.acc}.neighbor"
    
    # Skip if neighbor file already exists
    if [[ -f "$neighbor_file" ]] && [[ -s "$neighbor_file" ]]; then
        return
    fi
    
    local index=1
    
    while read -r line; do
        [[ -z "$line" ]] && continue
        
        local neighborhood_file="${acc_file}.${index}.gene.neighborhood"
        local base_acc="${line%.*}"
        
        echo "$line" >> "$neighborhood_file"
        # Use || true to prevent exit on fgrep failure
        fgrep -w -m 1 -A 5 -B 5 "$base_acc" "$DB" 2>/dev/null \
            | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' \
            >> "$neighborhood_file" || true
        
        index=$((index + 1))
    done < "$acc_file"
    
    # Combine neighborhoods
    if ls "${acc_file}."*.gene.neighborhood 1> /dev/null 2>&1; then
        cat "${acc_file}."*.gene.neighborhood \
            | sed 's/ /\t/g' \
            | awk -v i=9 'NR>1 && $i!=p { print "" } { p=$i } 1' \
            > "$neighbor_file"
    fi
}

# Process 3: Clean neighbor files with proper formatting
deduplicate_sort_clean() {
    local neighbor_file="$1"
    local clean_file="${neighbor_file}.clean"
    
    # Skip if clean file already exists and is not empty
    if [[ -f "$clean_file" ]] && [[ -s "$clean_file" ]]; then
        return
    fi
    
    # Check if neighbor file has actual data (more than just accession lines)
    local data_lines=$(awk -F'\t' 'NF > 1' "$neighbor_file" | wc -l)
    
    if [[ $data_lines -eq 0 ]]; then
        echo "    WARNING: No neighborhood data found in $neighbor_file, skipping..."
        return
    fi
    
    # Process the file to clean and format the output
    awk -F'\t' '
    NR > 1 && NF != 1 {
        if ($0 != "") {
            if (prev_blank) print ""
            print $0
            prev_blank = 0
        } else {
            prev_blank = 1
        }
    }
    END {
        if (prev_blank) print ""
    }
    ' "$neighbor_file" | sed '1d' | grep -v '^$' | sort -k10,10 -k6,6n -k1,1n -k2,2n | uniq | awk '
    BEGIN {FS=OFS="\t"}
    {
        if (prev != "" && $10 != prev) print "#"
        print
        prev = $10
    }
    ' > "$clean_file"
    
    # Remove clean file if it's empty
    if [[ ! -s "$clean_file" ]]; then
        rm -f "$clean_file"
    fi
}

# Process 4: Mark query sequences
mark_query_sequences() {
    local acc_file="$1"
    local clean_file="${acc_file%.acc}.neighbor.clean"
    
    [[ ! -f "$clean_file" ]] && return
    [[ ! -s "$clean_file" ]] && return
    
    # Remove existing markers (idempotency)
    sed -i 's/^==> //g' "$clean_file"
    
    # Mark query sequences
    while IFS= read -r accession; do
        base_acc="${accession%.*}"
        sed -i "/${base_acc}\./s/^/==> /" "$clean_file" || true
    done < "$acc_file"
}

# Process 5: Filter unmarked sections
filter_unmarked_sections() {
    local clean_file="$1"
    
    [[ ! -f "$clean_file" ]] && return
    [[ ! -s "$clean_file" ]] && return
    
    local temp_file="${clean_file}.filtered"
    
    awk '
    BEGIN {
        keep = 0
        buffer = ""
    }
    
    /^#/ {
        # If buffer is not empty and keep is true, write buffered content
        if (buffer != "" && keep == 1) {
            print buffer
        }
        # Reset buffer and keep flag
        buffer = ""
        keep = 0
        print $0
        next
    }
    
    /^==>/ {
        keep = 1
        buffer = buffer $0 "\n"
        next
    }
    
    {
        buffer = buffer $0 "\n"
    }
    
    END {
        # Check final buffer
        if (buffer != "" && keep == 1) {
            printf "%s", buffer
        }
    }
    ' "$clean_file" | grep -v "^$" > "$temp_file"
    
    mv "$temp_file" "$clean_file"
}

# Process 6: Collapse consecutive hash marks
collapse_consecutive_hashes() {
    local clean_file="$1"
    
    [[ ! -f "$clean_file" ]] && return
    [[ ! -s "$clean_file" ]] && return
    
    local temp_file="${clean_file}.collapsed"
    
    awk '
    BEGIN {
        prev_was_hash = 0
    }
    
    /^#/ {
        if (prev_was_hash == 0) {
            print $0
        }
        prev_was_hash = 1
        next
    }
    
    {
        print $0
        prev_was_hash = 0
    }
    ' "$clean_file" > "$temp_file"
    
    mv "$temp_file" "$clean_file"
}

# Process 7: Arrange by taxa
arrange_by_taxa() {
    echo "Step 7: Arranging results by taxa..."
    
    local arranged_count=0
    local skipped_count=0
    
    for file in *.neighbor.clean; do
        # Skip if glob didn't match
        [[ ! -f "$file" ]] && continue
        
        # Extract the base name without extension
        name="${file%%.*}"
        
        # Create a new file for the arranged content
        arranged_file="${name}.arranged.neighbor"
        
        # Skip if arranged file already exists and is not empty
        if [[ -f "$arranged_file" ]] && [[ -s "$arranged_file" ]]; then
            echo "  Skipping: $arranged_file (already exists)"
            skipped_count=$((skipped_count + 1))
            continue
        fi
        
        # Start with an empty arranged file
        > "$arranged_file"
        
        # Loop through each taxa and append its content to the arranged file
        for t in "${TAXA[@]}" "Invertebrate"; do
            lower_taxa=$(echo "$t" | tr '[:upper:]' '[:lower:]')
            # Construct the filename for the current taxa
            taxa_file="${name}.${lower_taxa}.neighbor.clean"
            
            # Check if the taxa file exists and is not empty
            if [[ -f "$taxa_file" ]] && [[ -s "$taxa_file" ]]; then
                # Write a clear and simple taxa header
                echo "---- ${t^^} ----" >> "$arranged_file"
                echo >> "$arranged_file"  # Add an empty row after the header
                
                # Append the content of the taxa file
                cat "$taxa_file" >> "$arranged_file"
                echo >> "$arranged_file"  # Add an extra empty row for separation
            fi
        done
        
        # Check if arranged file has content
        if [[ -s "$arranged_file" ]]; then
            echo "  Created: $arranged_file"
            arranged_count=$((arranged_count + 1))
        else
            rm -f "$arranged_file"
        fi
    done
    
    echo "  Summary: Created $arranged_count arranged files, skipped $skipped_count existing files"
}

# Main pipeline
main() {
    echo "Starting synteny analysis pipeline..."
    
    # Step 1: Extract accessions
    extract_accessions
    
    # Steps 2-6: Process each acc file
    echo "Step 2-6: Processing gene neighborhoods..."
    
    local acc_count=0
    local neighbor_created=0
    local neighbor_skipped=0
    local clean_created=0
    local clean_skipped=0
    
    for acc_file in *.acc; do
        # Skip if file doesn't exist (no glob match)
        [[ ! -f "$acc_file" ]] && continue
        
        # Skip empty files
        [[ ! -s "$acc_file" ]] && continue
        
        # Skip lineage files
        [[ "$acc_file" == *.taxid.lineage ]] && continue
        
        acc_count=$((acc_count + 1))
        
        echo "  Processing: $acc_file"
        
        neighbor_file="${acc_file%.acc}.neighbor"
        clean_file="${neighbor_file}.clean"
        
        # Check if final arranged file exists for this gene
        name="${acc_file%%.*}"
        arranged_file="${name}.arranged.neighbor"
        
        if [[ -f "$arranged_file" ]] && [[ -s "$arranged_file" ]]; then
            echo "    Final arranged file exists, skipping all processing"
            continue
        fi
        
        # Generate neighborhoods (skip if exists)
        if [[ -f "$neighbor_file" ]] && [[ -s "$neighbor_file" ]]; then
            echo "    Neighbor file exists, skipping generation"
            neighbor_skipped=$((neighbor_skipped + 1))
        else
            echo "    Generating neighborhoods..."
            cleanup
            generate_neighborhoods "$acc_file"
            if [[ -f "$neighbor_file" ]] && [[ -s "$neighbor_file" ]]; then
                neighbor_created=$((neighbor_created + 1))
            fi
        fi
        
        # Continue processing only if neighbor file exists
        if [[ -f "$neighbor_file" ]] && [[ -s "$neighbor_file" ]]; then
            # Check if clean file already exists
            if [[ -f "$clean_file" ]] && [[ -s "$clean_file" ]]; then
                echo "    Clean file exists, skipping processing"
                clean_skipped=$((clean_skipped + 1))
            else
                echo "    Deduplicating and cleaning..."
                deduplicate_sort_clean "$neighbor_file"
                
                if [[ -f "$clean_file" ]] && [[ -s "$clean_file" ]]; then
                    echo "    Marking query sequences..."
                    mark_query_sequences "$acc_file"
                    
                    echo "    Filtering unmarked sections..."
                    filter_unmarked_sections "$clean_file"
                    
                    echo "    Collapsing hash markers..."
                    collapse_consecutive_hashes "$clean_file"
                    
                    clean_created=$((clean_created + 1))
                fi
            fi
        else
            echo "    WARNING: No valid neighbor file for $acc_file"
        fi
    done
    
    # Clean up any empty files
    echo "  Cleaning up empty files..."
    remove_empty_files "*.neighbor"
    remove_empty_files "*.neighbor.clean"
    
    echo "  Summary: Processed $acc_count .acc files"
    echo "           Neighbor files - Created: $neighbor_created, Skipped: $neighbor_skipped"
    echo "           Clean files - Created: $clean_created, Skipped: $clean_skipped"
    
    # Final cleanup
    cleanup
    
    # Step 7: Arrange by taxa
    arrange_by_taxa
    
    echo "Pipeline complete!"
}

main "$@"
