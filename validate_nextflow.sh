#!/bin/bash
# Simple syntax validation for Nextflow files

echo "Validating Nextflow workflow syntax..."

# Check for basic syntax issues in .nf files
check_nf_syntax() {
    local file=$1
    echo "Checking $file..."
    
    # Check for DSL 2 declaration
    if grep -q "nextflow.enable.dsl = 2" "$file"; then
        echo "  ✓ DSL 2 enabled"
    else
        echo "  ✗ DSL 2 not found"
    fi
    
    # Check for process definitions
    if grep -q "^process " "$file"; then
        echo "  ✓ Process definitions found"
    else
        echo "  ✓ No process definitions (may be main workflow)"
    fi
    
    # Check for workflow definitions
    if grep -q "^workflow" "$file"; then
        echo "  ✓ Workflow definition found"
    else
        echo "  ✓ No workflow definition (may be process-only file)"
    fi
    
    # Basic bracket matching
    local open_braces=$(grep -o '{' "$file" | wc -l)
    local close_braces=$(grep -o '}' "$file" | wc -l)
    
    if [ "$open_braces" -eq "$close_braces" ]; then
        echo "  ✓ Braces balanced ($open_braces pairs)"
    else
        echo "  ✗ Braces unbalanced (open: $open_braces, close: $close_braces)"
    fi
    
    echo ""
}

# Check all .nf files
for file in main.nf workflows/*.nf; do
    if [ -f "$file" ]; then
        check_nf_syntax "$file"
    fi
done

# Check config file
if [ -f "nextflow.config" ]; then
    echo "Checking nextflow.config..."
    if grep -q "params {" nextflow.config; then
        echo "  ✓ Parameters section found"
    fi
    if grep -q "process {" nextflow.config; then
        echo "  ✓ Process configuration found"
    fi
    if grep -q "profiles {" nextflow.config; then
        echo "  ✓ Execution profiles found"
    fi
    echo ""
fi

echo "Basic syntax validation completed!"
echo "Note: This is a basic check. Full validation requires Nextflow installation."