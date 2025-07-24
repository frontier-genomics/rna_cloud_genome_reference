#!/bin/bash
# Comprehensive validation for Nextflow workflow

echo "🔍 Validating Nextflow workflow implementation..."
echo "=================================================="

# Check for basic syntax issues in .nf files
check_nf_syntax() {
    local file=$1
    echo "📄 Checking $file..."
    
    # Check for DSL 2 declaration
    if grep -q "nextflow.enable.dsl = 2" "$file"; then
        echo "  ✅ DSL 2 enabled"
    else
        echo "  ❌ DSL 2 not found"
        return 1
    fi
    
    # Check for process definitions
    if grep -q "^process " "$file"; then
        local process_count=$(grep -c "^process " "$file")
        echo "  ✅ Process definitions found ($process_count processes)"
    else
        echo "  ℹ️  No process definitions (may be main workflow)"
    fi
    
    # Check for workflow definitions
    if grep -q "^workflow" "$file"; then
        echo "  ✅ Workflow definition found"
    else
        echo "  ℹ️  No workflow definition (may be process-only file)"
    fi
    
    # Check for include statements
    if grep -q "^include " "$file"; then
        local include_count=$(grep -c "^include " "$file")
        echo "  ✅ Include statements found ($include_count includes)"
    else
        echo "  ℹ️  No include statements"
    fi
    
    # Basic bracket matching
    local open_braces=$(grep -o '{' "$file" | wc -l)
    local close_braces=$(grep -o '}' "$file" | wc -l)
    
    if [ "$open_braces" -eq "$close_braces" ]; then
        echo "  ✅ Braces balanced ($open_braces pairs)"
    else
        echo "  ❌ Braces unbalanced (open: $open_braces, close: $close_braces)"
        return 1
    fi
    
    # Check for publishDir directives
    if grep -q "publishDir" "$file"; then
        echo "  ✅ PublishDir directives found"
    else
        echo "  ℹ️  No publishDir directives"
    fi
    
    echo ""
    return 0
}

# Check workflow structure
check_workflow_structure() {
    echo "🏗️  Checking workflow structure..."
    
    # Check main workflow file
    if [ -f "main.nf" ]; then
        echo "  ✅ Main workflow file exists"
    else
        echo "  ❌ Main workflow file (main.nf) not found"
        return 1
    fi
    
    # Check workflows directory
    if [ -d "workflows" ]; then
        echo "  ✅ Workflows directory exists"
        local workflow_files=$(find workflows -name "*.nf" | wc -l)
        echo "  ✅ Found $workflow_files workflow files"
    else
        echo "  ❌ Workflows directory not found"
        return 1
    fi
    
    # Check configuration files
    if [ -f "nextflow.config" ]; then
        echo "  ✅ Main configuration file exists"
    else
        echo "  ❌ Main configuration file (nextflow.config) not found"
        return 1
    fi
    
    echo ""
    return 0
}

# Check configuration
check_config() {
    echo "⚙️  Checking configuration files..."
    
    if [ -f "nextflow.config" ]; then
        if grep -q "params {" nextflow.config; then
            echo "  ✅ Parameters section found"
        else
            echo "  ❌ Parameters section not found"
            return 1
        fi
        
        if grep -q "process {" nextflow.config; then
            echo "  ✅ Process configuration found"
        else
            echo "  ❌ Process configuration not found"
            return 1
        fi
        
        if grep -q "profiles {" nextflow.config; then
            echo "  ✅ Execution profiles found"
        else
            echo "  ❌ Execution profiles not found"
            return 1
        fi
        
        if grep -q "docker {" nextflow.config; then
            echo "  ✅ Docker configuration found"
        else
            echo "  ❌ Docker configuration not found"
            return 1
        fi
    fi
    
    if [ -f "test.config" ]; then
        echo "  ✅ Test configuration file exists"
    else
        echo "  ℹ️  Test configuration file not found"
    fi
    
    echo ""
    return 0
}

# Check dependencies and required files
check_dependencies() {
    echo "📦 Checking dependencies..."
    
    # Check for Python dependencies
    if [ -f "requirements.txt" ]; then
        echo "  ✅ Python requirements file exists"
    else
        echo "  ❌ Python requirements file not found"
        return 1
    fi
    
    # Check for original Python modules
    if [ -d "rnacloud_genome_reference" ]; then
        echo "  ✅ Python module directory exists"
    else
        echo "  ❌ Python module directory not found"
        return 1
    fi
    
    # Check for Dockerfile
    if [ -f "Dockerfile" ]; then
        echo "  ✅ Dockerfile exists"
    else
        echo "  ❌ Dockerfile not found"
        return 1
    fi
    
    # Check for convenience scripts
    if [ -f "run_nextflow.sh" ]; then
        echo "  ✅ Run script exists"
    else
        echo "  ❌ Run script not found"
    fi
    
    echo ""
    return 0
}

# Check documentation
check_documentation() {
    echo "📚 Checking documentation..."
    
    if [ -f "NEXTFLOW_README.md" ]; then
        echo "  ✅ Nextflow documentation exists"
    else
        echo "  ❌ Nextflow documentation not found"
    fi
    
    if [ -f "README.md" ]; then
        if grep -q -i "nextflow" README.md; then
            echo "  ✅ Main README mentions Nextflow"
        else
            echo "  ❌ Main README doesn't mention Nextflow"
        fi
    else
        echo "  ❌ Main README not found"
    fi
    
    echo ""
}

# Main validation
main() {
    local overall_status=0
    
    # Check workflow structure
    if ! check_workflow_structure; then
        overall_status=1
    fi
    
    # Check all .nf files
    for file in main.nf workflows/*.nf; do
        if [ -f "$file" ]; then
            if ! check_nf_syntax "$file"; then
                overall_status=1
            fi
        fi
    done
    
    # Check configuration
    if ! check_config; then
        overall_status=1
    fi
    
    # Check dependencies
    if ! check_dependencies; then
        overall_status=1
    fi
    
    # Check documentation
    check_documentation
    
    echo "=================================================="
    if [ $overall_status -eq 0 ]; then
        echo "🎉 Validation completed successfully!"
        echo "✅ Nextflow workflow appears to be properly structured."
    else
        echo "❌ Validation failed!"
        echo "🔧 Please fix the issues above before proceeding."
    fi
    
    echo ""
    echo "💡 Next steps:"
    echo "   1. Install Nextflow: curl -s https://get.nextflow.io | bash"
    echo "   2. Build Docker container: docker build -t rnacloud_runner ."
    echo "   3. Test workflow: nextflow run main.nf -c test.config --profile test -entry DOWNLOAD_ANNOTATION"
    echo "   4. Run full workflow: ./run_nextflow.sh"
    
    return $overall_status
}

# Run main function
main