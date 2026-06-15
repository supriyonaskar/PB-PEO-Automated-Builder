#!/bin/bash

if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    echo "======================================"
    echo " PB-PEO Test Runner Help "
    echo "======================================"
    echo "Usage: ./run_all_tests.sh"
    echo ""
    echo "This script automatically runs the generalized builder pipeline"
    echo "across 6 predefined extreme edge-case JSON configurations:"
    echo "  - Pure pb22peo14"
    echo "  - 10% Miglyol in pb22peo14"
    echo "  - 20% Cholesterol in pb22peo14"
    echo "  - Pure pb51peo27 (Dynamic topology extension)"
    echo "  - 5% Miglyol in pb51peo27"
    echo "  - 15% Cholesterol in pb100peo50"
    echo "======================================"
    exit 0
fi

echo "======================================"
echo " Running 6 Massive Automated Test Cases "
echo "======================================"

for config in examples/test_cases/*.json; do
    echo ""
    echo ">>> Running test case: $config"
    ./scripts/builder.sh "$config"
done

echo ""
echo "======================================"
echo " All 6 test cases generated successfully!"
echo "======================================"
