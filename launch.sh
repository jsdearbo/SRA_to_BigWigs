#!/bin/bash
set -euo pipefail

echo "Step 1: Downloading missing SRA files..."
bash /mnt/netfiles/jsdearbo/public_data/pipeline_v2/download_data.sh

echo "Step 2: Processing datasets..."
bash /mnt/netfiles/jsdearbo/public_data/pipeline_v2/processing_pipeline.sh

echo "Pipeline finished successfully!"