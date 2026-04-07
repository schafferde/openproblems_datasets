#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

echo "Running data preprocessing"
echo "  Make sure to run 'scripts/project/build_all_docker_containers.sh'!"

# generate a unique id
RUN_ID="datarun_$(date +%Y-%m-%d_%H-%M-%S)"
publish_dir="resources/datasets/${RUN_ID}"

# write the parameters to file
cat > /tmp/params.yaml << HERE
input_states: /local/dschaffe/openproblems_datasets/celegans_config.yaml
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

nextflow run . \
  -main-script target/nextflow/workflows/scrnaseq/process_local_h5ad/main.nf \
  -profile docker \
  -entry auto \
  -c labels_copy.config \
  -params-file /tmp/params.yaml \
  -resume
