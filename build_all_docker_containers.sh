#!/bin/bash

set -e

# Build all components in a namespace (refer https://viash.io/reference/cli/ns_build.html)
# and set up the container via a cached build
#viash ns build --parallel --setup cachedbuild
viash ns build --parallel --setup cachedbuild -n "normalization"
viash ns build --parallel --setup cachedbuild -n "processors"
viash ns build --parallel --setup cachedbuild -q "extract_dataset_meta"
viash ns build --parallel --setup cachedbuild -q "local_h5ad"

