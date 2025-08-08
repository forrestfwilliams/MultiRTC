#!/bin/bash --login
set -e
conda activate multirtc
exec python -um multirtc "$@"
