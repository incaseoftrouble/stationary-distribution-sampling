#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

rm -f sources.tar.gz
tar -czf sources.tar.gz \
  run.py eval.py \
  src/main/java \
  lib/models/src/main/java \
  lib/models/config/template-*.txt \
  models/ \
  lib/Jeigen-onefat.jar \
  \
  build.gradle.kts settings.gradle.kts \
  lib/models/build.gradle.kts lib/models/settings.gradle.kts