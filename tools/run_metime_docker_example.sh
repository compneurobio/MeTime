#!/usr/bin/env bash
set -euo pipefail

IMAGE="${1:-metime:local}"

if [[ "$IMAGE" == "metime:local" ]]; then
  echo "[MeTime] Building local image from Dockerfile..."
  docker build -t metime:local .
fi

echo "[MeTime] Running package smoke check in container: $IMAGE"
docker run --rm "$IMAGE" R -q -e "library(MeTime); data('humet_object'); print(packageVersion('MeTime')); print(names(humet_object@list_of_data))"

echo "[MeTime] Done."
