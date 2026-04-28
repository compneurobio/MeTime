# Docker Fallback Installation

If local R installation fails (system libraries, compiler issues, or dependency conflicts), use the Docker image.

## Option A: Pull from GitHub Container Registry (GHCR)

```bash
docker pull ghcr.io/compneurobio/metime:latest
```

Run a quick check:

```bash
docker run --rm ghcr.io/compneurobio/metime:latest
```

Open an interactive R session in the container:

```bash
docker run --rm -it ghcr.io/compneurobio/metime:latest R
```

## Option B: Build locally from Dockerfile

Clone the repository locally and from repository root:

```bash
docker build -t metime:local .
docker run --rm metime:local
```

## Example helper script

Use the helper script from repo root:

```bash
bash tools/run_metime_docker_example.sh ghcr.io/compneurobio/metime:latest
```

If no image argument is provided, the script will build `metime:local` from the local Dockerfile.
