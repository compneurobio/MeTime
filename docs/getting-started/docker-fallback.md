# Docker Fallback Installation

If local R installation fails (system libraries, compiler issues, or dependency conflicts), use the Docker image.

## Option A: Pull from GitHub Container Registry (GHCR)

Replace `<OWNER>` with your GitHub org/user.

```bash
docker pull ghcr.io/<OWNER>/metime:latest
```

Run a quick check:

```bash
docker run --rm ghcr.io/<OWNER>/metime:latest
```

Open an interactive R session in the container:

```bash
docker run --rm -it ghcr.io/<OWNER>/metime:latest R
```

## Option B: Build locally from Dockerfile

From repository root:

```bash
docker build -t metime:local .
docker run --rm metime:local
```

## Option C: Docker Hub (optional)

Current workflow publishes to GHCR. To publish to Docker Hub as well, add Docker Hub credentials to GitHub Actions and push a second image tag (for example `docker.io/<namespace>/metime:<tag>`).

## Example helper script

Use the helper script from repo root:

```bash
bash tools/run_metime_docker_example.sh ghcr.io/<OWNER>/metime:latest
```

If no image argument is provided, the script will build `metime:local` from the local Dockerfile.
