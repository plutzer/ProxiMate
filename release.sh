#!/usr/bin/env bash
# Build and publish the ProxiMate container for one release tag from a local
# machine.  Produces the same image the Release workflow builds: the tag is
# stamped in as PROXIMATE_VERSION and the annotation databases come from the
# Datasets/ directory, downloaded first if missing.
#
#   ./release.sh v0.2.0                       # build amd64+arm64, push to Docker Hub
#   ./release.sh v0.2.0 --platforms linux/amd64
#   ./release.sh v0.2.0 --no-push             # amd64 only, loaded into local Docker
#
# Requires Docker Desktop (buildx) and `docker login`.  Set GHCR=1 to also push
# to ghcr.io/plutzer/proximate after `docker login ghcr.io`.  The working tree
# must be clean and checked out at the tag, so the image matches the release.
set -euo pipefail

HUB_IMAGE=plutzer/proximate
GHCR_IMAGE=ghcr.io/plutzer/proximate

tag=${1:?usage: release.sh vX.Y.Z [--platforms list] [--no-push]}
shift
platforms=linux/amd64,linux/arm64
push=1
while [ $# -gt 0 ]; do
    case $1 in
        --platforms) platforms=$2; shift 2 ;;
        --no-push) push=0; platforms=linux/amd64; shift ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

cd "$(dirname "$0")"

[[ $tag =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]] || { echo "tag must look like v1.2.3, got $tag" >&2; exit 2; }
git rev-parse -q --verify "refs/tags/$tag" >/dev/null || { echo "tag $tag does not exist; create the release on GitHub first, then fetch" >&2; exit 2; }
[ -z "$(git status --porcelain)" ] || { echo "working tree is not clean" >&2; exit 2; }
[ "$(git rev-parse HEAD)" = "$(git rev-parse "$tag^{commit}")" ] || { echo "HEAD is not at $tag; check it out first" >&2; exit 2; }

if [ ! -s Datasets/build_info.txt ]; then
    echo "downloading annotation databases into Datasets/"
    docker run --rm -v "$(pwd -W 2>/dev/null || pwd):/work" -w /work python:3.12-slim \
        sh -c "pip install -q requests && python3 Scripts/setup_datasets.py --output-dir Datasets --skip corum"
fi
grep -q FAILED Datasets/build_info.txt && { echo "a database download failed; see Datasets/build_info.txt" >&2; exit 1; }

tags=(-t "$HUB_IMAGE:$tag" -t "$HUB_IMAGE:latest")
if [ "${GHCR:-0}" = 1 ]; then
    tags+=(-t "$GHCR_IMAGE:$tag" -t "$GHCR_IMAGE:latest")
fi
if [ $push = 1 ]; then
    output=--push
else
    output=--load
fi

docker buildx build --platform "$platforms" --build-arg "PROXIMATE_VERSION=$tag" \
    "${tags[@]}" $output .

if [ $push = 1 ]; then
    docker buildx imagetools inspect "$HUB_IMAGE:$tag"
fi
