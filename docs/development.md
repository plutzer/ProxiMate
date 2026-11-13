# Development

## Repository layout

| Path | Contents |
| --- | --- |
| `GUI/` | Shiny application, MCP server and the shared dataset backend |
| `Scripts/` | Parsing, scoring, imputation and annotation pipeline; `Scripts/tests/` is the pytest suite |
| `Scripts/GOGO/` | GO semantic similarity tool (Zhao and Wang 2018) |
| `saint/upstream/` | SAINTexpress 3.6.3 as distributed, with the Boost and NLopt sources it builds against |
| `saint/patches/` | ProxiMate's modified SAINTexpress intensity-model sources, overlaid on the upstream tree at image build |
| `Datasets/` | CORUM complexes; the other annotation databases are downloaded at build time |
| `examples/` | Public example datasets |
| `docs/` | Tutorials and reference pages |

## Build the image

```
docker build --build-arg PROXIMATE_VERSION=$(git rev-parse --short HEAD) \
  -t plutzer/proximate:latest .
```

The image contains no git repository. `PROXIMATE_VERSION` is the only source of
the version recorded in `run.json`; an unstamped build reports `unknown`.

## Run locally and test

```
pip install -r requirements.txt   # Python 3.11+
python GUI/app.py
pytest
```

The test suite needs neither the SAINT binary nor the databases.

## Pin a database snapshot

The Dockerfile copies `Datasets/` into the image, then runs
`python3 Scripts/setup_datasets.py --output-dir /Datasets --skip corum`, which
downloads only the files that are missing. So a file placed in `Datasets/` before
`docker build` takes precedence over the download.

To fetch the current set locally:

```
python3 Scripts/setup_datasets.py --output-dir Datasets
```

`--help` lists `--skip` and the other options. File names and sources are in
[Databases](databases.md#sources).

## Releases

Versions are `vMAJOR.MINOR.PATCH`. Images are published for `linux/amd64` and
`linux/arm64` to two registries; either name pulls the same image:

```
docker pull plutzer/proximate:v0.4.0          # Docker Hub
docker pull ghcr.io/plutzer/proximate:v0.4.0  # GitHub Container Registry
```

`latest` follows the newest release.

**To release:** publish a GitHub Release with the tag. That runs
`.github/workflows/release.yml`, which:

1. downloads the databases once,
2. builds both architectures natively from that snapshot,
3. joins them under the tag and `latest` on both registries,
4. appends the image digest and `build_info.txt` to the release notes.

**Without GitHub Actions**, `release.sh` runs the same build from a local checkout
of the tag:

```
./release.sh v0.4.0                  # both architectures, pushed to Docker Hub
./release.sh v0.4.0 --no-push        # amd64 only, loaded into local Docker
```

Each pull request adds one line under Unreleased in `CHANGELOG.md`.
