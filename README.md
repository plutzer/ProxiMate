# ProxiMate

Score and explore proximity-labeling proteomics data (BioID, TurboID, APEX) in a
browser. ProxiMate parses search-engine output, scores it with SAINTexpress and
CompPASS, annotates the interactions, and plots the results.

## Quickstart

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/).
2. Start ProxiMate:

   ```
   docker run -p 3838:3838 \
     --mount type=bind,source=<output_dir>,target=/Outputs \
     plutzer/proximate
   ```

   `<output_dir>` is a folder on your computer. Results are saved there; without
   the mount they are lost when the container is removed.
3. Open http://localhost:3838.
4. **Parse.** On the *Data & Scoring* tab, click **Download example dataset** and
   unzip it. Name the dataset, upload `proteinGroups.txt` and `ED.csv`, and click
   **Parse Data**.
5. **Score.** Click **Score Data**. The example takes about 20 seconds.
6. **Analyze.** Set score thresholds in the sidebar, then open the other tabs:

   | Tab | Use it to |
   | --- | --- |
   | Quality Control | check replicates by PCA and list passing interactions |
   | Feature Enrichment | find GO terms and domains enriched per bait |
   | Compare Baits | contrast two baits' networks |
   | Cytoscape | draw the network in Cytoscape desktop |
   | Export | download tables for other tools |

For your own data, upload MaxQuant, DIA-NN, Pioneer, FragPipe, MSstats or SAINT
files with an experimental design CSV. **Download ED template** gives the layout.

To pin a version, use a tag such as `plutzer/proximate:v0.4.0`. Each
[release](https://github.com/plutzer/ProxiMate/releases) lists its changes and
database dates. Images run on Intel/AMD and Apple Silicon.

## Tutorials

Step-by-step walkthroughs with screenshots, using data that ships with ProxiMate.

| Tutorial | You will |
| --- | --- |
| [1. Getting started](docs/tutorials/1-getting-started.md) | learn what the scores mean and how the sidebar thresholds work |
| [2. Parse and score](docs/tutorials/2-parse-and-score.md) | parse and score the example, then rescore it without one control |
| [3. Quality control](docs/tutorials/3-quality-control.md) | read the PCA, the threshold metrics and the SAINT scatter |
| [4. Exploring interactors](docs/tutorials/4-exploring-interactors.md) | run feature enrichment and compare two baits |
| [5. Networks and export](docs/tutorials/5-networks-and-export.md) | draw, style and cluster a network in Cytoscape, then export |

## More documentation

| Topic | Page |
| --- | --- |
| Every control in the GUI | [GUI reference](GUI/gui_docs.md) |
| Running without the GUI, on a cluster | [Command line](docs/command-line.md) |
| Connecting to Cytoscape | [Cytoscape setup](docs/cytoscape.md) |
| Driving ProxiMate from Claude Code or Codex | [Agent access (MCP)](docs/mcp.md) |
| Log files, `run.json`, environment variables | [Logs and provenance](docs/logs-and-provenance.md) |
| Annotation databases, Human Cell Map exclusion | [Databases](docs/databases.md) |
| Building, testing, releasing | [Development](docs/development.md) |
| Sample input files | [examples/](examples/README.md) |
| What changed in each version | [CHANGELOG.md](CHANGELOG.md) |

## Citation

Cite the version you used with the DOI Zenodo mints for each release, or with
`CITATION.cff` (GitHub's "Cite this repository" button). Please also cite the
scoring methods: SAINTexpress (Teo et al. 2014, J Proteomics 100:37-43) and
CompPASS (Sowa et al. 2009, Cell 138:389-403).

## License

MIT (`LICENSE`). Bundled SAINTexpress (GPL-3), Boost, NLopt and GOGO keep their
own licenses; see `THIRD_PARTY_LICENSES.md`.
