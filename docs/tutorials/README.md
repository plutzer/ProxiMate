# ProxiMate tutorials

Task-focused walkthroughs of the ProxiMate GUI. Tutorials 1 and 3 to 5 use the Human
Cell Map datasets that ship with the Docker image; tutorial 2 uses the small example
dataset downloadable from inside the GUI. Start the GUI with
`docker run -p 3838:3838 plutzer/proximate` and open http://localhost:3838.

1. [Getting started](1-getting-started.md): the scores in plain language and the sidebar
   (active dataset, thresholds, presets).
2. [Parse and score](2-parse-and-score.md): download the example dataset, parse and
   score it, then duplicate it with one control left out and score the copy.
3. [Quality control](3-quality-control.md): PCA, threshold metrics, the SAINT scatter and
   the passing interactions for one bait.
4. [Exploring interactors](4-exploring-interactors.md): GO and UniProt feature
   enrichment, and comparing two baits.
5. [Networks and export](5-networks-and-export.md): draw, restyle, select and cluster a
   network in Cytoscape, then export results.

For a reference description of every control, see [`GUI/gui_docs.md`](../../GUI/gui_docs.md).
