# ProxiMate tutorials

Task-focused walkthroughs of the ProxiMate GUI, using the Human Cell Map datasets that
ship with the Docker image. Start the GUI with `docker run -p 3838:3838 plutzer/proximate`
and open http://localhost:3838.

1. [Getting started](1-getting-started.md): the scores in plain language and the sidebar
   (active dataset, thresholds, presets).
2. [Parse and score](2-parse-and-score.md): turn a MaxQuant search and an experimental
   design into a scored dataset.
3. [Quality control](3-quality-control.md): PCA, threshold metrics, the SAINT scatter and
   the passing interactions for one bait.
4. [Exploring interactors](4-exploring-interactors.md): GO and UniProt feature
   enrichment, and comparing two baits.
5. [Networks and export](5-networks-and-export.md): draw, restyle, select and cluster a
   network in Cytoscape, then export results.

For a reference description of every control, see [`GUI/gui_docs.md`](../../GUI/gui_docs.md).
