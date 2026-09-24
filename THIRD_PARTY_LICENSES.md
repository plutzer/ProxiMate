# Third-party software

ProxiMate itself is released under the MIT license in `LICENSE`. The repository also
carries the source of the tools below, which keep their own licenses and copyrights.

| Component | Location | License | Notes |
| --- | --- | --- | --- |
| SAINTexpress 3.6.3 | `saint/upstream/SAINT-MRF-int`, `saint/upstream/SAINT-MRF-spc`, `saint/upstream/va` | GPL-3.0-or-later, Copyright (C) 2011 Hyungwon Choi and Damian Fermin | Teo et al. 2014, J Proteomics 100:37-43, doi:10.1016/j.jprot.2013.10.023. `saint/patches/SAINT-MRF-int` holds ProxiMate's modified copies of seven of these files and is distributed under the same GPL terms. |
| Boost 1.57 | `saint/upstream/boostsrc` | Boost Software License 1.0 | Vendored with SAINTexpress; only `program_options` is built. |
| NLopt 2.3 | `saint/upstream/nloptsrc` | LGPL-2.1-or-later for the library; individual algorithms carry their own notices in `nloptsrc/COPYING` and the per-directory `COPYRIGHT` files | Vendored with SAINTexpress. |
| GOGO | `Scripts/GOGO` | No license file is distributed by the authors; used with attribution | Zhao and Wang 2018, Sci Rep 8:15107, doi:10.1038/s41598-018-33219-y. `apcluster.c` implements affinity propagation (Frey and Dueck 2007). |
| CORUM | `Datasets/corum_humanComplexes.txt` | CC BY 4.0 | Tsitsiridis et al. 2023, Nucleic Acids Res 51:D539-D545. |

The BioGRID, UniProt and Human Protein Atlas snapshots are downloaded at build time; their terms of use are at the respective download pages listed in `README.md`.
