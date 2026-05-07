import os
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import parse  # noqa: E402


def _make_protein_level_data():
    runs = [f"run_{i:02d}" for i in range(1, 7)]
    proteins = ["BRD4_HUMAN", "ATX2L_HUMAN", "VATG1_HUMAN", "CHK1_HUMAN", "TP53_HUMAN"]
    rows = []
    rng = np.random.default_rng(0)
    for r_idx, run in enumerate(runs, start=1):
        bait = "BRD4" if r_idx <= 3 else "CHK1"
        # leave one (Protein, Run) cell missing per run to test sparse pivot
        skip_protein = proteins[(r_idx - 1) % len(proteins)]
        for p in proteins:
            if p == skip_protein:
                continue
            log_int = float(8.0 + rng.normal(0.0, 1.5))
            rows.append({
                "RUN": r_idx,
                "Protein": p,
                "LABEL": "L",
                "LogIntensities": log_int,
                "originalRUN": run,
                "GROUP": bait,
                "SUBJECT": f"{bait}_rep{(r_idx - 1) % 3 + 1}",
                "TotalGroupMeasurements": 100,
                "NumMeasuredFeature": 90,
                "MissingPercentage": 0.1,
                "more50missing": False,
                "NumImputedFeature": 5,
            })
    return pd.DataFrame(rows), runs


def _make_ed(runs):
    rows = []
    for run in runs:
        idx = int(run.split("_")[-1])
        is_test = idx <= 3
        rows.append({
            "Experiment Name": run,
            "Type": "T" if is_test else "C",
            "Bait": "BRD4" if is_test else "CHK1",
            "Replicate": ((idx - 1) % 3) + 1,
            "Bait ID": "BRD4_HUMAN" if is_test else "CHK1_HUMAN",
        })
    return pd.DataFrame(rows)


class TestParseMsstats(unittest.TestCase):
    def test_smoke(self):
        with tempfile.TemporaryDirectory() as tmp:
            pld, runs = _make_protein_level_data()
            ed = _make_ed(runs)

            pld_path = os.path.join(tmp, "ProteinLevelData.csv")
            ed_path = os.path.join(tmp, "ED.csv")
            out_dir = os.path.join(tmp, "out")
            pld.to_csv(pld_path, index=False)
            ed.to_csv(ed_path, index=False)

            n_exp, n_ctrl = parse.parse_msstats(pld_path, ed_path, out_dir)
            self.assertEqual(n_exp, 3)
            self.assertEqual(n_ctrl, 3)

            for fname in ["bait.txt", "prey.txt", "interaction.txt",
                          "to_CompPASS.csv", "msstats_qc.csv",
                          "ED.csv", "ProteinLevelData.csv", "proteinGroups.txt"]:
                self.assertTrue(os.path.exists(os.path.join(out_dir, fname)),
                                f"missing output: {fname}")

            interaction = pd.read_csv(
                os.path.join(out_dir, "interaction.txt"),
                sep="\t", header=None,
                names=["ExperimentName", "Bait", "Prey", "Intensity"],
            )

            self.assertGreater(len(interaction), 0)

            # Every ED experiment should appear; every prey writes one row per ED experiment
            n_proteins = len(pld["Protein"].unique())
            self.assertEqual(set(interaction["ExperimentName"]), set(runs))
            # write_interaction_file emits N_proteins * N_runs rows (zeros included)
            self.assertEqual(len(interaction), n_proteins * len(runs))

            nonzero = interaction[interaction["Intensity"] > 0]
            self.assertGreater(len(nonzero), 0)
            # Back-transform sanity: 2 ** N(8, 1.5) ≈ [1, 1e4]
            self.assertTrue((nonzero["Intensity"] >= 1.0).all())
            self.assertTrue((nonzero["Intensity"] <= 1e6).all())

            # Bait file matches ED rows
            bait = pd.read_csv(os.path.join(out_dir, "bait.txt"),
                               sep="\t", header=None,
                               names=["ExpName", "Bait", "Type"])
            self.assertEqual(len(bait), 6)
            self.assertEqual(set(bait["ExpName"]), set(runs))

            # Prey file: one line per protein, no header
            prey = pd.read_csv(os.path.join(out_dir, "prey.txt"),
                               sep="\t", header=None)
            self.assertEqual(len(prey), n_proteins)


if __name__ == "__main__":
    unittest.main()
