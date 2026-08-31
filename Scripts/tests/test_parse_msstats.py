import numpy as np
import pandas as pd
import pytest

import parse


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


@pytest.fixture
def msstats_outputs(tmp_path):
    """Run parse_msstats on synthetic inputs; yield its return value and output dir."""
    pld, runs = _make_protein_level_data()
    ed = _make_ed(runs)

    pld_path = tmp_path / "ProteinLevelData.csv"
    ed_path = tmp_path / "ED.csv"
    out_dir = tmp_path / "out"
    pld.to_csv(pld_path, index=False)
    ed.to_csv(ed_path, index=False)

    n_exp, n_ctrl = parse.parse_msstats(str(pld_path), str(ed_path), str(out_dir))
    return {
        "pld": pld,
        "runs": runs,
        "out_dir": out_dir,
        "n_exp": n_exp,
        "n_ctrl": n_ctrl,
        "n_proteins": len(pld["Protein"].unique()),
    }


def test_counts_test_and_control_experiments(msstats_outputs):
    assert msstats_outputs["n_exp"] == 3
    assert msstats_outputs["n_ctrl"] == 3


@pytest.mark.parametrize("fname", [
    "bait.txt", "prey.txt", "interaction.txt", "to_CompPASS.csv",
    "msstats_qc.csv", "ED.csv", "ProteinLevelData.csv", "proteinGroups.txt",
])
def test_output_file_written(msstats_outputs, fname):
    assert (msstats_outputs["out_dir"] / fname).exists(), f"missing output: {fname}"


def test_interaction_file_covers_every_protein_run_pair(msstats_outputs):
    interaction = pd.read_csv(
        msstats_outputs["out_dir"] / "interaction.txt",
        sep="\t", header=None,
        names=["ExperimentName", "Bait", "Prey", "Intensity"],
    )

    assert len(interaction) > 0
    # Every ED experiment should appear; every prey writes one row per ED experiment
    assert set(interaction["ExperimentName"]) == set(msstats_outputs["runs"])
    # write_interaction_file emits N_proteins * N_runs rows (zeros included)
    assert len(interaction) == msstats_outputs["n_proteins"] * len(msstats_outputs["runs"])


def test_interaction_intensities_are_back_transformed(msstats_outputs):
    interaction = pd.read_csv(
        msstats_outputs["out_dir"] / "interaction.txt",
        sep="\t", header=None,
        names=["ExperimentName", "Bait", "Prey", "Intensity"],
    )

    nonzero = interaction[interaction["Intensity"] > 0]
    assert len(nonzero) > 0
    # Back-transform sanity: 2 ** N(8, 1.5) ≈ [1, 1e4]
    assert (nonzero["Intensity"] >= 1.0).all()
    assert (nonzero["Intensity"] <= 1e6).all()


def test_bait_file_matches_ed_rows(msstats_outputs):
    bait = pd.read_csv(msstats_outputs["out_dir"] / "bait.txt",
                       sep="\t", header=None,
                       names=["ExpName", "Bait", "Type"])
    assert len(bait) == 6
    assert set(bait["ExpName"]) == set(msstats_outputs["runs"])


def test_prey_file_has_one_line_per_protein(msstats_outputs):
    prey = pd.read_csv(msstats_outputs["out_dir"] / "prey.txt",
                       sep="\t", header=None)
    assert len(prey) == msstats_outputs["n_proteins"]
