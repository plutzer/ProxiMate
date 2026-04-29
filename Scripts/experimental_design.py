import csv
import re
from ed_exceptions import (
    EDInvalidTypeError,
    EDMissingColumnError,
    EDFileNotFoundError,
    EDFileEmptyError,
    EDInvalidReplicateError,
    EDInvalidGroupError,
)
from log_config import get_logger

logger = get_logger(__name__)


GROUP_WILDCARD = "*"
_GROUP_INT_PATTERN = re.compile(r"^[1-9][0-9]*$")


def _parse_group_cell(raw, row_context=None):
    """
    Normalize a raw Group cell into a parsed spec.

    Returns one of:
      - None: row has no group declared (empty/NaN/whitespace)
      - "*": universal (control rows only; caller enforces)
      - frozenset({int, ...}): one or more explicit group numbers

    Raises EDInvalidGroupError on malformed content. `row_context` (a 1-based
    row number for user messages) is attached to the exception when provided.
    """
    if raw is None:
        return None
    s = str(raw).strip()
    if s == "":
        return None
    if s == GROUP_WILDCARD:
        return GROUP_WILDCARD

    tokens = [t.strip() for t in s.split(",")]
    has_wildcard = GROUP_WILDCARD in tokens
    explicit = [t for t in tokens if t != GROUP_WILDCARD and t != ""]

    if has_wildcard and explicit:
        raise EDInvalidGroupError(
            "control_wildcard_with_explicit_groups",
            [row_context] if row_context is not None else [],
        )

    values = set()
    for t in explicit:
        if not _GROUP_INT_PATTERN.match(t):
            raise EDInvalidGroupError(
                "invalid_group_value",
                [row_context] if row_context is not None else [],
            )
        values.add(int(t))

    if not values:
        # Only whitespace/commas after stripping
        raise EDInvalidGroupError(
            "invalid_group_value",
            [row_context] if row_context is not None else [],
        )

    return frozenset(values)


class Experiment:
    def __init__(self, header, values, row_number=None):
        self.attributes = dict(zip(header, values))
        if self.attributes["Type"] not in ["C", "T"]:
            raise EDInvalidTypeError([f"Type='{self.attributes['Type']}'"])

        if "Group" in self.attributes:
            self.group_spec = _parse_group_cell(
                self.attributes["Group"], row_context=row_number
            )
        else:
            self.group_spec = None


class ExperimentalDesign:
    def __init__(self, filepath):
        logger.info("Parsing experimental design file: %s", filepath)

        self.name2experiment = dict()

        try:
            with open(filepath, encoding='utf-8-sig') as csvfile:
                reader = csv.reader(csvfile)
                header = next(reader)

                # check that header has the minimum columns
                self.__check_column(header, "Experiment Name")
                self.__check_column(header, "Type")
                self.__check_column(header, "Bait")
                self.__check_column(header, "Replicate")

                for row_idx, row in enumerate(reader, start=2):
                    self.name2experiment[row[header.index("Experiment Name")]] = Experiment(
                        header, row, row_number=row_idx
                    )

        except FileNotFoundError:
            raise EDFileNotFoundError(filepath)
        except StopIteration:
            # Empty file or no rows after header
            raise EDFileEmptyError(filepath)

        self.__print_stats()

    def __check_column(self, header, column):
        if column not in header:
            raise EDMissingColumnError([column])

    def __print_stats(self):
        num_experiments = 0
        num_controls = 0
        baits = set()
        max_replicates = 0

        for experiment in self.name2experiment.values():
            if experiment.attributes["Type"] == "C":
                num_controls += 1
            else:
                num_experiments += 1
                baits.add(experiment.attributes["Bait"])
                try:
                    max_replicates = max(max_replicates, int(experiment.attributes["Replicate"]))
                except ValueError:
                    # Handle non-numeric replicate
                    raise EDInvalidReplicateError([experiment.attributes["Experiment Name"]])

        self.num_experiments = num_experiments
        self.num_controls = num_controls

    def is_grouped(self):
        """True iff any experiment has a non-None parsed group_spec."""
        return any(e.group_spec is not None for e in self.name2experiment.values())

    def get_groups(self):
        """Sorted unique group numbers drawn from test rows only."""
        groups = set()
        for e in self.name2experiment.values():
            if e.attributes["Type"] == "T" and isinstance(e.group_spec, frozenset):
                groups |= e.group_spec
        return sorted(groups)

    def get_experiments_for_group(self, group):
        """Return (test_experiments, control_experiments) that belong to `group`.

        Includes universal controls ('*'), explicit-match controls, and test
        rows whose group_spec contains `group`. Test rows with no group_spec
        are excluded (they are errors in grouped runs, handled by validation).
        """
        tests = []
        controls = []
        for e in self.name2experiment.values():
            spec = e.group_spec
            if e.attributes["Type"] == "T":
                if isinstance(spec, frozenset) and group in spec:
                    tests.append(e)
            else:  # C
                if spec == GROUP_WILDCARD:
                    controls.append(e)
                elif isinstance(spec, frozenset) and group in spec:
                    controls.append(e)
        return tests, controls
