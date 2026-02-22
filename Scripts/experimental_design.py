import csv
from ed_exceptions import (
    EDInvalidTypeError,
    EDMissingColumnError,
    EDFileNotFoundError,
    EDFileEmptyError,
    EDInvalidReplicateError
)


class Experiment:
    def __init__(self, header, values):
        self.attributes = dict(zip(header, values))
        if self.attributes["Type"] not in ["C", "T"]:
            raise EDInvalidTypeError([f"Type='{self.attributes['Type']}'"])


class ExperimentalDesign:
    def __init__(self, filepath):
        print("\nParsing experimental design file: " + filepath)

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

                for row in reader:
                    self.name2experiment[row[header.index("Experiment Name")]] = Experiment(header, row)

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

    # def get_num_baits(self):
    #     return self.num_baits
    
    # def get_num_ctrls(self):
    #     return len([experiment for experiment in self.name2experiment.values() if experiment.attributes["Type"] == "C"])
