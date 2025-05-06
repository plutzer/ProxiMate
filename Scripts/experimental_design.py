import csv


class Experiment:
    def __init__(self, header, values):
        self.attributes = dict(zip(header, values))
        if self.attributes["Type"] not in ["C", "T"]:
            print("Error: invalid experiment type: '" + self.attributes["Type"] + "' Expected: C or T")
            exit(1)


class ExperimentalDesign:
    def __init__(self, filepath):
        print("\nParsing experimental design file: " + filepath)

        self.name2experiment = dict()

        with open(filepath, encoding='utf-8-sig') as csvfile:
            reader = csv.reader(csvfile)
            header = next(reader)
            print("Header: ", header)

            # check that header has the minimum columns
            self.__check_column(header, "Experiment Name")
            self.__check_column(header, "Type")
            self.__check_column(header, "Bait")
            self.__check_column(header, "Replicate")

            for row in reader:
                self.name2experiment[row[header.index("Experiment Name")]] = Experiment(header, row)

        self.__print_stats()

    def __check_column(self, header, column):
        if column not in header:
            print("ERROR: missing required column: " + column)
            exit(1)

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
                max_replicates = max(max_replicates, int(experiment.attributes["Replicate"]))

        print("Num baits:", len(baits))
        print("Num bait experiments:", num_experiments)
        print("Max replicates per bait:", max_replicates)
        print("Num control experiments:", num_controls)

        self.num_experiments = num_experiments
        self.num_controls = num_controls

    # def get_num_baits(self):
    #     return self.num_baits
    
    # def get_num_ctrls(self):
    #     return len([experiment for experiment in self.name2experiment.values() if experiment.attributes["Type"] == "C"])
