from itertools import combinations
from check_complete_path import check_complete_path


class MicrobialInteractions:
    def __init__(self, species_annotation_dict, gpr_rules):
        self.species_annotation_dict = species_annotation_dict  # a dictionary containing all the KO's present in each specie {specie1: [ko1, ko3..], ...}
        self.gpr_rules = gpr_rules  # dictionary containing GPR_rules {ec_number: [kos]}

    def get_species_groups_performing_path(self, path_dict, number_of_species_in_each_combination, species_to_exclude=[]):
        """

        :param path_dict: dict - dictionary containing information about the pathway of interest {Reaction1: [ec_number1, ec_number2], reaction2: [ec_number3]}
        :param number_of_species_in_each_combination: int - number of species in each group
        :param species_to_exclude: list - list of species to exclude from the species_annotation_dict (species that are able to perform the complete path alone)
        :return: list - list of groups of species able to perform the complete path
        """

        # get list of species to use
        species_to_use = set(self.species_annotation_dict.keys()).difference(species_to_exclude)

        # get all possible groups of species combinations
        species_combinations = get_species_combinations(species_to_use, number_of_species_in_each_combination)
        # create an annotation dictionary for each one of the groups, to treat each group as if it was a single species
        fake_species_annotation_dict = self.create_fake_species(species_combinations)


        species_groups_with_complete_path = []

        # for each on of the calculated groups of species
        for species_group in fake_species_annotation_dict:
            # check if the complete path is present
            path_present = check_complete_path(fake_species_annotation_dict[species_group], path_dict, self.gpr_rules)
            if path_present:
                # if path is present add group to the result list to return
                species_groups_with_complete_path.append(species_group)

        return species_groups_with_complete_path

    def create_fake_species(self, species_combinations):
        res = {}
        for species_group in species_combinations:
            kos_set = set()
            for specie in species_group:
                kos_set.update(self.species_annotation_dict[specie])
            res[species_group] = list(kos_set)
        return res


def get_species_combinations(species_list, number_of_species):
    species_combinations = list(combinations(species_list, number_of_species))
    return species_combinations





