

# check if possibility can be performed by this species
def check_possibility(specie_annotation, possibility):
    present = True
    if isinstance(possibility, list):
        for ko in possibility:
            if ko not in specie_annotation:
                present = False
    elif possibility not in specie_annotation:
        present = False

    return present


# check is a complete path is present or not
def check_complete_path(specie_annotation, path, gpr_rules):
    """
    Check if the species has the complete path
    :param specie_annotation: dict - species annotation {specie1: [KO1, KO3, ...], specie2: [...], ...}
    :param path: dict - of the pathway to check {Reaction1: [enzyme1, enzyme2], reaction2: [enzyme3]}
    :param gpr_rules: dict - {e1: [KO1, KO2], e2: [KO3, [KO4,KO5]]}
    :return: boolean - True if the path is complete, False if is not
    """

    # for each reaction in the pahtway
    for reaction in path:
        # assume that the reaction is not present
        reaction_present = False

        # for each enzyme related with the reaction
        for enzyme in path[reaction]:
            # assume that the enzyme is not present
            enzyme_present = False

            # for each possibility in the gpr rules for the enzyme in question
            for possibility in gpr_rules[enzyme]:
                # check if the possibility is guaranteed
                if check_possibility(specie_annotation, possibility):
                    # enzyme is present, break the loop
                    enzyme_present = True
                    break

            if enzyme_present:
                # if the enzyme is present, reaction is present, break the loop
                reaction_present = True
                break

        if not reaction_present:
            # if we found a reaction is not present we can return False right way instead of calculating all the others
            return False

    return True
