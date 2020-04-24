

def get_species_annotation_dic(csv_file):
    '''
    Takes as input the binary information about species annotation from Species_Annotation_fake.csv (OrtAn results) and
    transforms it into a dictionary where to each species is assigned a list of the KOs present in the species genome
    :param csv_file: str - path to Species_Annotation_fake.csv
    :return: dict: {specie1: [KO1, KO3, ...], specie2: [...], ...}
    '''

    with open(csv_file) as f:
        lines = f.readlines()

    species_names_list = lines[0].split(';')[1:]
    species_names_list = [x.rstrip() for x in species_names_list]

    annotation_dic = {}

    # Create dic entry for each specie, starting with an empty list
    for specie in species_names_list:
        annotation_dic[specie] = []

    species_number = len(species_names_list)

    # Go through all lines of information in the csv file, each lines corresponds to one KO
    for line in lines[1:]:
        line = line.split(';')
        ko_name = line[0]
        specie = 0  # indicates the index of the species in the species_names_list
        specie_column = 1  # indicates the index of the column containing information for the current species

        # Go through all columns, each column represents a species
        while specie < species_number:
            # When we found the number 1 in the file cell, we must add the KO to that species
            if line[specie_column].strip() == str(1):
                annotation_dic[species_names_list[specie]].append(ko_name)
            specie += 1  # to go for the next column/species
            specie_column += 1
    return annotation_dic
