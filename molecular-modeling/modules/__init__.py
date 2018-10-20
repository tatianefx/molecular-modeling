from modules.chemistry.amino_acid import AminoAcid

amino_acid = input("Digite uma sequência de aminoácidos: ")

if AminoAcid.is_valid(amino_acid):
    AminoAcid.generate_amino_acid_molecule(amino_acid)


