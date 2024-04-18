

def find_closest_res_to_reference(prot1,prot2):

    # Initialize a list to store the indices of the closest residues
    closest_residues = []
    
    # Define the residues of interest in protein1
    residues_of_interest = [213, 251]
    
    # For each residue of interest, find the closest residue in protein2
    for resi in residues_of_interest:
        # Create a selection for the C-alpha atom of the residue of interest
        cmd.select('resi_interest', f'{prot1} and resi {resi} and name CA')
    
        # Get the list of all residues in protein2
        residues_in_protein2 = cmd.get_model(f'{prot2} and name CA').atom
    
        # Initialize the minimum distance to a large number
        min_distance = 9999
        nearest_resi = None
        nearest_resn = None
    
        # Iterate over all residues in protein2
        for atom in residues_in_protein2:
            # Calculate the distance between the current atom and the residue of interest
            distance = cmd.get_distance('resi_interest', f'{prot2} and resi {atom.resi} and name CA')
    
            # If the calculated distance is less than the minimum distance, update the minimum distance and the nearest residue
            if distance < min_distance:
                min_distance = distance
                nearest_resi = atom.resi
                nearest_resn = atom.resn
    
        # Print the result
        print(f'The closest residue to residue {resi} in {prot1} is residue {nearest_resn} {nearest_resi} in {prot2}')
    
        # Add the index of the closest residue to the list
        closest_residues.append(int(nearest_resi))
        
        return(closest_residues)

def get_alpha_helix_string(molecule, closest_residues, num_zeros):
    # Assign secondary structure
    cmd.dss(molecule)

    # Initialize the result string
    result = ""

    # Define a function to be called for each residue
    def check_residue(residue_name, ss, resi):
        nonlocal result
        if ss == 'H':
            result += '1'
        else:
            result += '0'

        # If the current residue is one of the closest residues, change the character to '2'
        if int(resi) in closest_residues:
            result = result[:-1] + '2'

    # Iterate over each residue
    cmd.iterate(f'{molecule} and name ca', 'check_residue(resn, ss, resi)', space=locals())

    # Store the original result
    original_result = result

    # Treat '2' as '1'
    result = result.replace('2', '1')

    # Check for '1's and change '0's to '1's if there are fewer than num_zeros '0's before and after
    result_list = list(result)
    for i in range(len(result_list)):
        if result_list[i] == '1':
            # Find the start and end of the cluster of '1's
            start = i
            while start > 0 and result_list[start - 1] == '1':
                start -= 1
            end = i
            while end < len(result_list) - 1 and result_list[end + 1] == '1':
                end += 1

            # Check the number of '0's before the cluster
            num_zeros_before = 0
            for j in range(start - 1, -1, -1):
                if result_list[j] == '0':
                    num_zeros_before += 1
                else:
                    break

            # Check the number of '0's after the cluster
            num_zeros_after = 0
            for j in range(end + 1, len(result_list)):
                if result_list[j] == '0':
                    num_zeros_after += 1
                else:
                    break

            # Change '0's to '1's if there are fewer than num_zeros '0's before and after the cluster
            if num_zeros_before < num_zeros:
                for j in range(start - num_zeros_before, start):
                    result_list[j] = '1'
            if num_zeros_after < num_zeros:
                for j in range(end + 1, end + 1 + num_zeros_after):
                    result_list[j] = '1'

    result = ''.join(result_list)

    # Change '2's back to '2'
    for i in range(len(result)):
        if original_result[i] == '2':
            result = result[:i] + '2' + result[i+1:]

    # Get the residue sequence of prot2
    residue_sequence = cmd.get_fastastr(molecule).split('\n', 1)[1].replace('\n', '')

    # Find clusters of '1's that contain a '2' and count how many '1's and '2's they contain
    helix_lengths = []
    i = 0
    while i < len(result):
        if result[i] == '2':
            start = i
            while start > 0 and result[start - 1] == '1':
                start -= 1
            end = i
            while end < len(result) - 1 and result[end + 1] == '1':
                end += 1
            helix_lengths.append(end - start + 1)
            i = end + 1
        else:
            i += 1

    # Write the results to a text file
    with open('Result_Helixes.txt', 'a') as f:
        f.write(f'{molecule}\n')
        f.write(f'{residue_sequence}\n')
        f.write(f'{original_result}\n')
        f.write(f'{result}\n')
        f.write(f'{helix_lengths}\n')


import pymol
from pymol import cmd
import os





def process_folder(folder_path):
    # Change to the specified folder
    os.chdir(folder_path)

    # Get a list of all .pdb files in the folder
    pdb_files = [f for f in os.listdir() if f.endswith(".pdb")]

    for pdb_file in pdb_files:
        # Load the PDB file and assign the name 'prot1'
        cmd.load(pdb_file, 'prot2')
        cmd.load('yourfolderpath\Reference.pdb",'prot1')
        
        prot1='prot1'
        prot2='prot2'
        
        cmd.super(prot1, prot2)
        
        get_alpha_helix_string(prot2, find_closest_res_to_reference(prot1,prot2), 6)
        
        
        # Reinitialize Pymol

        cmd.reinitialize()


folder_path = 'yourfolderpath'

# Call the function to process the folder
process_folder(folder_path)

