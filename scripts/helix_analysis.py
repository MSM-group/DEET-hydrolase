# Import PyMOL library
import pymol
import os

# Initialize PyMOL
pymol.finish_launching(['pymol', '-qxi'])

# Define the folder path containing the protein structure files
folder_path = '//eawag/userdata/felderfl/Desktop/alpha_helix_test'

# Define the name of the output text file
output_file_name = 'alpha_helix_counts.txt'

# Get the list of protein structure files in the folder
protein_files = [f for f in os.listdir(folder_path) if f.endswith('.pdb')]

# Open the text file for appending
with open(output_file_name, 'a') as output_file:
    # Loop through each protein file
    for protein_file in protein_files:
        # Clear the previous protein if any
        cmd.reinitialize()
        
        # Load the specified protein file
        cmd.load(os.path.join(folder_path, protein_file), 'protein')

        # Select alpha helices
        selection_name = 'helix'
        cmd.select(selection_name, '(protein and ss h and name CA)')

        # Get the count of selected residues
        selected_residues_count = cmd.count_atoms(selection_name)

        output_message = f'{protein_file}, {selected_residues_count}\n'

        # Print the output to the console
        print(output_message)

        # Write the output to the text file
        output_file.write(output_message)

# Close the text file
output_file.close()

