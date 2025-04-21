import re
import os
import numpy as np

def load_inverse_xsec_data(directory):
    """
    Load inverse cross-section data from files into a dictionary.
    """
    file_map = {
        "NuMu": "numu_AR23_0p1-20GeV_inverse_xsec_as_flux.data",
        "NuMubar": "numubar_AR23_0p1-20GeV_inverse_xsec_as_flux.data",
        "NuE": "nue_AR23_0p1-20GeV_inverse_xsec_as_flux.data",
        "NuEbar": "nuebar_AR23_0p1-20GeV_inverse_xsec_as_flux.data",
        "NuTau": "nutau_AR23_0p1-20GeV_inverse_xsec_as_flux.data",
        "NuTaubar": "nutaubar_AR23_0p1-20GeV_inverse_xsec_as_flux.data",
    }
    
    inverse_xsec_data = {}
    for key, filename in file_map.items():
        filepath = os.path.join(directory, filename)
        if os.path.exists(filepath):
            data = np.loadtxt(filepath)
            inverse_xsec_data[key] = {row[0]: row[1] for row in data}
    
    return inverse_xsec_data

def find_closest_energy(energy, energy_values):
    """
    Find the closest energy value in the given list.
    """
    return min(energy_values, key=lambda x: abs(x - energy))

def modify_flux_values(values, energy, inverse_xsec_data):
    """
    Replace flux values using the closest matching inverse cross-section values.
    """
    column_keys = ["NuMu", "NuMubar", "NuE", "NuEbar", "NuTau", "NuTaubar"]
    modified_values = values.copy()
    
    for i, key in enumerate(column_keys, start=1):  # Start at 1 to skip energy column
        if key in inverse_xsec_data:
            closest_energy = find_closest_energy(energy, inverse_xsec_data[key].keys())
            modified_values[i] = str(inverse_xsec_data[key][closest_energy])
    
    return modified_values

def process_flux_file(input_file, output_file, energy_threshold=20.0, data_dir="."):
    """
    Reads the atmospheric neutrino flux file, removes energy bins above the threshold,
    and writes the modified data to a new file.
    
    Parameters:
    - input_file: str, path to the input file.
    - output_file: str, path to the output file.
    - energy_threshold: float, maximum allowed neutrino energy (in GeV).
    - data_dir: str, directory containing inverse xsec files.
    """
    inverse_xsec_data = load_inverse_xsec_data(data_dir)
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        keep_writing = True
        for line in infile:
            if re.match(r'average flux in \[cosZ =', line):
                keep_writing = True  # Keep the angular bin header
                outfile.write(line)
                continue
            
            if re.match(r'Enu\(GeV\)', line):
                outfile.write(line)  # Keep the column headers
                continue
            
            try:
                parts = line.split()
                energy = float(parts[0])
                if energy > energy_threshold:
                    keep_writing = False  # Stop writing data rows for this bin
                if keep_writing:
                    modified_values = modify_flux_values(parts, energy, inverse_xsec_data)
                    outfile.write(" ".join(modified_values) + "\n")
            except (ValueError, IndexError):
                outfile.write(line)  # Keep any non-numeric lines (e.g., empty lines)

if __name__ == "__main__":
    input_filename = "hms-ally-20-12-solmax_3FlavOsc.d"  # Change to the actual input file name
    output_filename = "filtered-hms-ally-20-12_3Flav_0p1-20GeV.d"  # Output file
    process_flux_file(input_filename, output_filename, data_dir=".")
