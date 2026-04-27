import os
import pymzml
import numpy as np
import tkinter as tk
from tkinter import filedialog, messagebox
import logging  # Added for logging
import urllib.request
from tqdm import tqdm
from lock_mass_correction import LockMassCorrection



# Configure logging
logging.basicConfig(
    filename='processing.log',  # Log file name
    level=logging.DEBUG,  # Log level
    format='%(asctime)s - %(levelname)s - %(message)s'  # Log format
)

# Function to read default values from a configuration file
def read_config(file_path):
    defaults = {}
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):  # Ignore comments and blank lines
                    key, value = line.split('=')
                    defaults[key] = value
    except Exception as e:
        messagebox.showerror("Error", f"Could not read config file: {e}")
    return defaults

# Function to process mzML files
def process_files(input_dir, output_dir, mz_range, base_intensity_threshold, minimum_intensity_mz_output, decimal_places, filter_type, sn_threshold, intensity_threshold, sub_scan_filter_intensity_threshold, lock_mass_correction):
                  
    logging.basicConfig(filename='process_files.log', level=logging.DEBUG, 
                        format='%(asctime)s - %(levelname)s - %(message)s')

    helper_folder = os.path.join(output_dir, "helper files")
    os.makedirs(helper_folder, exist_ok=True)

    input_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.mzML')]
    
    for mzml_file in tqdm(input_files, desc="Processing mzML files", unit=" files"):
        logging.debug(f"Processing mzML file: {mzml_file}")
        run = pymzml.run.Reader(mzml_file)
        file_name = os.path.basename(mzml_file).split(".")[0]

        all_mz = []
        all_intensity = []

        helper_file_path = os.path.join(helper_folder, f"{file_name}_helper.txt")
        included_scan_count = 0

        logging.debug(f"Creating helper file at: {helper_file_path}")

        with open(helper_file_path, "w") as helper_file:
            helper_file.write(f"Processing file: {file_name}\n")
            helper_file.write(f"Minimum Intensity for m/z: {minimum_intensity_mz_output}\n")
            helper_file.write(f"Decimal Places: {decimal_places}\n")
            helper_file.write(f"Base Peak Intensity Threshold: {base_intensity_threshold}\n")
            helper_file.write(f"m/z Range: {mz_range[0]} - {mz_range[1]}\n")
            helper_file.write(f"Filter Type: {filter_type}\n")
            
            if filter_type == 'intensity':
                helper_file.write(f"Intensity Threshold: {sub_scan_filter_intensity_threshold}\n")
            elif filter_type == 'sn':
                helper_file.write(f"Signal-to-Noise (S/N) Threshold: {sn_threshold}\n")

            helper_file.write(f"\nScan\tBase Peak m/z\tBase Peak Intensity\tIncluded/Excluded\n")
            
            for spectrum in tqdm(run, desc=f"Processing spectra in {file_name}", unit=" spectra", leave=False, mininterval=1):
                if spectrum.ms_level == 1:
                    mz_values = spectrum.mz
                    intensity_values = spectrum.i

                    if len(intensity_values) == 0 or np.all(intensity_values == 0):
                        logging.error(f"Empty intensity values in spectrum ID: {spectrum.ID}")
                        helper_file.write(f"{spectrum.ID}\t-\t-\tExcluded (Empty intensity values)\n")
                        continue

                    try:
                        base_peak_intensity = max(intensity_values)
                    except ValueError as e:
                        logging.error(f"Error processing spectrum ID: {spectrum.ID}: {str(e)}")
                        helper_file.write(f"{spectrum.ID}\t-\t-\tExcluded (Error processing)\n")
                        continue

                    base_peak_mz = mz_values[np.argmax(intensity_values)]

                    if base_peak_intensity < base_intensity_threshold:
                        logging.debug(f"Spectrum ID: {spectrum.ID} excluded due to low base peak intensity ({base_peak_intensity:.2f} < {base_intensity_threshold}).")
                        helper_file.write(f"{spectrum.ID}\t{base_peak_mz:.{decimal_places}f}\t{base_peak_intensity:.2f}\tExcluded (Low base peak intensity)\n")
                        continue

                    # Apply m/z range filtering and sub-scan filtering
                    filtered_mz = []
                    filtered_intensity = []

                    for mz, intensity in zip(mz_values, intensity_values):
                        if mz_range[0] <= mz <= mz_range[1]:
                            if filter_type == 'intensity' and intensity >= sub_scan_filter_intensity_threshold:
                                filtered_mz.append(mz)
                                filtered_intensity.append(intensity)
                            elif filter_type == 'sn' and intensity / np.median(intensity_values) >= sub_scan_filter_sn_threshold:
                                filtered_mz.append(mz)
                                filtered_intensity.append(intensity)


                    # Append to overall list if there are valid m/z values and intensities above threshold
                    if filtered_mz and any(i >= sub_scan_filter_intensity_threshold for i in filtered_intensity):
                        all_mz.append(filtered_mz)
                        all_intensity.append(filtered_intensity)
                        included_scan_count += 1
                        logging.debug(f"Spectrum ID: {spectrum.ID} has valid m/z values and intensities above threshold.")
                        helper_file.write(f"{spectrum.ID}\t{base_peak_mz:.{decimal_places}f}\t{base_peak_intensity:.2f}\tIncluded\n")
                    else:
                        logging.debug(f"Spectrum ID: {spectrum.ID} excluded due to no valid m/z values or intensities below threshold.")
                        helper_file.write(f"{spectrum.ID}\t{base_peak_mz:.{decimal_places}f}\t{base_peak_intensity:.2f}\tExcluded (No valid m/z values or intensities below threshold)\n")


            # Log the total number of included scans
            logging.info(f"Total included scans for file {file_name}: {included_scan_count}")
            helper_file.write(f"\nTotal included scans: {included_scan_count}\n")


        # If any valid scans are included, calculate sums and save results
        if all_mz and all_intensity:  
            all_mz = np.concatenate(all_mz)  
            all_intensity = np.concatenate(all_intensity)


        # Calculate the summed intensity for each unique m/z value
        unique_mz = np.unique(all_mz)
        summed_intensity = []

        for mz in unique_mz:
            # Get indices of the current m/z value
            indices = np.where(all_mz == mz)
            # Calculate the summed intensity for the current m/z
            total_intensity = np.sum(all_intensity[indices]) if indices[0].size > 0 else 0
            summed_intensity.append(total_intensity)


        # Save the summed peak list to a text file
        output_file_path = os.path.join(output_dir, f"{file_name}_summed_spectrum.txt")

        with open(output_file_path, "w") as output_file:
            # Write headers
            if include_headers_var.get():
                output_file.write("m/z\tIntensity\n")


            # Apply filtering and write data to file
            data_written = False
            for mz, intensity in zip(unique_mz, summed_intensity):
                if intensity >= minimum_intensity_mz_output:
                    output_file.write(f"{mz:.{decimal_places}f}\t{intensity:.0f}\n")
                    data_written = True


            # Write "No data" message if no data found
            if not data_written:
                output_file.write("No data. Consider relaxing filters.\n")
                logging.info(f"No data found for file {file_name}. Consider relaxing filters.")

    logging.info("All files have been processed successfully.")
            
             
# Function to open file dialog for selecting an input directory
def select_input_directory():
    input_dir = filedialog.askdirectory(title="Select Input Directory")
    logging.info(f"Selected input directory: {input_dir}")
    input_dir_entry.delete(0, tk.END)
    input_dir_entry.insert(0, input_dir)

# Function to open file dialog for selecting the output directory
def select_output_directory():
    output_dir = filedialog.askdirectory(title="Select Output Directory")
    logging.info(f"Selected output directory: {output_dir}")
    output_dir_entry.delete(0, tk.END)
    output_dir_entry.insert(0, output_dir)

def submit_form():
    input_dir = input_dir_entry.get()
    output_dir = output_dir_entry.get()
    mz_min = float(mz_min_entry.get())
    mz_max = float(mz_max_entry.get())
    decimal_places = int(decimal_places_entry.get())
    minimum_intensity_mz_output = float(min_intensity_entry.get())
    base_intensity_threshold = float(base_peak_intensity_entry.get())

    filter_type = filter_type_var.get()
    if filter_type == 'sn':
        threshold = float(sn_threshold_entry.get())
    elif filter_type == 'intensity':
        threshold = float(intensity_threshold_entry.get())

    # Lock Mass Correction values
    lock_mass_correction = lock_mass_correction_checkbox.get()
    lock_mass_mz = float(lock_mass_mz_entry.get())
    lock_mass_window = float(lock_mass_window_entry.get())
    lock_mass_min_intensity = float(lock_mass_min_intensity_entry.get())

    if not input_dir or not output_dir:
        messagebox.showerror("Error", "Please select input and output directories.")
        return

    # Create the mz_range tuple here
    mz_range = (mz_min, mz_max)

    # Add logging here to track the input values
    logging.info(f"Input Directory: {input_dir}")
    logging.info(f"Output Directory: {output_dir}")
    logging.info(f"m/z range: {mz_range}, Filter Type: {filter_type}, Base Peak Intensity Threshold: {base_intensity_threshold}, "
                 f"Minimum Intensity for m/z: {minimum_intensity_mz_output}, Decimal Places: {decimal_places}")
    logging.info(f"Lock Mass Correction: {lock_mass_correction}, Lock Mass m/z: {lock_mass_mz}, Lock Mass Window: {lock_mass_window}, "
                 f"Lock Mass Minimum Intensity: {lock_mass_min_intensity}")

    # Pass checkbox state and general processing values to process_files
    threshold = float(sn_threshold_entry.get()) if filter_type == 'sn' else float(intensity_threshold_entry.get())
    process_files(input_dir, output_dir, mz_range, base_intensity_threshold, minimum_intensity_mz_output, decimal_places, filter_type,
                  threshold, threshold, threshold, lock_mass_correction)

    if lock_mass_correction:
        from lock_mass_correction import perform_lock_mass_correction
        perform_lock_mass_correction(input_dir, output_dir, lock_mass_mz, lock_mass_window, lock_mass_min_intensity, include_headers_var.get())

    messagebox.showinfo("Processing Complete", "All files have been processed.")

# Function to show the About box
def show_about():
    logging.info("About dialog displayed.")
    messagebox.showinfo("About", "mzSummer v5.0\n\nCreated by Nicholas Michael\n\n""Copyright 2024. All rights reserved. MichaelSoft\n\n"
            "This application processes mzML files and sums their MS1 spectra.")

# Function to show the Help box
def show_help():
    logging.info("Help dialog displayed.")
    messagebox.showinfo("Help", 
        """mzSummer Data Processing Workflow


Step 1: Input Directory Selection
User selects input directory containing mzML files. Software validates mzML files.

Step 2: Configuration Parameter Setup
User configures parameters:
Base Peak Intensity Threshold
m/z Range (min, max)
Filter Type (Intensity or S/N)
Sub-Scan Filter Threshold
Minimum Intensity for m/z Output
Decimal Places for m/z values
Lock Mass Correction (enabled/disabled)
Lock Mass m/z
Lock Mass Window
Lock Mass Minimum Intensity
Software validates user-input parameters.

Step 3: mzML File Reading
Software reads each mzML file, extracting:
Scan information (scan ID, retention time)
Mass spectrometry data (m/z values, intensities)

Step 4: Sub-Scan Selection
Software applies Base Peak Intensity Threshold filter.

Step 5: Sub-Scan Filtering
Software filters m/z values:
m/z Range filter
Filter Type (Intensity or S/N) threshold

Step 6: Lock Mass Correction (optional)
If enabled, software proportionally applies Lock Mass Correction across specified m/z range.

The script uses the following criteria to select the lock mass m/z:

Intensity-based selection: It chooses the m/z with the highest intensity within the specified lock mass window (lock_mass_window parameter).
Window boundaries: The lock mass window is defined as lock_mass_mz ± (lock_mass_window / 1e6) * lock_mass_mz. This means the window is centered around lock_mass_mz with a width proportional to lock_mass_mz.
Maximum intensity peak: Within this window, the script selects the m/z corresponding to the maximum intensity peak.

Important Note
Lock Mass Correction only applies within the specified m/z range. Filtering beyond this range will result in non-corrected .LMC files.

Step 7: Summed Intensity Calculation
Software calculates summed intensity for unique m/z values.

Step 8: Data Output
Software generates:
Text file (summed intensity data)
Helper file (base peak intensity, scan inclusion status)
LMC-corrected files (if Lock Mass Correction enabled)

Step 9: Logging
Software logs processing details:
Input file names
Processing parameters
Scans included/excluded
Lock Mass Correction status (if applied)

""")

            
# Tkinter GUI setup
root = tk.Tk()
root.title("mzSummer")
root.configure(bg='#B2EBF2')

# Create a variable to store the selected filter type
filter_type_var = tk.StringVar(value='intensity')

# Read default values from config file
config_path = "config.txt"  # Adjust path if necessary
config = read_config(config_path)

# Lock Mass Correction variables
lock_mass_correction_var = False
lock_mass_mz = 609.28066
lock_mass_window = 20
lock_mass_min_intensity = 10000

# Read config.txt to update Lock Mass Correction variables
with open('config.txt', 'r') as f:
    lines = f.readlines()
for line in lines:
    if line.startswith('lock_mass_correction_var'):
        lock_mass_correction_var = line.split('=')[1].strip().lower() == 'true'
    elif line.startswith('lock_mass_mz'):
        lock_mass_mz = float(line.split('=')[1].strip())
    elif line.startswith('lock_mass_window'):
        lock_mass_window = float(line.split('=')[1].strip())
    elif line.startswith('lock_mass_min_intensity'):
        lock_mass_min_intensity = int(line.split('=')[1].strip())

# GUI layout for input and output directories
tk.Label(root, text="Input Directory:", bg='#B2EBF2').grid(row=0, column=0, padx=10, pady=5, sticky="e")
input_dir_entry = tk.Entry(root, width=50)
input_dir_entry.grid(row=0, column=1, padx=10, pady=5)
tk.Button(root, text="Browse", command=select_input_directory).grid(row=0, column=2, padx=10, pady=5)

tk.Label(root, text="Output Directory:", bg='#B2EBF2').grid(row=1, column=0, padx=10, pady=5, sticky="e")
output_dir_entry = tk.Entry(root, width=50)
output_dir_entry.grid(row=1, column=1, padx=10, pady=5)
tk.Button(root, text="Browse", command=select_output_directory).grid(row=1, column=2, padx=10, pady=5)

# Sub-Scan Selector
tk.Label(root, text="**Sub-Scan Selector**", bg='#B2EBF2', font=("Arial", 12, "bold")).grid(row=2, column=0, columnspan=2, padx=10, pady=5)

tk.Label(root, text="Base Peak Intensity Threshold:", bg='#B2EBF2').grid(row=3, column=0, padx=10, pady=5, sticky="e")
base_peak_intensity_entry = tk.Entry(root, width=10)
base_peak_intensity_entry.grid(row=3, column=1, padx=10, pady=5, sticky="w")
base_peak_intensity_entry.insert(0, config.get('base_peak_intensity', '10000'))

# Sub-Scan Filter
tk.Label(root, text="**Sub-Scan Filter**", bg='#B2EBF2', font=("Arial", 12, "bold")).grid(row=4, column=0, columnspan=2, padx=10, pady=5)

tk.Label(root, text="m/z Range to Process (Min, Max):", bg='#B2EBF2').grid(row=5, column=0, padx=10, pady=5, sticky="e")
mz_min_entry = tk.Entry(root, width=10)
mz_min_entry.grid(row=5, column=1, padx=10, pady=5, sticky="w")
mz_min_entry.insert(0, config.get('mz_min', '800'))
mz_max_entry = tk.Entry(root, width=10)
mz_max_entry.grid(row=5, column=1, padx=10, pady=5, sticky="e")
mz_max_entry.insert(0, config.get('mz_max', '4000'))

tk.Label(root, text="Filter by:", bg='#B2EBF2').grid(row=6, column=0, padx=10, pady=5, sticky="e")

# Use column 1 for radio buttons, add padx to create small gap
tk.Radiobutton(root, text="Sub-Scan Filter Intensity", variable=filter_type_var, value='intensity', bg='#B2EBF2').grid(row=6, column=1, sticky="w")
tk.Radiobutton(root, text="Sub-Scan Filter S/N Threshold", variable=filter_type_var, value='sn', bg='#B2EBF2').grid(row=6, column=1, sticky="w", padx=(160, 0))

tk.Label(root, text="Sub-Scan Filter Intensity:", bg='#B2EBF2').grid(row=7, column=0, padx=10, pady=5, sticky="e")
intensity_threshold_entry = tk.Entry(root, width=10)
intensity_threshold_entry.grid(row=7, column=1, padx=10, pady=5, sticky="w")
intensity_threshold_entry.insert(0, config.get('intensity_threshold', '1000'))

tk.Label(root, text="Sub-Scan Filter S/N Threshold:", bg='#B2EBF2').grid(row=8, column=0, padx=10, pady=5, sticky="e")
sn_threshold_entry = tk.Entry(root, width=10)
sn_threshold_entry.grid(row=8, column=1, padx=10, pady=5, sticky="w")
sn_threshold_entry.insert(0, config.get('sn_threshold', '3'))

# Lock Mass Correction
tk.Label(root, text="**Lock Mass Correction**", bg='#B2EBF2', font=("Arial", 12, "bold")).grid(row=14, column=0, columnspan=2, padx=10, pady=5)

lock_mass_correction_checkbox = tk.BooleanVar(value=lock_mass_correction_var)
lock_mass_correction_checkbox_entry = tk.Checkbutton(root, text="Enable Lock Mass Correction", variable=lock_mass_correction_checkbox, bg='#B2EBF2')
lock_mass_correction_checkbox_entry.grid(row=15, column=0, columnspan=2, sticky="w")

tk.Label(root, text="Lock Mass m/z:", bg='#B2EBF2').grid(row=16, column=0, padx=10, pady=5, sticky="e")
lock_mass_mz_entry = tk.Entry(root, width=10)
lock_mass_mz_entry.grid(row=16, column=1, padx=10, pady=5, sticky="w")
lock_mass_mz_entry.insert(0, lock_mass_mz)

tk.Label(root, text="Lock Mass Window (ppm):", bg='#B2EBF2').grid(row=17, column=0, padx=10, pady=5, sticky="e")
lock_mass_window_entry = tk.Entry(root, width=10)
lock_mass_window_entry.grid(row=17, column=1, padx=10, pady=5, sticky="w")
lock_mass_window_entry.insert(0, lock_mass_window)

tk.Label(root, text="Lock Mass Minimum Intensity:", bg='#B2EBF2').grid(row=18, column=0, padx=10, pady=5, sticky="e")
lock_mass_min_intensity_entry = tk.Entry(root, width=10)
lock_mass_min_intensity_entry.grid(row=18, column=1, padx=10, pady=5, sticky="w")
lock_mass_min_intensity_entry.insert(0, lock_mass_min_intensity)

# Output Options
tk.Label(root, text="**Output Options**", bg='#B2EBF2', font=("Arial", 12, "bold")).grid(row=9, column=0, columnspan=2, padx=10, pady=5)

tk.Label(root, text="Minimum Intensity of m/z Output:", bg='#B2EBF2').grid(row=10, column=0, padx=10, pady=5, sticky="e")
min_intensity_entry = tk.Entry(root, width=10)
min_intensity_entry.grid(row=10, column=1, padx=10, pady=5, sticky="w")
min_intensity_entry.insert(0, config.get('min_intensity', '100'))

tk.Label(root, text="Decimal Places for m/z values:", bg='#B2EBF2').grid(row=11, column=0, padx=10, pady=5, sticky="e")
decimal_places_entry = tk.Entry(root, width=10)
decimal_places_entry.grid(row=11, column=1, padx=10, pady=5, sticky="w")
decimal_places_entry.insert(0, config.get('decimal_places', '4'))

include_headers_var = tk.BooleanVar(value=False)
include_headers_checkbox = tk.Checkbutton(root, text="Include Headers in Output", variable=include_headers_var, bg='#B2EBF2')
include_headers_checkbox.grid(row=12, column=0, columnspan=2, sticky="w")

# Submit button
tk.Button(root, text="Submit", command=submit_form).grid(row=12, column=0, columnspan=3, pady=10)

# About and Help buttons
tk.Button(root, text="About", command=show_about).grid(row=13, column=0, padx=10, pady=5)
tk.Button(root, text="Help", command=show_help).grid(row=13, column=1, padx=10, pady=5)

def update_filter_fields():
    filter_type = filter_type_var.get()
    if filter_type == 'sn':
        sn_threshold_entry.config(state='normal')
        intensity_threshold_entry.config(state='disabled')
    elif filter_type == 'intensity':
        sn_threshold_entry.config(state='disabled')
        intensity_threshold_entry.config(state='normal')
    else:
        print("Invalid filter type selected")

# Call the function initially
update_filter_fields()

# Update the fields when the radio button selection changes
filter_type_var.trace('w', lambda *args: update_filter_fields())

root.mainloop()