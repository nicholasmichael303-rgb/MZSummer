import os
import pandas as pd
import scipy.signal as signal
import logging
import glob

class LockMassCorrection:
    def __init__(self, input_dir, output_dir, lock_mass_mz, lock_mass_window, lock_mass_min_intensity):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.lock_mass_mz = lock_mass_mz
        self.lock_mass_window = lock_mass_window
        self.lock_mass_min_intensity = lock_mass_min_intensity

    def apply(self, file=None, include_headers=True):
        if file:
            file_path = os.path.join(self.output_dir, file)
            print(f"Processing file: {file_path}")
            
            data = pd.read_csv(file_path, sep="\t", header=None)
            
            if data.empty or len(data.columns) < 2:
                print(f"Skipping file {file}: No data points.")
                return
            
            self.process_file(data, file, include_headers)
        
        else:
            for file in os.listdir(self.input_dir):
                if file.endswith(".txt"):  
                    file_path = os.path.join(self.input_dir, file)
                    print(f"Processing file: {file_path}")
                    
                    data = pd.read_csv(file_path, sep="\t", header=None)
                    
                    if data.empty or len(data.columns) < 2:
                        print(f"Skipping file {file}: No data points.")
                        continue
                    
                    self.process_file(data, file, include_headers)

    def process_file(self, data, file, include_headers):
        helper_folder = os.path.join(self.output_dir, "LockMassRegionHelperFiles")
        os.makedirs(helper_folder, exist_ok=True)
        
        corrected_folder = os.path.join(self.output_dir, "Lock-mass corrected")
        os.makedirs(corrected_folder, exist_ok=True)
        
        data.columns = ['m/z', 'Intensity']
        
        original_precision = data['m/z'].apply(lambda x: len(str(x).split('.')[1])).max()
        
        # Extract lock mass region
        lock_mass_region_file_path = os.path.join(helper_folder, f"{file}_LockMassRegion.txt")
        lock_mass_region_data = self.get_lock_mass_region(data)
        
        if lock_mass_region_data.empty:
            with open(lock_mass_region_file_path, 'w') as f:
                f.write("No lock mass detected. No correction applied.\n")
            print(f"Lock mass region saved: {lock_mass_region_file_path}")
            
            # Save original data as corrected (no correction applied)
            file_name, file_extension = os.path.splitext(file)
            corrected_file_path = os.path.join(corrected_folder, f"{file_name}_LMC{file_extension}")
            if include_headers:
                data.to_csv(corrected_file_path, sep="\t", index=False, float_format=f'%.{original_precision}f')
            else:
                data.iloc[1:].to_csv(corrected_file_path, sep="\t", index=False, float_format=f'%.{original_precision}f', header=None)
            print(f"Corrected file saved: {corrected_file_path}")
        else:
            lock_mass_region_data.to_csv(lock_mass_region_file_path, sep="\t", index=False, float_format=f'%.{original_precision}f')
            
            with open(lock_mass_region_file_path, 'r+') as f:
                content = f.read()
                f.seek(0)
                lock_mass_mz_used = self.determine_lock_mass_mz(lock_mass_region_data)
                if lock_mass_mz_used in lock_mass_region_data['m/z'].values:
                    f.write(f"Lock mass correction successfully applied using m/z {lock_mass_mz_used:.6f}.\n")
                    lock_mass_error = self.calculate_lock_mass_error(lock_mass_mz_used)
                    f.write(f"Lock Mass m/z error = {'+' if lock_mass_error >= 0 else '-'}{abs(lock_mass_error):.6f} ppm\n")
                else:
                    f.write(f"Lock mass correction attempted but NOT applied (m/z {lock_mass_mz_used:.6f} not found).\n")
                f.write(content)
            
            print(f"Lock mass region saved: {lock_mass_region_file_path}")
            
            if lock_mass_mz_used in lock_mass_region_data['m/z'].values:
                # Apply lock mass correction
                corrected_data = self.correct(data)
                
                # Save corrected data
                file_name, file_extension = os.path.splitext(file)
                corrected_file_path = os.path.join(corrected_folder, f"{file_name}_LMC{file_extension}")
                if include_headers:
                    corrected_data.to_csv(corrected_file_path, sep="\t", index=False, float_format=f'%.{original_precision}f')
                else:
                    corrected_data.iloc[1:].to_csv(corrected_file_path, sep="\t", index=False, float_format=f'%.{original_precision}f', header=None)
                print(f"Corrected file saved: {corrected_file_path}")
                
    def correct(self, data):
        lock_mass_mz_used = self.determine_lock_mass_mz(self.get_lock_mass_region(data))
        lock_mass_error = self.calculate_lock_mass_error(lock_mass_mz_used)
        correction_factor = 1 - (lock_mass_error / 1e6)
        data['m/z'] = data['m/z'].apply(lambda x: x * correction_factor)
        return data
        
    def get_lock_mass_region(self, data):
        print("Extracting lock mass region...")
        
        data = data.loc[(data['m/z'] >= self.lock_mass_mz - 5) & 
                        (data['m/z'] <= self.lock_mass_mz + 5) &
                        (data['Intensity'] >= self.lock_mass_min_intensity)]  # Apply intensity threshold
        
        print(f"Lock mass region: {self.lock_mass_mz - 5} to {self.lock_mass_mz + 5}")
        
        if data.empty:
            print("No data points in lock mass region")
        
        return data
        
    def determine_lock_mass_mz(self, data):
        filtered_data = data[(data['m/z'] >= self.lock_mass_mz - self.lock_mass_window / 1e6 * self.lock_mass_mz) &
                        (data['m/z'] <= self.lock_mass_mz + self.lock_mass_window / 1e6 * self.lock_mass_mz)]

        if not filtered_data.empty:
            peak = filtered_data.loc[filtered_data['Intensity'].idxmax()]
            return peak['m/z']
        else:
            return self.lock_mass_mz

    def calculate_lock_mass_error(self, lock_mass_mz_used):
        lock_mass_error = (lock_mass_mz_used - self.lock_mass_mz) / self.lock_mass_mz * 1e6
        
        return lock_mass_error

    def apply_correction(self, file_path, lock_mass_error):
        print(f"Lock Mass Error: {lock_mass_error}")  
        
        if lock_mass_error is None or pd.isna(lock_mass_error):
            print("Lock Mass Error is None or NaN")
            return  
        
        data = pd.read_csv(file_path, sep="\t")
        
        data.loc[:, 'm/z'] = data['m/z'] * (1 + lock_mass_error / 1e6)
        
        data.to_csv(file_path, sep="\t", index=False)


def perform_lock_mass_correction(input_dir, output_dir, lock_mass_mz, lock_mass_window, lock_mass_min_intensity, include_headers=False):
    import logging
    
    logging.basicConfig(level=logging.DEBUG, filename='lock_mass_correction.log', format='%(asctime)s %(levelname)s: %(message)s')
    
    correction = LockMassCorrection(input_dir, output_dir, lock_mass_mz, lock_mass_window, lock_mass_min_intensity)
    
    os.makedirs(output_dir, exist_ok=True)
    
    files = glob.glob(output_dir + "/*_summed_spectrum.txt")
    
    logging.info(f"Files found: {files}")
    
    print(f"Found {len(files)} files in output_dir")
    
    helper_folder = os.path.join(output_dir, "LockMassRegionHelperFiles")
    os.makedirs(helper_folder, exist_ok=True)
    
    corrected_folder = os.path.join(output_dir, "Lock-mass corrected")
    os.makedirs(corrected_folder, exist_ok=True)
    
    helper_file_path = os.path.join(output_dir, "lock_mass_correction_helper.txt")
    
    with open(helper_file_path, "w") as helper_file:
        helper_file.write(f"Lock Mass m/z: {lock_mass_mz}\n")
        helper_file.write(f"Lock Mass Window: {lock_mass_window}\n")
        helper_file.write(f"Lock Mass Minimum Intensity: {lock_mass_min_intensity}\n")
        helper_file.write(f"Files processed: {len(files)}\n")
        helper_file.write("File Names:\n")
        for file_name in files:
            helper_file.write(f"- {os.path.basename(file_name)}\n")
    
    for file_name in files:
        file_basename = os.path.basename(file_name)
        try:
            correction.apply(file=file_basename, include_headers=include_headers)
            logging.info(f"Processed {file_name} successfully")
            
        except Exception as e:
            logging.error(f"Error processing {file_name}: {str(e)}")
            print(f"Error processing {file_name}: {str(e)}")
            
    logging.info("Lock mass correction completed")
    
    print("Helper files and folders created successfully.")