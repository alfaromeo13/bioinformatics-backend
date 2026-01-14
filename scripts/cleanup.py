import os
import zipfile

class Cleanup:
    def __init__(self, pdb_name):
        self.working_directory = os.getcwd()
        # Handle names safely even if not ending with ".pdb"
        self.base = os.path.splitext(os.path.basename(pdb_name))[0]
        self.zip_name = f"{self.base}_results.zip"

    def zip_files(self):
        file_list = os.listdir(self.working_directory)
        files_to_zip = []
        for filename in file_list:
            path = os.path.join(self.working_directory, filename)
            if not os.path.isfile(path):
                continue  # skip directories
            if filename == self.zip_name:
                continue  # don't include the zip we are creating

            # (joined*.pdb) OR (norm* | abs* | inter_ener*)
            if ((filename.startswith('joined') and filename.endswith('.pdb'))
                or filename.startswith(('norm', 'abs', 'inter_ener'))):
                files_to_zip.append(filename)

        if not files_to_zip:
            print("No matching files to zip.")
            return

        with zipfile.ZipFile(self.zip_name, 'w', compression=zipfile.ZIP_DEFLATED) as zip_file:
            for filename in files_to_zip:
                zip_file.write(os.path.join(self.working_directory, filename), arcname=filename)

        print(f"Results zipped to file {self.zip_name}")

    def delete_files(self):
        file_list = os.listdir(self.working_directory)
        for filename in file_list:
            path = os.path.join(self.working_directory, filename)
            if not os.path.isfile(path):
                continue  # skip directories
            if filename == self.zip_name:
                continue  # don't delete the archive we just made

            # Match your original intent, with clearer grouping
            if (
                filename.startswith('joined')
                or (len(filename) == 5 and filename.endswith('.pdb'))  # e.g., 'a.pdb'
                or filename.startswith(('pro', 'process'))
                or filename.endswith((
                    '.inp', '.jpg', '.svg', '.dat', '.out', '.rtf', '.log', '.prm'
                ))
            ):
                try:
                    os.remove(path)
                except OSError as e:
                    print(f"Could not delete {filename}: {e}")
