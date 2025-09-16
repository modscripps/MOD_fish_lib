import os
import shutil
import sys

def create_deployment(deployment_name, fish_number):
    raw_dir_path = '/Volumes/Epsidrive/COMCEPTS/DATA/'  
    py_script_path = '/Volumes/Epsidrive/COMCEPTS/python_scripts/'  
    yml_script_path = '/Volumes/Epsidrive/COMCEPTS/data_processing_scripts/' 

    # Create a new directory with the provided deployment name
    deployment_path = os.path.join(raw_dir_path, deployment_name)
    os.makedirs(deployment_path, exist_ok=True)

    # Determine the file to copy based on fish_number
    if fish_number == 1:
        minnow_file = 'start_minnow1.py'
        yml_file = 'comcepts_minnow1_template.yml'
    elif fish_number == 3:
        minnow_file = 'start_minnow3.py'
        yml_file = 'comcepts_minnow3_template.yml'
    else:
        print("Invalid fish_number. Please provide 1 or 3.")
        return

    # Construct the full path of the minnow file and copy start_minnow script to the new directory
    minnow_file_path = os.path.join(py_script_path, minnow_file)
    shutil.copy(minnow_file_path, deployment_path)

    # Also copy stop_minnow script to the new directory
    minnow_file_path = os.path.join(py_script_path, 'stop_minnow.py')
    shutil.copy(minnow_file_path, deployment_path)

    # Construct the full path of the yml file
    yml_file_path = os.path.join(yml_script_path, yml_file)

    # Define the target path with the new name right away
    new_yml_path = os.path.join(deployment_path, f"{deployment_name}.yml")

    # Copy and rename in one step
    shutil.copy(yml_file_path, new_yml_path)

    # Modify the YAML file to include the deployment name after "deployment_name:"
    with open(new_yml_path,'r') as f:
        lines = f.readlines()

    with open(new_yml_path,'w') as f:
        for line in lines:
            if line.strip().startswith("deployment_name:"):
                f.write(f"deployment_name: {deployment_name}\n")
            else:
                f.write(line)


    # Make a subdirectory called 'raw' where you will store the raw data files
    os.makedirs(os.path.join(deployment_path, "raw"), exist_ok=True)

    print(f"Deployment '{deployment_path}' created successfully.")

    return deployment_path  # Return the path for later use


if __name__ == "__main__":
    deployment_name = input("Enter deployment name: ")
    fish_number = int(input("Enter fish number (1 or 3): "))

    created_deployment_path = create_deployment(deployment_name, fish_number)
    print(created_deployment_path)
