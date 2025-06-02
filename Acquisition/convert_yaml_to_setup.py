# Example usage
#convert_yaml_to_setup('setup.yaml')

import yaml
import os
import shutil
import argparse

def convert_yaml_to_setup(yaml_file_path):

    with open(yaml_file_path, 'r') as yaml_file:
        data = yaml.safe_load(yaml_file)

    # Define the processed data directory for saving the YAML file
    processed_data_path = '/Volumes/EPSI_PROCESSING/Current_Cruise/Processed/'
    survey_name = data.get('survey_name')
    yaml_output_directory = os.path.join(processed_data_path, survey_name)

    # Create the directory for YAML output if it does not exist
    os.makedirs(yaml_output_directory, exist_ok=True)
    
    # Save the original YAML file in the specified directory
    yaml_output_path = os.path.join(yaml_output_directory, f"{survey_name}.yml")
    shutil.copy(yaml_file_path, yaml_output_path)

    # Define the path for the hardcoded Setup file location
    setup_file_hardcoded_path = '/Users/Shared/Software_current_cruise/MOD_fish_lib/Acquisition/fctd_epsi_acq/build/fctd_epsi/Build/Products/Debug/Setup'

    with open(setup_file_hardcoded_path, 'w') as setup_file:
        setup_file.write("%TCPIP port for send and receive data\n")
        # TCPIP and UDP setup based on 'in_use' flag
        ipad_tcip_in_use = data.get('comms', {}).get('ipad', {}).get('tcip_in_use', 0)
        ipad_udp_in_use = data.get('comms', {}).get('ipad', {}).get('udp_in_use', 0)

        # TCIP line
        if ipad_tcip_in_use == 1: 
            setup_file.write(f"TCPIPSocket.portnum={data['comms']['ipad']['tcip']}\n")
        else:
            setup_file.write(f"%TCPIPSocket.portnum={data['comms']['ipad']['tcip']}\n")

        # UDP line    
        if ipad_udp_in_use == 1: 
            setup_file.write(f"UDPSocket.portnum={data['comms']['ipad']['udp']}\n")
        else:
            setup_file.write(f"%UDPSocket.portnum={data['comms']['ipad']['udp']}\n") 

        # SENSOR SETUP
        setup_file.write("%SENSOR SETUP\n")
        setup_file.write(f"% survey name\nCTD.survey='{data['survey_name']}'\n")
        setup_file.write(f"% experiment name\nCTD.experiment='{data['experiment_name']}'\n")
        setup_file.write(f"% cruise name\nCTD.cruise='{data['mission_name']}'\n")
        setup_file.write(f"% vehicle name\nCTD.vehicle='{data['vehicle_name']}'\n")
        setup_file.write(f"% vehicle pressure case\nCTD.fish_pc='{data['pressure_case']}'\n")
        setup_file.write(f"% CTD.fishflag\nCTD.fishflag='{data['fish_flag']}'\n")

        # CTD serial number
        sn = data.get('sn', {})
        setup_file.write(f"% CTD SerialNum\nCTD.SerialNum='{sn.get('ctd', '0000')}'\n")

        # Calibration paths
        cal_paths = data.get('paths_from_acquisition_machine',{}).get('calibrations',{})
        setup_file.write(f"CTD.shear_Probecal_path='{cal_paths.get('shear','')}'\n")
        setup_file.write(f"CTD.FPO7_Probecal_path='{cal_paths.get('fpo7','')}'\n")
        setup_file.write(f"CTD.CTD_cal_path='{cal_paths.get('ctd','')}'\n")

        # Probe serial numbers
        setup_file.write(f"CTD.ch1_sn='{sn.get('t1', '')}'\n")
        setup_file.write(f"CTD.ch2_sn='{sn.get('t2', '')}'\n")
        setup_file.write(f"CTD.ch3_sn='{sn.get('s1', '')}'\n")
        setup_file.write(f"CTD.ch4_sn='{sn.get('s2', '')}'\n")
        setup_file.write(f"CTD.MicroCond='{sn.get('microcond','')}'\n")
        setup_file.write(f"CTD.Compass='{sn.get('compass','')}'\n")

        # Display
        display = data.get('display',{})
        setup_file.write(f"CTD.printData={display['print_data']}\n")
        setup_file.write(f"CTD.engDispRate = {display['engineer_display_rate']}\n")

        # Communication ports and speeds
        comms = data.get('comms', {})
        setup_file.write(f"%Serial port for Fish and PCode\nCTD.CTDPortName='{comms['data']['port']}'\n")
        setup_file.write(f"CTD.CommandPortName='{comms['commands']['port']}'\n")
        setup_file.write(f"CTD.speed={comms['commands']['speed']}\n")
        
        # Length and format of CTD data
        ctd = data.get('ctd',{})
        setup_file.write(f"CTD.CTDlength={ctd['length']}\n")
        setup_file.write(f"CTD.output_format={ctd['output_format']}\n")

        # GPS settings based on 'in_use' flag
        gps_in_use = comms.get('gps', {}).get('in_use', 0)
        if gps_in_use == 1:
            setup_file.write(f"PCode.PCodePortName='{comms['gps']['port']}'\n")
            setup_file.write(f"PCode.speed={comms['gps']['speed']}\n")
        else:
            setup_file.write(f"%PCode.PCodePortName='{comms['gps']['port']}'\n")
            setup_file.write(f"%PCode.speed={comms['gps']['speed']}\n")

        # ASCII file settings
        paths = data.get('paths_from_acquisition_machine', {})
        raw_files = data.get('raw_files', {})

        setup_file.write(f"%File settings\n")
        setup_file.write(f"Ascii_dataFile.runname='{raw_files.get('name_prefix', '')}'\n")
        setup_file.write(f"Ascii_dataFile.path='{paths.get('raw_incoming', '')}'\n")
        setup_file.write(f"Ascii_dataFile.dataFileSize={raw_files['filesize']}\n")

        # AFE settings
        afe = data.get('afe', {})
        setup_file.write(f"%AFE configuration\nAFE.channels = {list(afe.get('channels', {}).keys())}\n")

    print(f"YAML file has been saved at: {yaml_output_path}")
    print(f"Setup file has been created at: {setup_file_hardcoded_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert YAML to Setup format")
    parser.add_argument("yaml_file", type=str, help="Path to the YAML file to convert")
    args = parser.parse_args()

    # Run the conversion function with the provided YAML file
    convert_yaml_to_setup(args.yaml_file)

