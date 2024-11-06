#!/bin/bash

while true
do
date +"copy RAW to DEV3 +%c"
rsync -auq mod_admin@192.168.1.168:/Users/Shared/FCTD_EPSI_DATA/Current_Cruise/ /Users/Shared/EPSI_PROCESSING/Current_Cruise/Realtime_RAW/raw
rsync -auq mod_admin@192.168.1.168:/Users/Shared/FCTD_EPSI_DATA/Current_Cruise/ /Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/RAW_full_cruise
rsync -auq mod_admin@192.168.1.168:/Users/Shared/FCTD_EPSI_DATA/Simulated_Data/Incoming_Data/ /Users/Shared/EPSI_PROCESSING/Simulated_Data/Realtime_RAW/raw
rsync -auq mod_admin@192.168.1.168:/Users/Shared/FCTD_EPSI_DATA/Simulated_Data/Incoming_Data/ /Users/Shared/EPSI_PROCESSING/Simulated_Data/Processed/RAW_full_cruise
sleep 1 
done


