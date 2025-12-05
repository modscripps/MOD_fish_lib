#!/bin/sh

echo "starting rsync to A-DRIVE at $(date)" >> /Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/backup_logs/mod_backup.log

rsync -auq --delete /Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/ /Volumes/A-DRIVE/MOTIVE2025/05_processed_data/fctd_epsi/
# rsync -auq --delete /Users/Shared/EPSI_PROCESSING/Processed/ /Volumes/B-DRIVE/TFO2023/05_processed_data/fctd_epsi

echo "rsync to A-DRIVE finished at $(date)" >> /Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/backup_logs/mod_backup.log

# Save this script on DEV3
#
# To make the file executable
# chmod a+x copy_raw_to_DEV3
#
# To add to the crontab
# crontab -e
#
# Type the following and save with :w
# 0 * * * * /Users/Shared/FCTD_EPSI/copy_raw_to_DEV3
# (The asterisks say when and how often to run the file. This means run every hour)
# crontab.guru - website to help you figure out how to schedule crontabs 
