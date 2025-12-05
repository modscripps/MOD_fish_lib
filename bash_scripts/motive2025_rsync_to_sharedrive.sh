#!/bin/sh

# write a line to log file
echo "starting rsync to science share at $(date)" >> /Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/backup_logs/mod_backup.log

# sync raw data
rsync -aq --delete /Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/Raw_full_cruise /Volumes/science/ScienceShare/MOD/03_raw_mod_data/

# sync processed data :: exclude raw and twist count files as they can always be regenerated
rsync -aq --delete --exclude 'RAW_full_cruise/' --exclude 'RAW_full_cruise_twist_counter/' --exclude 'MAT_full_cruise_twist_counter/' --exclude 'ROT_full_cruise_twist_counter/' /Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/ /Volumes/science/ScienceShare/MOD/05_processed_data/

# write a line to log file
echo "rsync to science share finished at $(date)" >> /Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/backup_logs/mod_backup.log

# Save this script on DEV3
#
# To make the file executable
# chmod a+x rsync_to_sharedrive_motive2025
#
# To add to the crontab
# crontab -e
#
# Type the following and save with :w
# 0 0/6 * * * /Volumes/DEV1_HD/Users/Shared/Software_current_cruise/MOD_fish_lib/bash_scripts/rsync_to_sharedrive_motive2025.sh

# Alternatively, use this line to have a lock file to prevent several backup processes at once:
# 2/5 * * * * /usr/bin/lockf -t 0 -k /tmp/rsync_sharedrive.lock /Volumes/DEV1_HD/Users/Shared/Software_current_cruise/MOD_fish_lib/bash_scripts/rsync_to_sharedrive_motive2025.sh
#
# (The asterisks say when and how often to run the file. This means run every day at 0600 1200 1800 0000.)
# crontab.guru - website to help you figure out how to schedule crontabs
