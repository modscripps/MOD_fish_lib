#!/bin/sh

rsync -auq --delete /Users/Shared/EPSI_PROCESSING/Current_Cruise/Processed/ /Volumes/sci/shipside/SKQ202417S/MOD/05_processed_data/

# Save this script on DEV3
#
# To make the file executable
# chmod a+x copy_raw_to_DEV3
#
# To add to the crontab
# crontab -e
#
# Type the following and save with :w
# 0 0/6 * * * /Volumes/MOD HD/Users/Shared/Software_current_cruise/MOD_fish_lib/bash_scripts/rsync_to_sharedrive.sh
#
# (The asterisks say when and how often to run the file. This means run every day at 0600 1200 1800 0000.)
# crontab.guru - website to help you figure out how to schedule crontabs
