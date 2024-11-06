#!/bin/sh

rsync -auq --delete /Volumes/Software_TFO2024/ /Volumes/science_party_share/RR2410/MOD-SIO/08_software_used_in_cruise/

# Save this script on DEV3
#
# To make the file executable
# chmod a+x copy_raw_to_DEV3
#
# To add to the crontab
# crontab -e
#
# Type the following and save with :w
# 0 0/6 * * * /Users/Shared/FCTD_EPSI/copy_raw_to_DEV3
#
# (The asterisks say when and how often to run the file. This means run every day at 0600 1200 1800 0000.)
# crontab.guru - website to help you figure out how to schedule crontabs
