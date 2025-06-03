#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  9 11:29:33 2020
modified on Sat 07/30/2022
@author: aleboyer
"""

import numpy as np
import serial
import sys
import time

###### user input ########

port_name    = 'tty.usbserial-FTU7BXC3'
mission_name = "TFO2024"
vehicle_name = "FCTD1"
##########################

coef=int(input("Neutral is 60 \r\n Percent extension [0-100]:"))

## THIS IS FOR FCTD A

ser_write = serial.Serial('/dev/%s' % port_name,230400)  # open serial port
time.sleep(.1)
print("##connect to: \r\n")
print(ser_write.name)         # check which port was really used

utc_time=time.time()
str_time=time.strftime("%D %T", time.gmtime(utc_time))
print(f"\t!! {str_time} UTC: Setting actuator to {coef}% !!")
print("##setting actu \r\n")


str_cmd  = "actu.set %i\r\n" % coef
bin_cmd=str_cmd.encode()
ser_write.write(bin_cmd)
print(bin_cmd)
# read echo from epsi
ser_write.readline()
print("##done\r\n")

    
    
