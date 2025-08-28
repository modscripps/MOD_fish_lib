#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 12:12:19 2020

@author: aleboyer
"""

#virtual epsi
import keyboard


import struct
import serial
import time
import numpy as np
import collections
import matplotlib.pyplot as plt

#def checksum(barray):
#    return functools.reduce(lambda x,y:x+y, barray)
def checksum(array):
    sum=bytes()
    for x in array:
        sum=sum or x
    return sum



ser = serial.Serial('COM12',115200)  # open serial port
ser.flush()


efe={
     "nbchannel":7,
     "channels":["t1","t2","s1","s2","a1","a2","a3"],
     "ADCbytes":3,
     "sampleperframe":10,
     "sampling_frequency":325
     }

efe.update({"delay":4*efe["sampleperframe"]/efe["sampling_frequency"]})
print("Delay: ", efe['delay'])

# efe.update({"delay":2}) # Debug: Reduce time of samples sent

EFE_sample={
        "t1":0x000000,\
        "t2":0x000000,\
        "s1":0x000000,\
        "s2":0x000000,\
        "a1":0x000000,\
        "a2":0x000000,\
        "a3":0x000000,\
        }

EFE_sample_struct=struct.Struct('3s3s3s3s3s3s3s')


EFE={
     "tag":b"$EFE",
     "timestamp":str.encode(hex(int(time.time()*1e6)),"utf8"),
     "payload_size":str.encode(hex(efe["sampleperframe"] * efe["ADCbytes"] * efe["nbchannel"]),"utf8"),
     "header_checksum":str.encode("*"+("%02x" % 0),"utf8"),
     "payload":bytearray(efe["nbchannel"]*efe["ADCbytes"]*efe["sampleperframe"]),
     "checksum":str.encode("*"+("%02x" % 0),"utf8"),
     "endblock":b"\r\n",
     }

frameformat=str(len(EFE["tag"]))+'s'+\
            str(len(EFE["timestamp"]))+'s'+\
            str(len(EFE["payload_size"])) + 's' +\
            str(len(EFE['header_checksum'])) + 's' +\
            str(len(EFE["payload"]))+'s'+\
            str(len(EFE["checksum"]))+'s'+\
            str(len(EFE["endblock"]))+'s'
                 
efe.update({"frameformat":frameformat})
print("Frame format: %s" % efe['frameformat'])
# time.sleep(3)
#EFE=(b"$$EFE,",
#     b"0x1e,",
#     str.encode(str(int(time.time()*1e6))+",","utf8"),
#     bytearray(efe["nbchannel"]*efe["ADCbytes"]*efe["sampleperframe"]),
#     str.encode("*"+("%02x" % 0),"utf8"),
#     b"\r\n",
#     )
#efe_block_struct = struct.Struct('6s5s16s3360s3s2s')
#efe_block_struct.pack(*EFE)


signal_freq=10
phase0=0
phase1=np.pi
y_sin=0

dt=.5
t_0=time.time()
sample=0
store_flag=1

line = 0
current_time = time.time()
previoustime = current_time
sampleCount = 0
varyingSignals = True

while 1:
        current_time = time.time()
        while store_flag==1:
            for cha in efe["channels"]:
                if (cha == "a1"):
                    EFE_sample[cha]= ((line % 2000) + 3000).to_bytes(3, "big")
                elif (cha == "a2") and varyingSignals:
                     EFE_sample[cha]= ((line % 2000) + 4000).to_bytes(3, "big")
                elif (cha == "t1") and varyingSignals:
                     EFE_sample[cha] = (int(1000*np.sin(2*np.pi*signal_freq*sample/efe["sampling_frequency"]
                                  +180)+5000) + int(abs(np.random.normal(scale=100)))).to_bytes(3, "big")
                elif (cha == "s1") and varyingSignals:
                     EFE_sample[cha] = int(5000 + np.random.normal(scale=1000)).to_bytes(3, "big")
                else:
                    EFE_sample[cha]=y_sin.to_bytes(3,'big')
                
            local_sample=EFE_sample_struct.pack(EFE_sample["t1"],EFE_sample["t2"],\
                                                EFE_sample["s1"],EFE_sample["s2"],\
                                                EFE_sample["a1"],EFE_sample["a2"],\
                                                EFE_sample["a3"])
            index_sample=(sample % efe["sampleperframe"])*efe["nbchannel"]*efe["ADCbytes"]
            EFE["payload"][index_sample:index_sample+ \
                            (efe["nbchannel"]*efe["ADCbytes"])]=local_sample

            y_sin=int(1000*np.sin(2*np.pi*signal_freq*sample/efe["sampling_frequency"]
                                  +phase0)+5000)
            
            line += 1

            sample=(sample+1)
            if sample % efe["sampleperframe"]==0:
                store_flag=0
      
        # _input = input("Blocking... Send bad packet? (y - send c - custom) \n")
        if (keyboard.is_pressed('y')):
             ser.write(b'*04241e2$24fsss')
        elif (keyboard.is_pressed('c')):
            badpacket = input("Type bad packet:")
            badpacket = badpacket.encode()
            ser.write(badpacket)
        EFE["timestamp"] = str.encode(hex(int(time.time()*1e6))+",","utf8")
        EFE["checksum"]  = str.encode("*"+("%02x" % checksum(EFE["payload"])),"utf8")
        EFE["header_checksum"] = str.encode("*"+("%02x" % checksum(EFE["timestamp"])), "utf8") # Checksum includes tag to end of payload size
        efe_block_struct = struct.pack(efe["frameformat"],\
                           EFE["tag"],\
                           EFE["timestamp"],\
                           EFE["payload_size"],\
                           EFE['header_checksum'],\
                           EFE["payload"],\
                           EFE["checksum"],\
                           EFE["endblock"]) 
        # with open(r'C:\Users\figir\OneDrive\Desktop\testpacket.txt', "wb") as file:
        #      file.write(efe_block_struct)
        for x in EFE:
             print("Key: ", x, " Value: ", EFE[x])
        # if current_time - previoustime >= 1:
        #      previoustime = current_time
        #      if not sampleCount:
        #         print("Samples: ", sampleCount / (current_time - previoustime))
        #      sampleCount = 0
        ser.write(efe_block_struct)

        # sampleCount += 1
        

        # ser.flush()
        time.sleep(efe["delay"])
        # time.sleep(.1)
        store_flag=1



        # if len(cb) >= 200:
        #     fig = plt.figure()
        #     ax1 = fig.add_subplot(111, label = "temp")
        #     ax1.plot(cb)
        #     plt.show()
        
        
        
        
        
        
        
        
