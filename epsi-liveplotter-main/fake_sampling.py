# attempt to import pyepsi
import serial 
import numpy as np
import struct 
import time

# import matplotlib.pyplot as plt
from collections import deque
import pyqtgraph as pg
from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import QApplication
from scipy.signal import welch, detrend
import cProfile



"""
- ! Work on aspect ratio subplots
- Channel 5,6,7 on top
- Channel 1,2
- Channel 3,4
- Fix X-Axis
Use multiprocessing to to reduce packet loss
Trust the sampling frequency - it's known, interpolate the data
    Using the timestamp given a certain number of points 
    Subtract the first timestamp, 

Look at the header for the real data,
! Detrend the data
fix sbplot sizes (fft focs)
Magic Numbers: nfft=1024, nperseg=1024
Offload the data processing (sample conversion, welch computation)
    Print data solely for updating the graph and gui with modified data
    24, 20, 16 bit noise floor lines in fft graph 
        Fn    = .5*FS;  % Nyquist frequency
        FR    = 2.5;    % Full range in Volts
        def_noise=@(x)((FR/2^x)^2 /Fn);

! add to version control (github)

"""



class epsi_stream():
    # XOR checksum
    def checksum(self, arr):
        checksum = 0
        for i in range(len(arr)):
            checksum ^= arr[i]
        return checksum

    def __init__(self,
                 baud_rate = 230400,
                 channels = ["time","t1","t2","s1","s2","a1","a2","a3"],
                 frame_format = '5s16s8s3s2320s3s2s',
                 sample_rate=320, 
                 bytes_per_sample = 3
                 ):
        
        # Debug / profiler params
        self.counter = 0

        # Set up serial communication
        self._ser = serial.Serial("COM13", baud_rate, timeout=2)
        self._ser.flush()
        # Strip the frame format into individual packets
        # [0]: Tag 
        # [1]: timestamp 
        # [2]: payload_size
        # [3]: header_checksum
        # [4]: payload
        # [5]: checksum
        # [6]: endblock
        numberOfChanels = len(channels)
        self.frame_format = list(map(int, frame_format.strip("s").split("s")))
        sample_format = '8s' + (str(bytes_per_sample) + 's') * (numberOfChanels-1)
        self.packet_header = b"$EFE4"
        
        # Initialize parameters
        self._sample_format = struct.Struct(sample_format) #8s3s3s3s3s3s3s3s timestamp, 7 channels, 3 bytes each 
        self._framestruct = struct.Struct(frame_format)
        self.channels = channels
        self.channel_buffers = {channel: deque(maxlen=(15*sample_rate)) for channel in self.channels}

        self.samplerate = sample_rate
        self.window = sample_rate*10
        
        # Set state values:
        # State 0: Reading / no packet
        # State 1: Header validated
        # State 2: Payload Validated
        # State 3: Bad packet
        self.state = 0

        # Plotting stuff
        self.win = pg.GraphicsLayoutWidget(show=True)
        self.p1 = self.win.addPlot(title='Accel. Data') #row=0, col=0, colspan=1)
        self.curve1 = self.p1.plot(name='curve1', pen='yellow')
        self.curve2 = self.p1.plot(name='curve2', pen='crimson')
        self.curve3 = self.p1.plot(name='curve3',pen='coral')
        self.win.nextRow()

        self.p2 = self.win.addPlot(title='Strain') # , row=0, col=1, colspan=1
        self.curve4 = self.p2.plot(name='curve4', pen='orange')
        self.curve5 = self.p2.plot(name='curve5', pen='red')
        self.win.nextRow()

        self.p3 = self.win.addPlot(title='Temperature') # , row=0, col=2, colspan=1
        self.curve6 = self.p3.plot(name='curve6', pen='darkolivegreen')
        self.curve7 = self.p3.plot(name='curve7', pen='lime')
        self.win.nextRow()

        self.p4 = self.win.addPlot(title='FFT') # , row=3, rowspan=3, colspan=3
        self.p4.setLogMode(x=True, y=True) # set fft to log mode
        self.curve8 = self.p4.plot(name='curve8', pen='yellow')
        self.curve9 = self.p4.plot(name='curve9', pen='crimson')
        self.curve10 = self.p4.plot(name='curve10', pen='coral')
        self.curve11 = self.p4.plot(name='curve11', pen='orange')
        self.curve12 = self.p4.plot(name='curve12', pen='red')
        self.curve13 = self.p4.plot(name='curve13', pen='darkolivegreen')
        self.curve14 = self.p4.plot(name='curve14', pen='lime')

        self.curves = [self.curve6, self.curve7, self.curve4, self.curve5, self.curve1, self.curve2, self.curve3]
        self.psd = [self.curve13, self.curve14, self.curve11, self.curve12, self.curve8, self.curve9, self.curve10]


    def validate_header(self): # Checksum includes tag to end of payload size
        header_sum = self.packet[(sum(self.frame_format[0:3])):(sum(self.frame_format[0:4]))]
        if len(header_sum) <= 2  or ord(b'*') != header_sum[0] or not header_sum[1:].isalnum() or (self.packet[:self.frame_format[0]] != self.packet_header and len(self.packet) >= self.frame_format[1]):
            if len(self.packet) > sum(self.frame_format[0:4]): # Check length of packet to the sum of the frame
                self.state = 3
                print('Bad packet failed len sum check in header')
            return False
        arr = self.packet[0:sum(self.frame_format[0:3])]
        if int(header_sum[1:],16) == self.checksum(arr): 
            self.state = 1
            return True
        else:
            print('Bad packet failed checksum matching header')
            self.state = 3
    
    def validate_payload(self):
        payload_sum = self.packet[(sum(self.frame_format[0:5])):(sum(self.frame_format[0:6]))]
        if len(payload_sum) <= 2 or ord(b'*') != payload_sum[0] or not payload_sum[1:].isalnum():
            if len(self.packet) > sum(self.frame_format[0:6]):
                print('Bad packet failed len sum check in payload')
                self.state = 3
            return False
        arr = self.packet[sum(self.frame_format[0:4]):sum(self.frame_format[0:5])]
        if int(payload_sum[1:],16) == self.checksum(arr): 
            self.state = 2
            return True
        else:
            print('Bad packet failed checksum matchin in payload')
            self.state = 3

    def validate_data(self):
        if (self.state < 1):
            self.validate_header()
            return False
        if (self.state < 2):
            self.validate_payload()
            return False
        if len(self.packet) >= sum(self.frame_format) and self.state < 3:
            self._data = self.packet[:sum(self.frame_format)]
            self.parse_data()
            self.state = 0
            return True
        return False

    def print_data(self): # Uses idx-1 as a bandage. Refactoring probably good idea
        for idx, x in enumerate(self.channel_buffers):    # this is not very efficient find way to optimize
          if (x == 'time'): # Quick fix 
              continue
          dat = np.zeros(len(self.channel_buffers[x]))
          for index, samples in enumerate(self.channel_buffers[x]):
            dat[index] = int.from_bytes(samples)
          self.curves[idx-1].setData(y=dat)
          dat = dat[-self.window:]
          freq, psd = welch(x=dat, fs=self.samplerate, nfft=1024, nperseg=1024) # np.array(dat[-self.window:])
          self.psd[idx-1].setData(x=freq, y=psd)


        # detrending the data
        # dat = detrend()
        
        self.counter+=1
        print(self.counter)
        if (self.counter >= 1500):
            quit()
        QApplication.processEvents()

    def run(self):
        self.full_packet = False
        # self.count = 0
        self.packet = bytearray()
        # Read a line from the serial object
        while (not self.full_packet):
            self._data = self._ser.read(self._ser.in_waiting)   #self._ser.in_waiting
            if self.packet_header in self._data:
                self.packet = bytearray()
                idx = self._data.find(b'$EFE4') # $ or eader "EFE4"
                self.packet.extend(self._data[idx:])
                while (not self.validate_data()):
                    if (self.state >= 3):
                        self.state = 0
                        break
                    self._data = self._ser.read(self._ser.in_waiting) #self._ser.in_waiting
                    self.packet.extend(self._data)
                    self.print_data()

    def parse_data(self):
        self._data = self._framestruct.unpack(self._data)
        frame = self._sample_format.iter_unpack(self._data[4])
        for data in frame:
            for idx, data_value in enumerate(data):
                self.channel_buffers[self.channels[idx]].append(data_value)

    
    def stream(self):
        thread = QtCore.QThread()
        thread.started.connect(self.run)
        thread.start()
        pg.QtGui.QGuiApplication.instance().exec_()

# epsi = epsi_stream()
# print("Epsi_stream initialized")
# epsi.stream()


epsi = epsi_stream()
print("Epsi_stream initialized")
epsi.run()

# # profile_run()
# cProfile.run('profile_run()', sort='cumulative')
