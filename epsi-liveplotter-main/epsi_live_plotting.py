# attempt to import pyepsi
import serial 
import serial.tools.list_ports
import numpy as np
import struct 
import time
import pg_axisitem
# import matplotlib.pyplot as plt
from collections import deque
import pyqtgraph as pg
from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import QApplication
from scipy.signal import welch, detrend


"""
If I could refactor this, I would use multiprocessing to offload computation
Trust the sampling frequency - it's known, interpolate the data
    Using the timestamp given a certain number of points 
    Subtract the first timestamp, 

Look at the header for the real data,
! Detrend the data
Magic Numbers: nfft=1024, nperseg=1024
Offload the data processing (sample conversion, welch computation)
    Print data solely for updating the graph and gui with modified data
    24, 20, 16 bit noise floor lines in fft graph 
        Fn    = .5*FS;  % Nyquist frequency
        FR    = 2.5;    % Full range in Volts
        def_noise=@(x)((FR/2^x)^2 /Fn);
        (2.5 / 2^24)**2 / (2.5/2)

! add to version control (github)


[ ] RMS Value of data in voltage on the screen
[ ] Work on displaying fake data in order to help visualize what the software does 
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

        # Set up serial communication (Default)
        # self._ser = serial.Serial("/dev/cu.usbserial-FTE20C01", baud_rate, timeout=2)
        # Set up serial communication by promopting user
        ports = serial.tools.list_ports.comports()
        for idx, port in enumerate(sorted(ports)):
                print("%s: %s" % (idx, port))
        
        port_name = sorted(ports)[int(input("Select a port: "))]
        self._ser = serial.Serial(port_name.device, baud_rate, timeout=2)
        

        self._ser.flush()
        self._ser.write(b'\r')
        # self._ser.write(b'som.stop \r')  # Clear an instance if it's running
        # self._ser.write(b'som.epsi \r')
        self._ser.write(b'som.start \r')

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

        # set up curves for plotting
        QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
        self.win = pg.GraphicsLayoutWidget(show=True)
        self.win.ci.layout.setColumnStretchFactor(0,2)
        self.win.ci.layout.setColumnStretchFactor(1,3)
        self.p1 = self.win.addPlot(title='Accel. Data', labels={'bottom':'Count', 'left':'Volts'}, row=0, col=0, colspan=1)
        self.p1.showGrid(x=True, y=True, alpha=.5)
        self.p1.addLegend()
        self.curve1 = self.p1.plot(name='a1', pen='yellow')
        self.curve2 = self.p1.plot(name='a2', pen='crimson')
        self.curve3 = self.p1.plot(name='a3',pen='coral')
        # self.p1.setYRange(1, 10) # Sets the Y initial Y range of the program

        self.win.nextRow()

        self.p2 = self.win.addPlot(title='Shear', labels={'bottom':'Count', 'left':'Volts'}, row=1, col=0)
        self.p2.showGrid(x=True, y=True, alpha=.5)
        self.p2.addLegend()
        self.curve4 = self.p2.plot(name='s1', pen='lightgreen')
        self.curve5 = self.p2.plot(name='s2', pen='green')
        self.win.nextRow()
        # self.p2.setYRange(1, 10) # Sets the Y initial Y range of the program

        self.p3 = self.win.addPlot(title='Temperature', labels={'bottom':'Count', 'left':'Volts'}, row=2, col=0)
        self.p3.showGrid(x=True, y=True, alpha=.5)
        self.p3.addLegend()
        self.curve6 = self.p3.plot(name='t1', pen='lightblue')
        self.curve7 = self.p3.plot(name='t2', pen='blue')
        self.win.nextRow()
        # self.p3.setYRange(1, 10) # Sets the Y initial Y range of the program

        self.p4 = self.win.addPlot(title='FFT', labels={'bottom':'Frequency', 'left':'V² / Hz'} , row=0, rowspan=3, col=1)
        self.p4.setAxisItems({'bottom':pg_axisitem.CustomAxisItem('bottom')})
        self.p4.addLegend()
        self.legend = self.p4.legend
        self.p4.showGrid(x=True, y=True, alpha=.5)
        # self.p4.getAxis('bottom').setTickDensity(.1)
        
        self.p4.setLogMode(x=True, y=True) # set fft to log mode
        self.curve8 = self.p4.plot(name='a1', pen='yellow')
        self.curve9 = self.p4.plot(name='a2', pen='crimson')
        self.curve10 = self.p4.plot(name='a3', pen='coral')
        self.curve11 = self.p4.plot(name='s1', pen='lightgreen')
        self.curve12 = self.p4.plot(name='s2', pen='green')
        self.curve13 = self.p4.plot(name='t1', pen='lightblue')
        self.curve14 = self.p4.plot(name='t2', pen='blue')
        self.p4.setYRange(-5, -18, padding=.1)
        # self.bitnoise24 = self.p4.plot(name='24 Bit Noise', pen=pg.mkPen('w', width=3, style=QtCore.Qt.DashLine))
        # self.bitnoise20 = self.p4.plot(name='20 Bit Noise', pen=pg.mkPen('w', width=3, style=QtCore.Qt.DashLine))
        # self.bitnoise16 = self.p4.plot(name='16 Bit Noise', pen=pg.mkPen('w', width=3, style=QtCore.Qt.DashLine))


        # Noise floor
        self.bitnoise24Value = np.log10((2.5 / 2**24)**2 / (320/2))
        self.bitnoise20Value = np.log10((2.5 / 2**20)**2 / (320/2))
        self.bitnoise16Value = np.log10((2.5 / 2**16)**2 / (320/2))
        self.ADXLnoise = np.log10((20e-6)**2)

        # Add noise plots to graph
        # https://github.com/pyqtgraph/pyqtgraph/issues/22 There is a bug with infiniteline, hence the log10 for noise floor calculations
        vline = pg.InfiniteLine(angle=0, name='16 Bit Noise', pen=pg.mkPen('w', width=1, style=QtCore.Qt.DashLine), pos=self.bitnoise16Value, 
                                label='')
        self.p4.addItem(vline, ignoreBounds=True)
        vline = pg.InfiniteLine(angle=0, name='20 Bit Noise', pen=pg.mkPen('w', width=1, style=QtCore.Qt.DashLine), pos=self.bitnoise20Value,
                                label='')
        self.p4.addItem(vline, ignoreBounds=True)
        vline = pg.InfiniteLine(angle=0, name='24 Bit Noise', pen=pg.mkPen('w', width=1, style=QtCore.Qt.DashLine), pos=self.bitnoise24Value,
                                label='')
        self.p4.addItem(vline, ignoreBounds=True)
        vline = pg.InfiniteLine(angle=0, name='ADXL Accel Noise', pen=pg.mkPen('w', width=1, style=QtCore.Qt.DashLine), pos=self.ADXLnoise,
                                label='')
        self.p4.addItem(vline, ignoreBounds=True)

        

        self.curves = [self.curve6, self.curve7, self.curve4, self.curve5, self.curve1, self.curve2, self.curve3]
        self.psd = [self.curve13, self.curve14, self.curve11, self.curve12, self.curve8, self.curve9, self.curve10]

        # # Debug for GUI so you don't have to have a real serial port
        # pg.QtGui.QGuiApplication.instance().exec_()

        # Set up constants
        self.acc_offset = 1.8/2
        self.acc_factor = .4
        self.gain = 1

    def close(self):
        self._ser.write(b'som.stop \r')
        self._ser.close()
        print('Terminating...')

    def calculate_rms_value(self, data, idx):
        """
        Calculates the RMS values of each channel and updates the legend to display the values
        Must be RMS for time series
        """
        return
        rms = np.sqrt(np.mean(data**2))
        self.legend.items[idx-1][1].setText("%s- %s" % (self.channels[idx], rms))

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
        self.counter += 1 # Computation counter for profiling 
        for idx, x in enumerate(self.channel_buffers):    # this is not very efficient find way to optimize
          if (x == 'time'): # Quick fix 
              continue
          dat = np.zeros(len(self.channel_buffers[x]))
          for index, samples in enumerate(self.channel_buffers[x]):
            dat[index] = samples
          self.curves[idx-1].setData(y=dat)
          dat = dat[-self.window:]
          freq, psd = welch(x=dat, fs=self.samplerate, nfft=1024, nperseg=1024) # np.array(dat[-self.window:])
          if (self.counter % self.samplerate*3) == 0:
            self.calculate_rms_value(psd, idx)
          self.psd[idx-1].setData(x=freq, y=psd)
        # detrending the data
        # dat = detrend()
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
                idx = self._data.find(b'$EFE4') # $ or header "EFE4"
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
        printing_data = True
        for data in frame:
            for idx, data_value in enumerate(data):
                # Conversion to bytes and bipolar/unipolar processing
                if (self.counter % (self.samplerate / 16)) == 0 and printing_data:
                    print(" " + str(self.channels[idx]) + " : " + data_value.hex(), end ="")
                    if(self.channels[idx] == self.channels[-1]):
                        print("")
                        printing_data = False
                data_value = int.from_bytes(data_value)
                if self.channels[idx] in {'t1', 't2'}: # Unipolar
                    data_value = 2.5 / self.gain * data_value / 2**24 # FR = 2.5
                elif self.channels[idx] in {'s1', 's2'}: # Bipolar
                    data_value = 2.5 / self.gain * (data_value / 2**(23) - 1)
                elif self.channels[idx] in {'a1', 'a2', 'a3'}: # Convert for g
                     data_value = 1.8 / self.gain * data_value / 2**24 # Convert to unipolar
                    #  data_value = (data_value - self.acc_offset) / self.acc_factor     # Acc_factor = .4 | acc_offset = 1.8/2 convert for g
                self.channel_buffers[self.channels[idx]].append(data_value)
        
    
    def stream(self):
        thread = QtCore.QThread()
        thread.started.connect(self.run)
        thread.start()
        pg.QtGui.QGuiApplication.instance().exec_()

try:
    epsi = epsi_stream()
    print("Epsi_stream initialized")
    epsi.run()
except KeyboardInterrupt:
    epsi.close()


# # profile_run()
# cProfile.run('profile_run()', sort='cumulative')
