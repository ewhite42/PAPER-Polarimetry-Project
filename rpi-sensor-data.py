#rpi-sensor-data.py

import gpiozero
import socket
import time

def main():
    fftsize = 8192.0
    srate = 4000000.0 
    #decimation = 15000.0

    INT_TIME = 300 #this should be the same as burst length in the dsp program (in seconds)
    HEADERSIZE = 10

    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    s.bind(('192.168.10.130', 1234))
    s.listen(5)

    #read in voltages from the appropriate channels
    temp_0 = gpiozero.MCP3208(channel=0, differential=True, max_voltage=5.15) 
    hum_v = gpiozero.MCP3208(channel=2, differential=True, max_voltage=5.15)
    temp_3 = gpiozero.MCP3208(channel=4, differential=True, max_voltage=5.15) 

    while True:
        # now our endpoint knows about the OTHER endpoint.
        clientsocket, address = s.accept()
        print("Connection from {0} has been established.".format(address))
        print("\nSending Raspberry Pi sensor data to client...")

        while True:

            initial_time = time.time()
            tot_time = initial_time
            total_temp0_voltage = 0
            total_temp3_voltage = 0
            total_humidity_voltage = 0
            i = 1

            while time.time() - initial_time <= INT_TIME:
                #tot_time += time.time()

                diff0_v = temp_0.voltage #get temp sensor voltage
                total_temp0_voltage += diff0_v

                diff3_v = temp_3.voltage
                total_temp3_voltage += diff3_v

                curr_hum_v = hum_v.voltage #get humidity voltage
                total_humidity_voltage += curr_hum_v

                i = i + 1

            #average sensor voltages over the amount of time specified
            temp0_voltage = round((total_temp0_voltage / i), 5)
            temp3_voltage = round((total_temp3_voltage / i), 5)

            #convert temperature voltage to actual temperature
            #using crude calibration
            
            temp0_voltage = round(temp0_voltage*100, 5)
            temp3_voltage = round(temp3_voltage*100, 5)
            hum_voltage = round((total_humidity_voltage / i), 5)

            time_stamp = time.time() #time.ctime(tot_time / i)

            msg = "{0:<20}{1:<15}{2:<20}{3:<20}".format(temp0_voltage, temp3_voltage,
hum_voltage, time_stamp)
            msg = "{0:<{1}}".format(len(msg), HEADERSIZE)+msg

            clientsocket.send(bytes(msg,"utf-8"))


if __name__ == "__main__":
  main()

