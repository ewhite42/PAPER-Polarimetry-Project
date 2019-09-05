#! /usr/bin/python
import socket
import datetime
import time

def rpiClient():
    HEADERSIZE = 10
    interval = 86400

    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect(('192.168.10.130', 1234))
    #today = time.gmtime(time.clock())
    fname_date = time.strftime("%m_%d_%Y_%H-%M-%S", time.gmtime())

    ofName = "RPI_DATA_"+fname_date+".txt"
    outFile = open(ofName, 'w')

    print("\nAccepting data from the RPI, saving to file {}...\n".format(ofName))
    print("\n{0:<20}{1:<15}{2:<20}{3:<20}".format("Balun Temp (K)","Rx Temp (K)","Humidity Voltage","Time Stamp"))
    start_time = time.time()
    while True:
        full_msg = ''
        new_msg = True
        while (time.time() - start_time) < interval:
            msg = s.recv(128)
            if new_msg:
                #print("new msg len:",msg[:HEADERSIZE])
                msglen = int(msg[:HEADERSIZE])
                new_msg = False

            #print("full message length: {0}".format(msglen))

            full_msg += msg.decode("utf-8")

            #print(len(full_msg))

            if len(full_msg)-HEADERSIZE == msglen:
                #print("full msg recvd")
                print(full_msg[HEADERSIZE:])
                outFile.write(full_msg[HEADERSIZE:]+"\n")
                new_msg = True
                full_msg = ""

    print "Observation complete!"

def main():
    rpiClient()

if __name__ == "__main__":
    main()
