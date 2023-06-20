import datetime
import time
import thread
import pexpect

s1 = "A0:E6:F8:AE:14:01"
s2 = "54:6C:0E:52:FF:04"
s3 = "A0:E6:F8:AE:8F:86"
s4 = "54:6C:0E:53:12:64"
s5 = "54:6C:0E:78:C3:06"
s6 = "24:71:89:07:75:03"
s7 = "A0:E6:F8:AE:AA:04"

classifier = []

stream1 = pexpect.spawn('gatttool -I')
stream2 = pexpect.spawn('gatttool -I')
stream3 = pexpect.spawn('gatttool -I')
stream4 = pexpect.spawn('gatttool -I')
stream5 = pexpect.spawn('gatttool -I')
stream6 = pexpect.spawn('gatttool -I')
stream7 = pexpect.spawn('gatttool -I')

log = open("log.txt","w")

def input_thread(classifier):
        while True:
                c =  raw_input()
                classifier.append(c)


def ConnectSensors():
        
        stream1.sendline('connect '+s1)
        time.sleep(1)
        stream2.sendline('connect '+s2)
        time.sleep(1)
        stream3.sendline('connect '+s3)
        time.sleep(1)
        stream4.sendline('connect '+s4)
        time.sleep(1)
        stream5.sendline('connect '+s5)
        time.sleep(1)
        stream6.sendline('connect '+s6)
        time.sleep(1)
        stream7.sendline('connect '+s7)
        time.sleep(1)

        try:
                stream1.expect('Connection successful',timeout=5)
                stream2.expect('Connection successful',timeout=5)
                stream3.expect('Connection successful',timeout=5)
                stream4.expect('Connection successful',timeout=5)
                stream5.expect('Connection successful',timeout=5)
                stream6.expect('Connection successful',timeout=5)
                stream7.expect('Connection successful',timeout=5)                
                
        except pexpect.TIMEOUT, e:
                log.write("\n***ERROR***\n")
                log.write(e[0])
                log.write("\n")
                index = e[0].find("connect")+8
                print"There was an issue connecting to sensor",e[0][index:index+17]
                return False
        print"All sensors connected"

def DisconnectSensors():
        stream1.sendline('disconnect '+s1)
        time.sleep(1)
        stream2.sendline('disconnect '+s2)
        time.sleep(1)
        stream3.sendline('disconnect '+s3)
        time.sleep(1)
        stream4.sendline('disconnect '+s4)
        time.sleep(1)
        stream5.sendline('disconnect '+s5)
        time.sleep(1)
        stream6.sendline('disconnect '+s6)
        time.sleep(1)
        stream7.sendline('disconnect '+s7)

def SetParameters():
        stream1.sendline('char-write-cmd 0x3c 38:02') #turn on x,y,z accelerometers, magnetometer, gyro range = 8G
        time.sleep(.5)
        stream2.sendline('char-write-cmd 0x3c 38:02') #turn on x,y,z accelerometers, range = 8G
        time.sleep(.5)
        stream3.sendline('char-write-cmd 0x3c 38:02') #turn on x,y,z accelerometers, range = 8G
        time.sleep(.5)
        stream4.sendline('char-write-cmd 0x3c 38:02') #turn on x,y,z accelerometers, range = 8G
        time.sleep(1)
        stream5.sendline('char-write-cmd 0x3c 38:02') #turn on x,y,z accelerometers, magnetometer, gyro range = 8G
        time.sleep(.5)
        stream6.sendline('char-write-cmd 0x3c 38:02') #turn on x,y,z accelerometers, range = 8G
        time.sleep(.5)
        stream7.sendline('char-write-cmd 0x3c 38:02') #turn on x,y,z accelerometers, range = 8G
        time.sleep(.5)        

        stream1.sendline('char-write-cmd 0x3e 0A') 
        stream2.sendline('char-write-cmd 0x3e 0A') 
        stream3.sendline('char-write-cmd 0x3e 0A') 
        stream4.sendline('char-write-cmd 0x3e 0A')
        stream5.sendline('char-write-cmd 0x3e 0A') 
        stream6.sendline('char-write-cmd 0x3e 0A') 
        stream7.sendline('char-write-cmd 0x3e 0A') 
        time.sleep(1)

def TurnOnNotifs():
        stream1.sendline('char-write-cmd 0x3a 01:00') 
        stream2.sendline('char-write-cmd 0x3a 01:00') 
        stream3.sendline('char-write-cmd 0x3a 01:00') 
        stream4.sendline('char-write-cmd 0x3a 01:00')
        stream5.sendline('char-write-cmd 0x3a 01:00') 
        stream6.sendline('char-write-cmd 0x3a 01:00') 
        stream7.sendline('char-write-cmd 0x3a 01:00') 

def TurnOffNotifs():
        stream1.sendline('char-write-cmd 0x3a 00:00') 
        stream2.sendline('char-write-cmd 0x3a 00:00') 
        stream3.sendline('char-write-cmd 0x3a 00:00') 
        stream4.sendline('char-write-cmd 0x3a 00:00')
        stream5.sendline('char-write-cmd 0x3a 00:00') 
        stream6.sendline('char-write-cmd 0x3a 00:00') 
        stream7.sendline('char-write-cmd 0x3a 00:00')        
        
def main():
        print"QMACP Recording Script v2.1"
        print"created by John Rattray\n"

        print"Attempting to connect to sensors .."
        if (ConnectSensors()==False):
                log.close()
                return
        
                
        SetParameters()
        

        print "3"
        time.sleep(1)
        print "2"
        time.sleep(1)
        print "1"
        time.sleep(1)
        print "Recording..."
        print "Press Ctrl+C to end recording"

        TurnOnNotifs()        

        file = open("raw_data.txt","w")
        start = time.time()
        count = 0
        c = '1'
        thread.start_new_thread(input_thread, (classifier,))
        try:
                while True:
                        
                        if classifier:
                                c  = classifier[0]
                                classifier.pop()
                                                
                        try:
                                stream1.expect("value:", timeout=1)
                                ts1 = datetime.datetime.now().strftime("%H:%M:%S:%f")     
                                stream2.expect("value:", timeout=1)
                                ts2 = datetime.datetime.now().strftime("%H:%M:%S:%f")
                                stream3.expect("value:", timeout=1)
                                ts3 = datetime.datetime.now().strftime("%H:%M:%S:%f")
                                stream4.expect("value:", timeout=1)
                                ts4 = datetime.datetime.now().strftime("%H:%M:%S:%f")
                                stream5.expect("value:", timeout=1)
                                ts5 = datetime.datetime.now().strftime("%H:%M:%S:%f")     
                                stream6.expect("value:", timeout=1)
                                ts6 = datetime.datetime.now().strftime("%H:%M:%S:%f")
                                stream7.expect("value:", timeout=1)
                                ts7 = datetime.datetime.now().strftime("%H:%M:%S:%f")                                
                                file.write("1 "+stream1.buffer[1:54]+" "+c+" "+ts1+ "\n")
                                file.write("2 "+stream2.buffer[1:54]+" "+c+" "+ts2+ "\n")
                                file.write("3 "+stream3.buffer[1:54]+" "+c+" "+ts3+ "\n")
                                file.write("4 "+stream4.buffer[1:54]+" "+c+" "+ts4+ "\n")
                                file.write("5 "+stream5.buffer[1:54]+" "+c+" "+ts5+ "\n")
                                file.write("6 "+stream6.buffer[1:54]+" "+c+" "+ts6+ "\n")
                                file.write("7 "+stream7.buffer[1:54]+" "+c+" "+ts7+ "\n")                                
                        except pexpect.TIMEOUT, e:
                                log.write("\n***ERROR***\n")
                                log.write(e[0])
                                log.write("\n")
                                index = e[0].find("LE")-26
                                print"Connection was dropped with",e[0][index:index+17]
                                TurnOffNotifs()
                                DisconnectSensors()
                                print"Attempting to reconnect"                        
                                while (ConnectSensors()==False):
                                        print"Attempting to reconnect"
                                        DisconnectSensors()
                                SetParameters()
                                TurnOnNotifs()
                                print"Recording..."
                                
                        
        except KeyboardInterrupt:
                pass


        print time.time() - start
        TurnOffNotifs()
        
        file.close()
        log.close()

        ##convert all saved data
        file = open("raw_data.txt",'r')
        file_out = open("data.txt",'w')
        for line in file:

                gx = int(line [5:7]+line[2:4],16)
                if gx > 0x7FFF:
                        gx-=0x10000
                gx = gx/65536.0/500.0
                gy = int(line [11:13]+line[8:10],16)
                if gy > 0x7FFF:
                        gy-=0x10000
                gy = gy/65536.0/500.0
                gz = int(line[17:19]+line[14:16],16)
                if gz > 0x7FFF:
                        gz-=0x10000
                gz = gz/65536.0/500.0

                x = int(line [23:25]+line[20:22],16)
                if x > 0x7FFF:
                        x-=0x10000
                x = x/4096.0 #32768/2 since range is 8G
                y = int(line [29:31]+line[26:28],16)
                if y > 0x7FFF:
                        y-=0x10000
                y = y/4096.0 #32768/2 since range is 8G
                z = int(line[35:37]+line[32:34],16)
                if z > 0x7FFF:
                        z-=0x10000
                z = z/4096.0#32768/2 since range is 8G
                
                mx = int(line[41:43]+line[38:40],16)
                if mx > 0x7FFF:
                        mx-=0x10000

                my = int(line[47:49]+line[44:46],16)
                if my > 0x7FFF:
                        my-=0x10000

                mz = int(line[53:55]+line[50:52],16)
                if mz > 0x7FFF:
                        mz-=0x10000

                file_out.write (line[:2]+str(gx)+" "+str(gy)+" "+str(gz)+" "+str(x)+" "+str(y)+" "+str(z)+" "+str(mx)+" "+str(my)+" "+str(mz)+line[55:57]+" "+line[58:73]+"\n")

        file.close()
        file_out.close()        


main()
