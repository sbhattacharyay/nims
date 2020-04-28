import datetime
import time
import thread
import pexpect

s3 = "A0:E6:F8:AE:8F:86"  #RW

classifier = []

stream3 = pexpect.spawn('gatttool -I')

log = open("logRW.txt","w")

def input_thread(classifier):
        while True:
                c =  raw_input()
                classifier.append(c)


def ConnectSensors():
        
        stream3.sendline('connect '+s3)
        time.sleep(1)

        try:
                stream3.expect('Connection successful',timeout=5)
                
        except pexpect.TIMEOUT, e:
                log.write("\n***ERROR***\n")
                log.write(e[0])
                log.write("\n")
                index = e[0].find("connect")+8
                print"There was an issue connecting to sensor",e[0][index:index+17]
                return False
        print"RW sensor connected"

def DisconnectSensors():
        stream3.sendline('disconnect '+s3)
        time.sleep(1)

def SetParameters():
        stream3.sendline('char-write-cmd 0x3c 38:02') #turn on x,y,z accelerometers, range = 8G
        time.sleep(.5)

        stream3.sendline('char-write-cmd 0x3e 0A') 
        time.sleep(1)

def TurnOnNotifs():
        stream3.sendline('char-write-cmd 0x3a 01:00') 

def TurnOffNotifs():
        stream3.sendline('char-write-cmd 0x3a 00:00') 
        
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

        file = open("raw_dataRW.txt","w")
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
                                stream3.expect("value:", timeout=1)
                                ts3 = datetime.datetime.now().strftime("%H:%M:%S:%f")
                                file.write("3 "+stream3.buffer[1:54]+" "+c+" "+ts3+ "\n")
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
        file = open("raw_dataRW.txt",'r')
        file_out = open("dataRW.txt",'w')
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
