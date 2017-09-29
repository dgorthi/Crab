#! /usr/bin/python

import corr
import time


#HOST = 'rpi-crab.casper.pvt'

HOST = 'rpicrab'
#HOST = "192.168.1.51"
DEST_IP   = (10<<24) + (10<<16) + (12<<8) + (35<<0)
DEST_MAC  = 0xe41d2d7dbac1
DEST_PORT = 10000
#ACC_LEN   = 2**18 # 2**18 is ~1 second
#ACC_LEN   = 2**8 #998.6438095238us, 2**8*8192/2100000000.0
ACC_LEN   = 2**11 #7989.150476us,2**11*8192/2100000000.0
#BOFFILE   = 'adc5g_spec_2017-02-01_1034.bof'
#BOFFILE   = 'adc5g_spec_2017-02-06_1610.bof'
#BOFFILE   = 'crab_20170208.bof'

BOFFILE   = 'crab_20170221.bof'
SRC_IP   = (10<<24) + (10<<16) + (12<<8) + (100<<0)
SRC_MAC  = 0x001122334455

print 'Connecting to board:', HOST
fpga = corr.katcp_wrapper.FpgaClient(HOST)
time.sleep(0.1)

print 'Programming!'
fpga.progdev(BOFFILE)
print 'done'

print 'Board clock (in MHz):', fpga.est_brd_clk()

# reset the boards
fpga.write_int('rst', 1)

# configure tge stuff
fpga.write_int('eth_dest_ip', DEST_IP)
fpga.write_int('eth_dest_port', DEST_PORT)


arp_table = [DEST_MAC] * 256

fpga.config_10gbe_core('eth_ten_Gbe_v2', SRC_MAC, SRC_IP, DEST_PORT, arp_table)

fpga.write_int('filterbank_fft_shift',8191)
# Set accumulation length
fpga.write_int('vacc_acc_len', ACC_LEN)

# deassert reset
fpga.write_int('rst', 0)

# trigger PPS
fpga.write_int('pps_arm', 0)
fpga.write_int('pps_arm', 1)
fpga.write_int('pps_arm', 0)

fpga.write_int('sw_pps', 0)
fpga.write_int('sw_pps', 1)
fpga.write_int('sw_pps', 0)


# enable 10gbe output
fpga.write_int('eth_tge_en', 1)

#fpga.write_int('tvg_cmult_en',1)

#fpga.write_int('vacc_shift',56)
fpga.write_int('vacc_shift',48)





