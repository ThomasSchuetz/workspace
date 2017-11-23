# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 12:25:04 2017

@author: tsz
"""

import socket

import json

UDP_IP = "127.0.0.1"
UDP_PORT = 5005

sock = socket.socket(socket.AF_INET, # Internet
                     socket.SOCK_DGRAM) # UDP
sock.bind((UDP_IP, UDP_PORT))

while True:
    data, addr = sock.recvfrom(1024) # buffer size is 1024 bytes
    print "received message:", json.loads(data)
    test_dict = json.loads(data)
    
    print test_dict["number"]