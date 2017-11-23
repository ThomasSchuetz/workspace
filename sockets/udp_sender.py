# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 12:25:30 2017

@author: tsz
"""

import socket

import numpy as np
import json


test_dict = {"bool_true": str(True),
             "bool_false": str(False),
             "number": 3.5,
             "list": [1, 2, "a", "c", 5],
             "array": list(np.random.rand(10)),
             "string": "Hello World!",
             1: 15}

UDP_IP = "127.0.0.1"
UDP_PORT = 5005
#MESSAGE = "Hello, World!"
MESSAGE = json.dumps(test_dict)

print "UDP target IP:", UDP_IP
print "UDP target port:", UDP_PORT
print "message:", str(MESSAGE)

sock = socket.socket(socket.AF_INET, # Internet
                     socket.SOCK_DGRAM) # UDP
sock.sendto(MESSAGE, (UDP_IP, UDP_PORT))

