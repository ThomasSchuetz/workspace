#!/usr/bin/env python
# -*- coding: utf-8 -*-

import smtplib

def send_message(toaddrs, msg, 
                 fromaddr='sender@address.com', 
				 username='sender@address.com', 
				 password='yourPassword'):
	server = smtplib.SMTP('smtpAddress') # for gmail: smtp.gmail.com:587
	server.ehlo()
	server.starttls()
	server.login(username,password)
	server.sendmail(fromaddr, toaddrs, msg)
	server.quit()

if __name__ == "__main__":
	toaddrs  = 'recipient@address.com'

	msg = 'Your message'

	fromaddr = 'sender@address.com'
	username = 'sender@address.com'
	password = 'yourPassword'
	
	send_message(toaddrs, msg, fromaddr, username, password)