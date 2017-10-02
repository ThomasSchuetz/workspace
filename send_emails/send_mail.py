#!/usr/bin/env python
# -*- coding: utf-8 -*-

import smtplib

def send_message(toaddrs, msg, 
                 fromaddr='sender@address.com', 
				 username='sender@address.com', 
				 password='yourPassword'):
    """
    This function sends a message from a certain email address to a recipient
    
    Parameters
    ----------
    toaddrs : string
        Recipient's address, e.g. joesmith@yahoo.com
    msg : string
        The message you want to send to 'toaddrs'. 
    fromaddr : string
        Sender's address, e.g. janesmith@gmail.com
    username : string
        Sender's identification, e.g. janesmith@gmail.com
    password : string
        Sender's email password, e.g. ***********
        This function does not forward passwords, etc. to a third party. I do 
        not know if something like this is done within smtplib, but this 
        function DOES NOT!
    """
    # Create a new server connection
	server = smtplib.SMTP('smtpAddress') # for gmail: smtp.gmail.com:587
	server.ehlo()
	server.starttls()
    
    # Log in
	server.login(username,password)
	
    # Send the email
    server.sendmail(fromaddr, toaddrs, msg)
	
    # Close the server connection
    server.quit()

if __name__ == "__main__":
    # This is a short example (of course with fake email addresses.)
    # Parameter definition
	toaddrs  = 'recipient@address.com'
	msg = 'Your message'
	fromaddr = 'sender@address.com'
	username = 'sender@address.com'
	password = 'yourPassword'

    # Send email
	send_message(toaddrs, msg, fromaddr, username, password)