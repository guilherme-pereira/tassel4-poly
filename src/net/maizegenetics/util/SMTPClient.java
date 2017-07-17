package net.maizegenetics.util;

import java.io.File;
import java.util.*;
import javax.mail.*;
import javax.mail.internet.*; 
import javax.activation.*;


/**
 * This class can be used to send email.  When the host is "appsmtp.mail.cornell.edu", any machine within the Cornell
 * network can be used to send an email.
 * Reference: https://confluence.cornell.edu/login.action?os_destination=https%3A%2F%2Fconfluence.cornell.edu%2Fdisplay%2Fcitapps%2FHow%2Bto%2Bconfigure%2Bapplications%2Bto%2Bsend%2Bmail%2Bthrough%2Bthe%2BCornell%2Bmail%2Bservers.
 */
public class SMTPClient { 

    private MimeMessage message;
    
    public SMTPClient(String host, String[] toAddresses){
        
        Properties properties = System.getProperties();
        properties.setProperty("mail.smtp.host", host);
        
        // Get the default Session object.
        Session session = Session.getDefaultInstance(properties);
        message = new MimeMessage(session); 

        for(int i = 0; i < toAddresses.length; i++){
            try{
                message.addRecipient(Message.RecipientType.TO, new InternetAddress(toAddresses[i]));
                message.setFrom(new InternetAddress(toAddresses[0]));  // just use first address as From
            }catch(javax.mail.internet.AddressException ae){
                ae.printStackTrace();
            }
             catch(javax.mail.MessagingException me){
                me.printStackTrace();
             }
        }
    }

    public void sendMessageWithAttachment(String subject, String msg, String fileAttachment) throws javax.mail.MessagingException{
        
        MimeBodyPart messageBodyPart = new MimeBodyPart();
        messageBodyPart.setText(msg);
        message.setSubject(subject);
        
        Multipart multipart = new MimeMultipart();
        multipart.addBodyPart(messageBodyPart);



        File aFileAttachment = new File(fileAttachment);
        if(aFileAttachment.canRead()){

            messageBodyPart = new MimeBodyPart();
            DataSource source = new FileDataSource(fileAttachment);
            messageBodyPart.setDataHandler(new DataHandler(source));
            messageBodyPart.setFileName(aFileAttachment.getName());
            multipart.addBodyPart(messageBodyPart);
        }else{
            messageBodyPart.setText(msg + "\n\nCould not attach file:\n" + fileAttachment);
        }
        // Put parts in message
        message.setContent(multipart);

        Transport.send(message);
        
    }

    public void sendMessage(String subject, String msg) throws javax.mail.MessagingException{
        message.setSubject(subject);
        message.setText(msg);
        Transport.send(message);
    }
}
