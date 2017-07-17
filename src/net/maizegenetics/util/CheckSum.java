package net.maizegenetics.util;

import java.io.*;
import java.security.MessageDigest;
import org.apache.log4j.Logger;

public class CheckSum {

    private static final Logger myLogger = Logger.getLogger(CheckSum.class);

    private CheckSum() {
        // utility class
    }

    public static String getMD5Checksum(String filename) {
        return getChecksum(filename, "MD5");
    }

    public static String getChecksum(String filename, String protocol) {

        try {
            InputStream inputStream = Utils.getInputStream(filename);
            MessageDigest digester = MessageDigest.getInstance(protocol);

            byte[] buffer = new byte[8192];
            int numOfBytesRead;
            while ((numOfBytesRead = inputStream.read(buffer)) > 0) {
                digester.update(buffer, 0, numOfBytesRead);
            }
            byte[] hashValue = digester.digest();
            return convertBytesToHex(hashValue);
        } catch (Exception ex) {
            myLogger.error(ex.getMessage());
        }

        return null;

    }

    private static String convertBytesToHex(byte[] bytes) {
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < bytes.length; i++) {
            builder.append(String.format("%02x", bytes[i] & 0xff));
        }
        return builder.toString();
    }
}
