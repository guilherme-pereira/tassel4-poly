/*
 * ExceptionUtils.java
 *
 * Created on May 25, 2003, 2:11 AM
 */

package net.maizegenetics.util;


import org.apache.log4j.Logger;
import org.apache.log4j.Level;


/**
 *
 * @author  terryc
 */
public final class ExceptionUtils {
    
    /**
     * Default maximun number of entries.
     */
    public static final int DEFAULT_MAX_NUMBER_ENTRIES = 50;
    
    
    /**
     * ExceptionUtils Constructor.
     */
    private ExceptionUtils() {
    }
    
    
    /**
     * This sends an exception's stack trace to the specified
     * logger.  The max parameter defines the maximun number of
     * entries to send to the logger.  This logs the default number
     * of entries defined by DEFAULT_MAX_NUMBER_ENTRIES at the
     * default priorty of debug.
     *
     * @param exception exception to log
     * @param logger logger to send entries
     */
    public static void logStackTrace(Exception exception, Logger logger) {
        logStackTrace(exception, DEFAULT_MAX_NUMBER_ENTRIES, logger, Level.DEBUG);
    }
    
    
    /**
     * This sends an exception's stack trace to the specified
     * logger.  The max parameter defines the maximun number of
     * entries to send to the logger.  Setting max to a negative
     * number sends all the entries in the stack track to the logger.
     *
     * @param exception exception to log
     * @param max maximun number of entries to log
     * @param logger logger to send entries
     * @param priorty priorty of log message
     */
    public static void logStackTrace(Exception exception, int max, Logger logger, Level priorty) {
        logger.log(priorty, getStackTrace(exception, max));
    }
    
    
    /**
     * This sends an exception's stack trace to the specified
     * logger.  The max parameter defines the maximun number of
     * entries to send to the logger.  Setting max to a negative
     * number sends all the entries in the stack track to the logger.
     *
     * @param exception exception to log
     * @param max maximun number of entries to log
     *
     * @return stack trace as string
     */
    public static String getStackTrace(Exception exception, int max) {
        
        StackTraceElement [] elements = exception.getStackTrace();
        
        int numEntries = elements.length;
        
        if ((max >= 0) && (max < numEntries)) {
            numEntries = max;
        }
        
        StringBuffer buffer = new StringBuffer();
        for (int i=0; i<numEntries; i++) {
            buffer.append(elements[i].toString());
            buffer.append("\n");
        }
        
        return buffer.toString();
        
    }
    
    
    /**
     * This sends an exception's cause history to the
     * specified logger.
     *
     * @param exception exception to log
     * @param logger logger to send entries
     * @param priorty priorty of log message
     */
    public static void logExceptionCauses(Throwable exception, Logger logger, Level priorty) {
        logger.log(priorty, getExceptionCausesWithClassNames(exception));
    }
    
    
    /**
     * This gets an exception's cause history.  Only the message
     * text is included (e.g. the exception's class names are not output).
     *
     * @param exception exception to log
     *
     * @return exception's cause history as a string
     */
    public static String getExceptionCauses(Throwable exception) {
        
        StringBuffer buffer = new StringBuffer();
        buffer.append(exception.getMessage());
        
        Throwable cause = exception.getCause();
        while (cause != null) {
            buffer.append("\n");
            buffer.append(cause.getMessage());
            Throwable temp = cause.getCause();
            cause = temp;
        }
        
        String str = buffer.toString();
        
        return str;
    }
    
    
    /**
     * This gets an exception's cause history.
     *
     * @param exception exception to log
     *
     * @return exception's cause history as a string
     */
    public static String getExceptionCausesWithClassNames(Throwable exception) {
        
        StringBuffer buffer = new StringBuffer();
        buffer.append(exception.toString());
        
        Throwable cause = exception.getCause();
        while (cause != null) {
            buffer.append("\n");
            buffer.append(cause.toString());
            Throwable temp = cause.getCause();
            cause = temp;
        }
        
        String str = buffer.toString();
        
        return str;
    }

}
