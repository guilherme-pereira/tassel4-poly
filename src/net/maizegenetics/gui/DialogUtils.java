/*
 * DialogUtils.java
 *
 * Created on April 17, 2006
 */

package net.maizegenetics.gui;


import java.awt.Component;

import javax.swing.JOptionPane;

import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;




/**
 *
 * @author terryc
 */
public class DialogUtils {
    
    private DialogUtils() {
    }
    
   
    public static void showError(String str, Component parent) {
        JOptionPane.showMessageDialog(parent, str, "Error", JOptionPane.ERROR_MESSAGE);
    }
    
    
    public static void showError(Throwable e, Component parent) {
        String str = Utils.shortenStrLineLen(ExceptionUtils.getExceptionCauses(e), 50, 7);
        JOptionPane.showMessageDialog(parent, str, "Error", JOptionPane.ERROR_MESSAGE);
    }
    
    
    public static void showErrorCause(Throwable e, Component parent) {
        
        Throwable temp = e.getCause();
        
        if (temp != null) {
            showError(temp, parent);
        }
        
        else {
            showError(e, parent);
        }
        
    }
    
}
