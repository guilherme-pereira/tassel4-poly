/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.annotation;

import java.lang.annotation.Documented;

/**
 * Experimental algorithms or classes that should be used in production code.  These
 * are new ideas that have not been validated.  No one should use these classes without
 * reading and understanding the underlying code.  Please contact the owner and
 * creator to find out more
 * 
 * @author Ed Buckler
 */
@Documented
public @interface Experimental {
    String owner();
    String value();
    String createdBy();
    String lastChanged();
}
