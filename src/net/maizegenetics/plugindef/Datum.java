/*
 * Datum.java
 *
 */
package net.maizegenetics.plugindef;

import java.io.Serializable;

/**
 * This wraps data elements used as input or output to
 * Tassel modules.
 *
 * @author terryc
 */
public class Datum implements Serializable {

    private static final long serialVersionUID = -5197800047652332969L;
    
    private String myName;
    private final Object myData;
    private String myComment;

    /**
     * Creates a new instance of Datum
     *
     * @param name name of this data instance.
     * @param data data element
     * @param comment optional comment
     *
     */
    public Datum(String name, Object data, String comment) {

        if ((name == null) || (name.length() == 0)) {
            throw new IllegalArgumentException("Datum: init: name must be assigned.");
        }

        if (data == null) {
            throw new IllegalArgumentException("Datum: init: data can not be null.");
        }

        myName = name;
        myData = data;
        myComment = comment;

    }

    public String getName() {
        return myName;
    }

    public void setName(String s) {
        myName = s;
    }

    public Object getData() {
        return myData;
    }

    public String getComment() {
        return myComment;
    }

    public Class getDataType() {
        return myData.getClass();
    }

    public String toString() {
        return myName;
    }
}
