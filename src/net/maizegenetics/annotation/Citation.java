/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.annotation;

import java.lang.annotation.Documented;

/**
 * Provides the publication citation information for an algorithm.  Ultimately, this
 * citation should be passed on to the resulting datasets, so that inventors of cool
 * algoritms get credit for their work.
 * @author edbuckler
 */
@Documented
public @interface Citation {
     String value() default "Publication in process please contact TASSEL group for citation";
//    String journalCitation();
//    boolean pubEmbargoed();
//    String pubEmbargoDate();
    
    
}
