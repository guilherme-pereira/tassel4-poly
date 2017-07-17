
package net.maizegenetics.pal.ids;

import com.google.common.collect.HashMultimap;

/**
 * Identifiers normally only track a name, while the AnnotatedIdentifier provides support
 * for numerous other characteristics about the taxon (sample).  This could be used
 * to store information on secondary IDs, species, genus, parents, family, other sets, 
 * expected inbreeding coefficients, etc.
 * 
 * @author  Ed Buckler
 * @deprecated taxa package has the annotated replacements
 */
@Deprecated
public class AnnotatedIdentifier extends Identifier  {
    private HashMultimap<String, String> textAnnotations=null;
    private HashMultimap<String, Double> quantAnnotations=null;

    public AnnotatedIdentifier(String name) {
        super(name);
        textAnnotations=HashMultimap.create(10,2);
        quantAnnotations=HashMultimap.create(10,2);
    }
    
    public AnnotatedIdentifier(String name, HashMultimap<String, String> textAnnotations, 
            HashMultimap<String, Double> quantAnnotations) {
        super(name);
        this.textAnnotations=textAnnotations;
        this.quantAnnotations=quantAnnotations;
    }
    
    public void addAnnotation(String annoName, String value) {
        textAnnotations.put(annoName, value);
    }
    
    public void addAnnotation(String annoName, double value) {
        quantAnnotations.put(annoName, value);
    }


    public Object[] getAnnotation(String annoName) {
       if(textAnnotations.containsKey(annoName)) return textAnnotations.get(annoName).toArray(new String[0]);
       if(quantAnnotations.containsKey(annoName)) return quantAnnotations.get(annoName).toArray(new Double[0]);
       return null;
    }

    public String[] getTextAnnotation(String annoName) {
        return textAnnotations.get(annoName).toArray(new String[0]);
    }

    public Double[] getQuantAnnotation(String annoName) {
        return quantAnnotations.get(annoName).toArray(new Double[0]);
    }
    
    public String toString() {
        StringBuilder sb=new StringBuilder(super.toString());
        for (String key : textAnnotations.keySet()) {
            for (String value : textAnnotations.get(key)) {
               sb.append("\t");
               sb.append("<");
               sb.append(key);
               sb.append(">");
               sb.append(value); 
            }
            
        }
        for (String key : quantAnnotations.keySet()) {
            for (Double value : quantAnnotations.get(key)) {
            sb.append("\t");
            sb.append("<#");
            sb.append(key);
            sb.append(">");
            sb.append(value); 
            }
        }
        return sb.toString();
    }
    
    
}
