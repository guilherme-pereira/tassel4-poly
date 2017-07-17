package net.maizegenetics.pd;

import net.maizegenetics.pal.alignment.Locus;

import java.util.List;
import java.util.Set;
import java.util.TreeMap;

/**
 * User: dkroon
 * Date: 6/26/13
 * Time: 2:07 PM
 */
public class SiteAnnotation {

    private int site;
    private int chromosome;
    private int position;
    private TreeMap<String, Double> gwasResults;
    public TreeMap<String, Double> annotations;

    public SiteAnnotation(int aSite, int aChromosome, int aPosition ){
        site = aSite;
        chromosome = aChromosome;
        position = aPosition;
    }


    public void setGwasResults(TreeMap<String, Double> gwasResultsIn){
        gwasResults = gwasResultsIn;
    }

    public void setAnnotations(TreeMap<String, Double> annotationsIn){
        annotations = annotationsIn;
    }

    /**
     * Obtain complete listing of all available traits for which GWAS traits may be available
     * @return
     */
    public List<String> getGWASTraits(){
        List<String> list = null;
        Set<String> traits = this.gwasResults.keySet();

        if(traits != null){
            list.addAll(traits);
        }
        return list;
    }

    /**
     * Obtain GWAS result value for given trait
     * @param trait
     * @return
     */
    public Double getGWASResult(String trait){
        Double result = null;
        if(this.gwasResults.containsKey(trait)){
            result = this.gwasResults.get(trait);
        }
        return result;
    }

    /**
     * Obtain complete listing of all annotations types available. Examples: Minor Allele Frequency
     * @return
     */
    public List<String> getAnnotationTypes(){
        List<String> annotationTypes = null;

        Set<String> annotations = this.annotations.keySet();

        if(annotations != null){
            annotationTypes.addAll(annotations);
        }
        return annotationTypes;
    }


    /**
     *
     * @param annotationType
     * @return
     */
    public Double getAnnotationValue(String annotationType){
        Double value = null;
        if(this.annotations.containsKey(annotationType)){
            value = annotations.get(annotationType);
        }
        return value;
    }


}
