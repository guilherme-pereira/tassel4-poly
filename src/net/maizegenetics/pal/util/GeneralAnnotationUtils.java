package net.maizegenetics.pal.util;

import java.util.ArrayList;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: edbuckler
 * Date: 8/23/13
 * Time: 4:55 PM
 * To change this template use File | Settings | File Templates.
 */
public final class GeneralAnnotationUtils {
    public static final Object[] EMPTYOBJ=new Object[0];
    public static final String[] EMPTYSTR=new String[0];
    public static final double[] EMPTYdouble=new double[0];

    private GeneralAnnotationUtils() {
    }

    public static final Object[] getAnnotation(Map.Entry<String, String>[] myVariantsAndAnno, String annoName) {
        ArrayList<Object> result=new ArrayList<>(1);
        for (Map.Entry<String, String> me: myVariantsAndAnno) {
            if(me.getKey().equals(annoName)) result.add(me.getValue());
        }
        return result.toArray(EMPTYOBJ);
    }

    public static final String[] getTextAnnotation(Map.Entry<String, String>[] myVariantsAndAnno, String annoName) {
        ArrayList<String> result=new ArrayList<>(1);
        for (Map.Entry<String, String> me: myVariantsAndAnno) {
            if(me.getKey().equals(annoName)) result.add(me.getValue());
        }
        return result.toArray(EMPTYSTR);
    }

    public static final double[] getQuantAnnotation(Map.Entry<String, String>[] myVariantsAndAnno, String annoName) {
        try{ArrayList<Double> result=new ArrayList<>(1);
            for (Map.Entry<String, String> me: myVariantsAndAnno) {
                if(me.getKey().equals(annoName)) {
                    result.add(Double.parseDouble(me.getValue()));
                }
            }
            if(result.size()==0) return EMPTYdouble;
            double[] d=new double[result.size()];
            for (int i = 0; i < result.size(); i++) {
               d[i]=result.get(i);
            }
            return d;
        }catch(Exception e) {
            return EMPTYdouble;
        }
    }

    public static final String getConsensusAnnotation(Map.Entry<String, String>[] myVariantsAndAnno, String annoName) {
        return null;
        //  return myGA.getConsensusAnnotation(annoName);
    }

    public static final double getAverageAnnotation(Map.Entry<String, String>[] myVariantsAndAnno, String annoName) {
        return -1;
        //   return myGA.getAverageAnnotation(annoName);
    }
}
