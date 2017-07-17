package net.maizegenetics.pd;

import ch.systemsx.cisd.hdf5.*;

import java.util.ArrayList;
import java.util.List;

/**
 * User: dkroon
 * Date: 6/25/13
 */
public class PDReport {

    private String polyDescFile;
    public PDReport(String polymorphismDescriptorFile){
        this.polyDescFile = polymorphismDescriptorFile;

        IHDF5Reader reader = HDF5Factory.openForReading(polymorphismDescriptorFile);
        //  int[] variableSites = reader.readIntArray(HapMapHDF5Constants.POSITIONS);
        String delimiter="\t";
        List<HDF5LinkInformation> fields=reader.getAllGroupMemberInformation("chr9", true);

        List<HDF5LinkInformation> fields2=new ArrayList(fields);
        for (HDF5LinkInformation is : fields) {
            //if(is.isGroup()==false) continue;
            if(is.isGroup())fields2.addAll(reader.getAllGroupMemberInformation(is.getPath(), true));
        }
        float[][] fa=new float[20][];
        String[] fNames=new String[20];
        int[][] ia=new int[20][];
        String[] iNames=new String[20];
        int currentFA=0;
        int currentIA=0;
        for (HDF5LinkInformation is : fields2) {
            System.out.println(is.getPath().toString()+"::"+reader.getObjectType(is.getPath()).toString());
            if(is.isDataSet()==false) continue;
            HDF5DataSetInformation info=reader.getDataSetInformation(is.getPath());
            if(info.getTypeInformation().getDataClass()== HDF5DataClass.FLOAT) {
                fNames[currentFA]=is.getName();
                fa[currentFA]=reader.readFloatArray(is.getPath());
                System.out.println("currentFA = " + currentFA);
                currentFA++;
            } else if(info.getTypeInformation().getDataClass()==HDF5DataClass.INTEGER) {
                iNames[currentIA]=is.getName();
                ia[currentIA]=reader.readIntArray(is.getPath());
                currentIA++;
            }

            System.out.println(is.getPath().toString()+"::"+reader.getDataSetInformation(is.getPath()).toString());
        }
        StringBuilder sb=new StringBuilder("Site"+delimiter);
        for (int fi = 0; fi < currentFA; fi++) {sb.append(fNames[fi]); sb.append(delimiter);}
        for (int ii = 0; ii < currentIA; ii++) {sb.append(iNames[ii]); sb.append(delimiter);}
        System.out.println(sb.toString());
        for (int i = 0; i < fa[0].length; i++) {
            sb=new StringBuilder();

            for (int fi = 0; fi < currentFA; fi++) {
                sb.append(fa[fi][i]);sb.append(delimiter);
            }
            for (int ii = 0; ii < currentIA; ii++) {
                sb.append(ia[ii][i]);sb.append(delimiter);
            }
        }
        sb.toString();
    }




    // because everything is organized by chromosome, iterate through each chromosome unless locus is specified


    public void reportAllSitesGWASStatus(String outFile, double minAlleleFreq, double maxAlleleFreq, char delimiter){

    }


    public static void main(String[] args) {


        String PDfile = "C:\\Documents and Settings\\dkroon\\My Documents\\PD\\out\\anno_testPD.h5";
        new PDReport(PDfile);
    }
}
