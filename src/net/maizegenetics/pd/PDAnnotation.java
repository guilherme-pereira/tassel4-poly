package net.maizegenetics.pd;

import ch.systemsx.cisd.hdf5.*;
import net.maizegenetics.pal.alignment.ImportUtils;

import java.io.*;
import java.util.*;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.HapMapHDF5Constants;
import net.maizegenetics.util.Utils;

public class PDAnnotation {

    private static final String HAS_DATA = "HasData"; // summary index where if any trait has a value at that location, value is set to 1
    private boolean myIsSBit = true;

    private String hapMapFile_prefix = "/maizeHapMapV2_B73RefGenV2_201203028_";
    private String hapMapFile_suffix = ".hmp.txt.gz";
//    private String hapMapFile_suffix = "h.hmp.txt.gz";    // truncated hapmap file for  chr 9
    private static final int PHYSICAL_POSITION_COLUMN = 0;
    private static final int MINOR_ALLELE_FREQUENCY_COLUMN = 1;
    private static final int COLUMN_OFFSET = 1; // one for physical position column and another for minor allele frequency column

    private static final int TRAIT_INDEX = 0;
    private static final int CHR_INDEX = 1;
    private static final int PHYS_POS_INDEX = 2;
    private static final int RESULT_INDEX = 5;

    private String[] allTraits;
    private int[][] allPositions;
    private float[][] allResults;
    private int[][] featurePositions;
    private String[] allFeatures;



    public PDAnnotation(String hapMapPath, String pathGwas, String annoPath, String outputFile,
            int startChr, int endChr) {
        File aHapMapDir = new File(hapMapPath);
        File aGwasDir = new File(pathGwas);

        // Ed
//        loadGWAS(aGwasDir, outputFile, startChr);        // previous version - for comparison/testing

        // Dallas
        createAnnoHDF5WithHapMap(aHapMapDir, outputFile, startChr, endChr);
        loadGWAS(aGwasDir, outputFile, "\t", startChr, endChr);  // new version
//        File annoFile = new File(annoPath);
        //instatiate annotations in the HDF5 file - GENE (String), DistToGene (Integer), 
        //downstream_gene_variant, 3_prime_UTR_variant, missense_variant, synonymous_variant
//        loadAnnotationsByFile(annoFile, "\t");
    }



    public void createAnnoHDF5WithHapMap(File hapMapDir, String outputFile, int startChr, int endChr) {
        IHDF5WriterConfigurator config = HDF5Factory.configure(outputFile);
        config.overwrite();
        IHDF5Writer writer = config.writer();
        int[] hasData = null;  //recorder if there is any GWAS data for a site; TODO perhaps increment
        for (int currChr = startChr; currChr <= endChr; currChr++) {
            String chromosomeFile = hapMapDir + hapMapFile_prefix + "chr" + currChr + hapMapFile_suffix;
            System.out.println("Loading:" + chromosomeFile);
            Alignment bna = ImportUtils.readFromHapmap(chromosomeFile, myIsSBit, null /*progressListener*/);
            //System.out.printf("Sites:%d StartPosition:%d EndPosition:%d %n", bna.getSiteCount(), bna.getPositionInLocus(0), bna.getPositionInLocus(bna.getSiteCount() - 1));
            bna.optimizeForSites(null);
            int siteCnt = bna.getSiteCount();
            int[] alignmentPhysPos = bna.getPhysicalPositions();
            hasData = new int[siteCnt];
            float[] maf = new float[siteCnt];
            byte[] mjAllele = new byte[siteCnt];
            byte[] mnAllele = new byte[siteCnt];
            for (int j = 0; j < siteCnt; j++) {
                mjAllele[j] = bna.getMajorAlleleAsString(j).getBytes()[0];
                mnAllele[j] = bna.getMinorAlleleAsString(j).getBytes()[0];
                maf[j] = (float) bna.getMinorAlleleFrequency(j);
            }

            //write positions to hdf "pos"+chromosome
            String chrGroup = "chr" + currChr + "/";
            writer.createGroup(chrGroup);
            writer.setIntAttribute(chrGroup, HapMapHDF5Constants.NUM_SITES, siteCnt);
            HDF5IntStorageFeatures features= HDF5IntStorageFeatures.createDeflation(2);
            writer.createIntArray(chrGroup + HapMapHDF5Constants.POSITIONS, alignmentPhysPos.length, features);
            Arrays.sort(alignmentPhysPos);
            writer.writeIntArray(chrGroup + HapMapHDF5Constants.POSITIONS, alignmentPhysPos);

            //write alleles to hdf "allele"+chromosome
            // which version? String[][] ?
            writer.createByteArray(chrGroup + HapMapHDF5Constants.MAJOR_ALLELE, mjAllele.length);
            writer.writeByteArray(chrGroup + HapMapHDF5Constants.MAJOR_ALLELE, mjAllele);
            writer.createByteArray(chrGroup + HapMapHDF5Constants.MINOR_ALLELE, mnAllele.length);
            writer.writeByteArray(chrGroup + HapMapHDF5Constants.MINOR_ALLELE, mnAllele);
            // write minor allele frequencies
            HDF5FloatStorageFeatures floatFeatures= HDF5FloatStorageFeatures.createDeflateAndFloatScaling(2,3);
            writer.createFloatArray(chrGroup + HapMapHDF5Constants.MAF, maf.length);
            writer.writeFloatArray(chrGroup + HapMapHDF5Constants.MAF, maf);

            writer.createGroup(chrGroup + HapMapHDF5Constants.GWAS);
            writer.createGroup(chrGroup + HapMapHDF5Constants.GENOMIC);
            writer.createGroup(chrGroup + HapMapHDF5Constants.POP_GEN);
        }

        writer.close();
    }

    /**
     * For the provided chromosome, writes out all positions and gwas results on a trait-by-trait basis
     *
     * @param gwasFileIn
     * @param outputFile
     * @param delimiter
     * @param startChr
     * @param endChr
     */
    private void loadGWAS(File gwasFileIn, String outputFile,  String delimiter, int startChr, int endChr){
        IHDF5Writer writer = HDF5Factory.open(outputFile);

        //1. ArrayList<String> traitsInFile=getTraitListFromGWASInputFile();
        //2. Evaluate whether the traits already exist in the HDF5 file, if not create
        //3. Add GWAS results to the HDF5 file

        for(int currChr = startChr; currChr<= endChr; currChr++){
        System.out.printf("Loading GWAS by chromosome:%d %n", currChr);
        loadDataByChromosome(gwasFileIn, currChr,  delimiter);
        String chrGroup = "chr" + currChr + "/";

        int[] positions = writer.readIntArray(chrGroup + HapMapHDF5Constants.POSITIONS);
        Arrays.sort(positions);

        try{
        write(positions);
        }catch(IOException ioe){
            ioe.printStackTrace();
        }

            for (int i = 0; i < allTraits.length; i++) {
                int posMatch = 0, posMisMatch = 0;
                float[] rmip = new float[positions.length];
                for (int j = 0; j < allPositions[i].length; j++) {
//                if(allPositions[i][j]>3600000) continue;  // TODO: remove after testing

                    // for the current traits, transfer result values to array
                    int[] aInt = allPositions[i];
                    int site = Arrays.binarySearch(positions, allPositions[i][j]);
                    if (site < 0) {
                        System.out.println("Chr: \t" + currChr + " \tTrait: \t" + allTraits[i] + " \tPosition not found: \t" + allPositions[i][j] );
                        posMisMatch++;
                    } else {
                        posMatch++;
                        rmip[site] = allResults[i][j];
//                    System.out.printf("Hit Chr:%d Trait:%s Position:%d site:%d rmip:%f %n ", currChr, allTraits[i], allPositions[i][j], site, allResults[i][j]);
                    }
                }
                System.out.printf("Chr: %d Trait: %s Position matches:%d errors:%d %n", currChr, allTraits[i], posMatch, posMisMatch);
                String dataSetName = chrGroup + HapMapHDF5Constants.GWAS + "/" + allTraits[i];
                writer.createFloatArray(dataSetName, rmip.length);
                writer.writeFloatArray(dataSetName, rmip);
            }
        }
        writer.close();
    }

    public  void write ( int[]x) throws IOException{
        BufferedWriter outputWriter = null;
        outputWriter = new BufferedWriter(new FileWriter("/home/local/CORNELL/dek29/Documents/BucklerLab/PD/gwas/test/positions.txt"));
        for (int i = 0; i < x.length; i++) {
            // Maybe:
            outputWriter.write(x[i]+"\n");
            // Or:
//            outputWriter.write(Integer.toString(x[i]));
            outputWriter.newLine();
        }
        outputWriter.flush();
        outputWriter.close();
    }

    // Original - deprecated
    private void loadGWAS(File gwasFileIn, String outputFile, int currChr) {
        IHDF5Writer writer = HDF5Factory.open(outputFile);

        String[] traits =getGWASTraits(gwasFileIn, TRAIT_INDEX, "\t");

        String chrGroup = "chr" + currChr + "/";
        //read in all chromosome position
        //create a method to hold this memory


        int[] positions = writer.readIntArray(chrGroup + HapMapHDF5Constants.POSITIONS);

        for (int j = 0; j < traits.length; j++) {
           
            System.out.println(gwasFileIn.toString());
            // pull out the physical location and p-value
            int posMatch = 0, posMisMatch = 0;
            float[] rmip = new float[positions.length];
            Arrays.fill(rmip, Float.NaN);
            try {
                BufferedReader fileIn = Utils.getBufferedReader(gwasFileIn, 1000000);
                String s;
                while ((s = fileIn.readLine()) != null) {
                    String[] fields = s.split("\t");
                    try {
                        int theChr = Integer.parseInt(fields[CHR_INDEX]);
                        int position = Integer.parseInt(fields[PHYS_POS_INDEX]);
                        float rmipValue = Float.parseFloat(fields[RESULT_INDEX]);
                        if(theChr!=9) continue;
//                        if(position>7600000) continue;
                        //int site = Arrays.binarySearch(allPositions[theChr-1], position);
                        int site = Arrays.binarySearch(positions, position); 
                        if (site < 0) {
                            System.out.println("Error Position not found:"+position);
                            System.out.println("s = " + s);
                            posMisMatch++;
                        } else {
                            posMatch++;
                            rmip[site] = rmipValue;
//                            System.out.printf("Hit Chr:%d Trait:%s Position:%d site:%d rmip:%d %n ",theChr, traits[j], position, site, rmipValue);
                        }

                    } catch (Exception e) {
                        //                     System.out.println("Header");
                    }
                }
            } catch (IOException e) {
                System.out.println("IOError");
                e.printStackTrace();
            }
            System.out.printf("Position matches:%d errors:%d %n", posMatch, posMisMatch);
            String dataSetName = chrGroup + HapMapHDF5Constants.GWAS + "/" + traits[j];
            writer.createFloatArray(dataSetName, rmip.length);
            writer.writeFloatArray(dataSetName, rmip);
        } // end of traits loop
    }

    // Only used in deprecated version of loadGWAS
    private String[] getGWASTraits(File gwasResults, int traitIndex, String delimiter){
        BufferedReader br = Utils.getBufferedReader(gwasResults, 1000000);
        String line = null;
        Set<String> aSet = new HashSet();
        try{
            while((line =  br.readLine()) != null){
                String[] fields = line.split(delimiter);

                aSet.add(fields[traitIndex]);
            }
        }catch(IOException ioe){
            ioe.printStackTrace();
        }
    
        String[] result = new String[aSet.size()];
        aSet.toArray(result);
        return result;
    }

    /**
     * For a given chromosome, loads all positions and results for all traits
     *
     * @param gwasFileIn
     * @param currChr
     * @param delimiter
     */
    private void loadDataByChromosome(File gwasFileIn, int currChr, String delimiter){
        BufferedReader br = Utils.getBufferedReader(gwasFileIn, 1000000);
        String line = null;

        Map<String, List> traitPosition = new HashMap<String, List>();
        Map<String, List> traitResult = new HashMap<String, List>();
        try{
            while(( line = br.readLine()) != null) {
                String[] fields = line.split(delimiter);

                // for the current chromosome and trait, hold the positions
                try {
                    int chromosome = Integer.parseInt(fields[CHR_INDEX]);
                    if (currChr != chromosome) continue;

                    String aTrait = fields[TRAIT_INDEX].trim();
                    if (traitPosition.containsKey(aTrait)) {
                        traitPosition.get(aTrait).add(fields[PHYS_POS_INDEX]);
                        List l = traitResult.get(aTrait);
//                        System.out.println( "trait:" + aTrait + "loading result: " + fields[RESULT_INDEX]);
                        l.add(fields[RESULT_INDEX]);
                        traitResult.put(aTrait, l);
                    } else {
                        List<String> position = new ArrayList();
                        position.add(fields[PHYS_POS_INDEX]);
                        traitPosition.put(aTrait, position);
                        List<String> result = new ArrayList();
                        result.add(fields[RESULT_INDEX]);
                        traitResult.put(fields[TRAIT_INDEX], result);
                    }
                } catch (Exception e) {
//                    System.out.println("Header");
                }
            }
        }catch(IOException ioe){
            ioe.printStackTrace();
        }

        allTraits = new String[traitPosition.size()];

        // create a two-dimensional array of positions for each trait
        // first dimension index is shared with allTraits
        traitPosition.keySet().toArray(allTraits);

        int traitCount = allTraits.length;
        allPositions = new int[traitCount][];
        allResults = new float[traitCount][];
        for(int i = 0; i < traitCount; i++){
            List posList = traitPosition.get(allTraits[i]);
            int posCount = posList.size();
            allPositions[i] = new int[posCount];
            Iterator posIt = posList.iterator();
            int count = 0;
            while(posIt.hasNext()){
                allPositions[i][count++] = Integer.parseInt((String)posIt.next());
            }

            List resList = traitResult.get(allTraits[i]);
            int resCount = resList.size();
            allResults[i] = new float[resCount];
            Iterator resIt = resList.iterator();
            count = 0;
            while(resIt.hasNext()){
                allResults[i][count++] = Float.parseFloat((String)resIt.next());
            }
        }
    }

    // annotations files have been organized by chromosome
    private void loadAnnotationsByFile(File annoFile, String delimiter){

        int snpIdIndex = 0;
        int locationIndex = 1;      // location is  specified as <chr>:<position>, e.g., 9:30
                                    // TODO: how to handle range locations? e.g., 9:513883-513884
        String locationDelimiter = ":";
        int featureIndex = 6;     // Feature_typeConsequence

        Map<String, List> featurePosition = new HashMap<String, List>();

        BufferedReader br = Utils.getBufferedReader(annoFile, 1000000);
        String line = null;
        try {
            while ((line = br.readLine()) != null) {
                String[] fields = line.split(delimiter);

                try {
                    String aLoc =  fields[locationIndex];
                    int aPosition = getPosition(aLoc, locationDelimiter);
                    String feature = fields[featureIndex].trim();
                    if (featurePosition.containsKey(feature)) {
                        List l = featurePosition.get(feature);
                        l.add(aPosition);
                        featurePosition.put(feature, l);
                    } else {
                        List<Integer> l = new ArrayList<Integer>();
                        l.add(aPosition);
                        featurePosition.put(feature, l);
                    }

                } catch (Exception e) {
                    System.out.println("Header");
                }
            }
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

        // convert to a two dimensional array
        int featureCount = featurePosition.size();
        featurePositions = new int[featureCount][];
        allFeatures = new String[featureCount];
        featurePosition.keySet().toArray(allFeatures);
        for(int i = 0; i < featureCount; i++){
            List posList = featurePosition.get(allFeatures[i]);
            int posCount = posList.size();
            featurePositions[i] = new int[posCount];
            Iterator<Integer> iterator = posList.iterator();
            for(int j = 0; j < posCount; j++){
                featurePositions[i][j] = iterator.next().intValue();
            }
        }
    }


    private int getPosition(String location, String delimiter){
        String[] fields = location.split(delimiter);
        int position = Integer.parseInt(fields[1]);
        return position;
    }


    public static void main(String[] args) {
        // Dallas
        String hapMapPath = "/home/local/CORNELL/dek29/Documents/BucklerLab/PD/HapMap/compressed";
        String pathGwas = "/home/local/CORNELL/dek29/Documents/BucklerLab/PD/gwas/gwas_hits_all_sorted.txt";
        String PDfilePrefix = "/home/local/CORNELL/dek29/Documents/BucklerLab/PD/out/20130913_chr";
        String PDfileSuffix = ".h5";
        String pdFile = "/home/local/CORNELL/dek29/Documents/BucklerLab/PD/out/20130923_pd.h5";
        String annoPath = "/home/local/CORNELL/dek29/Documents/BucklerLab/PD/Annotations/20130522_SnpAnnotations_FromJason/maizeHapMapV2_B73RefGenV2_201203028_chr9h.WorstPerSnp.vcf";


        // Ed
//         String hapMapPath = "/Volumes/LaCie/HapMapV2/compressed/";
//        String pathGwas = "/Volumes/LaCie/PolymorphismDescriptors/gwas_hits_all.txt";
//        String PDfile = "/Volumes/LaCie/PolymorphismDescriptors/XtestPD.h5";
//        String annoPath = "/Volumes/LaCie/PolymorphismDescriptors/maizeHapMapV2_B73RefGenV2_201203028_chr9h.WorstPerSnp.vcf";

//        for(int i = 1; i < 11; i++){
//
//            String pdFile = null;
//            if(i<10){
//                pdFile = PDfilePrefix + "0" + i + PDfileSuffix;
//            }else{
//                pdFile = PDfilePrefix + i + PDfileSuffix;
//            }
            new PDAnnotation(hapMapPath, pathGwas, annoPath, pdFile, 1, 10);
//        }


    }
}
