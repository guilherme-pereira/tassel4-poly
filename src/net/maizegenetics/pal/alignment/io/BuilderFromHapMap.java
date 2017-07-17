package net.maizegenetics.pal.alignment.io;

import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.pal.alignment.genotype.Genotype;
import net.maizegenetics.pal.alignment.genotype.GenotypeBuilder;
import net.maizegenetics.pal.site.*;
import net.maizegenetics.pal.taxa.AnnotatedTaxon;
import net.maizegenetics.pal.taxa.TaxaList;
import net.maizegenetics.pal.taxa.TaxaListBuilder;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;

/**
 * Create an alignment based on HapMap format file (either .txt or compressed).  Alleles are set as global major and
 * global minor.
 * e.g. code <p></p>
 * {@code
 * AlignmentNew a=BuilderFromHapMap.getBuilder(infileName).build();
 * }
 * <p></p>
 * TODO:  Add filtering while reading, provide an option to define the alleles as reference and alternate
 *
 * @author Ed Buckler
 */
public class BuilderFromHapMap {
    private static final Logger myLogger = Logger.getLogger(BuilderFromHapMap.class);
    private static final Pattern WHITESPACE_PATTERN = Pattern.compile("\\s");
    private static final int NUM_HAPMAP_NON_TAXA_HEADERS = 11;
    private final String infile;

    private BuilderFromHapMap(String infile) {
        this.infile=infile;
    }

    public static BuilderFromHapMap getBuilder(String infile) {
        return new BuilderFromHapMap(infile);
    }

    //TODO provide options on caching to use, read only some sites, etc.

    public AlignmentNew build() {
        long time=System.nanoTime();
        AlignmentNew result=null;
        try {
            int numThreads = Runtime.getRuntime().availableProcessors();
            ExecutorService pool = Executors.newFixedThreadPool(numThreads);
            BufferedReader r=Utils.getBufferedReader(infile, -1);
            TaxaList taxaList=processTaxa(r.readLine());
            String currLine;
            int linesAtTime=1<<12;  //this is a critical lines with 20% or more swings.  Needs to be optimized with transposing
          //  int linesAtTime=1<<8;  //better for with lots of taxa.
            ArrayList<String> txtLines=new ArrayList<>(linesAtTime);
            ArrayList<ProcessHapMapBlock> pbs=new ArrayList<>();
            int lines=0;
            while((currLine=r.readLine())!=null) {
                txtLines.add(currLine);
                lines++;
                if(lines%linesAtTime==0) {
                    ProcessHapMapBlock pb=ProcessHapMapBlock.getInstance(pbs.size(), taxaList.getTaxaCount(), txtLines);
                    pbs.add(pb);
               //     pb.run();
                    pool.execute(pb);
                    txtLines=new ArrayList<>(linesAtTime);
                }
            }
            r.close();
            if(txtLines.size()>0) {
                ProcessHapMapBlock pb=ProcessHapMapBlock.getInstance(pbs.size(), taxaList.getTaxaCount(), txtLines);
                pbs.add(pb);
                pool.execute(pb);
            }
            pool.shutdown();
            if (!pool.awaitTermination(60, TimeUnit.SECONDS)) {
                throw new IllegalStateException("BuilderFromHapMap: processing threads timed out.");
            }
            int currentSite=0;
            PositionArrayList.Builder posBuild=new PositionArrayList.Builder();
            GenotypeBuilder gb = GenotypeBuilder.getUnphasedNucleotideGenotypeBuilder(taxaList.getTaxaCount(), lines);
            for(ProcessHapMapBlock pb: pbs) {
                posBuild.addAll(pb.getBlkPosList());
                byte[][] bgTS=pb.getGenoTS();
                for (int t=0; t<bgTS.length; t++) {
                    gb.setBaseRangeForTaxon(t,currentSite,bgTS[t]);
                }
                currentSite+=pb.getSiteNumber();
            }
            if(posBuild.validateOrdering()==false) {
                throw new IllegalStateException("BuilderFromHapMap: Ordering incorrect HapMap must be ordered by position");
            }
            Genotype g = gb.build();
            result=new CoreAlignment(g,posBuild.build(),taxaList,null,null);
        } catch (IOException|InterruptedException e) {
            e.printStackTrace();
        }
        long totalTime=System.nanoTime()-time;
        System.out.printf("BuilderFromHapMap data timing %gs %n", totalTime/1e9);
        return result;
    }

    private static TaxaList processTaxa(String readLn) {
        String[] header = WHITESPACE_PATTERN.split(readLn);
        int numTaxa = header.length - NUM_HAPMAP_NON_TAXA_HEADERS;
        TaxaListBuilder tlb = new TaxaListBuilder();
        for (int i = 0; i < numTaxa; i++) {
            AnnotatedTaxon at = new AnnotatedTaxon.Builder(header[i + NUM_HAPMAP_NON_TAXA_HEADERS])
                    .build();
            tlb.add(at);
        }
        return tlb.build();
    }

}

class ProcessHapMapBlock implements Runnable {
    private static final Pattern WHITESPACE_PATTERN = Pattern.compile("\\s");
    private static final Pattern SLASH_PATTERN = Pattern.compile("/");
    private static final int NUM_HAPMAP_NON_TAXA_HEADERS = 11;
    private static final int GENOIDX = NUM_HAPMAP_NON_TAXA_HEADERS;
    private static final int SNPID_INDEX= 0;
    private static final int VARIANT_INDEX= 1;
    private static final int CHROMOSOME_INDEX = 2;
    private static final int POSITION_INDEX= 3;
    private final int order;
    private final int taxaN;
    private final int siteN;
    private ArrayList<String> txtL;
    private byte[][] gTS;
    private final ArrayList<Position> blkPosList;
    private final byte[] convert;
    private final boolean isOneLetter; //true e.g. A,R, false=AA,CT

    private ProcessHapMapBlock(int order, int taxaN, ArrayList<String> txtL) {
        this.order=order;
        this.taxaN=taxaN;
        this.siteN=txtL.size();
        this.txtL=txtL;
        blkPosList=new ArrayList<>(siteN);
        String[] tokens = WHITESPACE_PATTERN.split(txtL.get(0),NUM_HAPMAP_NON_TAXA_HEADERS+1);
        double avg=(double)(tokens[GENOIDX].length()+1)/(double)taxaN;
        if((avg>1.99)&&(avg<2.01)) {isOneLetter=true;}
        else if((avg>2.99)&&(avg<3.01)) {isOneLetter=false;}
        else{throw new IllegalStateException("ProcessHapMapBlock: Genotype coded wrong use one or 2 letters per genotype");}

        //todo Move to Nucleotide constants if best way
        convert = new byte[128];
        for (int i = 0; i < convert.length; i++) {
            try {
                convert[i] = NucleotideAlignmentConstants.getNucleotideDiploidByte((char) i);
            } catch (IllegalArgumentException e) {
                convert[i] = Alignment.UNKNOWN_DIPLOID_ALLELE;
            }
        }
    }

    public static ProcessHapMapBlock getInstance(int order, int taxaN, ArrayList<String> txtL) {
        return new ProcessHapMapBlock(order,taxaN, txtL);
    }


    //@Override
    public void run() {
        Map<String, Chromosome> chromosomeLookup = new HashMap<>();
        gTS=new byte[taxaN][siteN];
        for (int s=0; s<siteN; s++) {
            String input=txtL.get(s);
            int[] tabPos=new int[NUM_HAPMAP_NON_TAXA_HEADERS];
            int tabIndex=0;
            int len=input.length();
            for (int i = 0; (tabIndex<NUM_HAPMAP_NON_TAXA_HEADERS) && (i<len); i++) {
                if(input.charAt(i)=='\t') tabPos[tabIndex++]=i;
            }
            String chrName=input.substring(tabPos[CHROMOSOME_INDEX-1]+1,tabPos[CHROMOSOME_INDEX]);
            Chromosome currChr = chromosomeLookup.get(chrName);
            if (currChr == null) {
                currChr = new Chromosome(new String(chrName));
                chromosomeLookup.put(chrName, currChr);
            }
            String variants=input.substring(tabPos[VARIANT_INDEX-1]+1,tabPos[VARIANT_INDEX]);
            GeneralPosition.Builder apb= new GeneralPosition.Builder(currChr,
                    Integer.parseInt(input.substring(tabPos[POSITION_INDEX-1]+1,tabPos[POSITION_INDEX])))
                    .snpName(input.substring(0,tabPos[SNPID_INDEX]))
                    .knownVariants(variants)
                            //TODO                    strand, variants,
                    ;
            try{
                byte glbMajor=convert[variants.charAt(0)];
                apb.allele(Position.Allele.GLBMAJ,glbMajor);
                if(variants.length()==3) {
                    byte glbMinor=convert[variants.charAt(2)];
                    apb.allele(Position.Allele.GLBMIN,glbMinor);
                }
            } catch (IllegalArgumentException e) {
                //for the indels that cannot be converted correctly now
                // System.out.println("Error Parsing this variant"+Arrays.toString(variants));
            }
            blkPosList.add(apb.build());
            int offset=tabPos[NUM_HAPMAP_NON_TAXA_HEADERS-1]+1;
            if(isOneLetter) {
                for (int i = offset; i < len; i+=2) {
                    gTS[(i-offset)/2][s] = convert[input.charAt(i)];
                }
            } else {
                for(int i = offset; i < len; i+=3) {
                    //System.out.println(i+":"+input.charAt(i+1)+input.charAt(i));
                    //there is a phasing conflict with the existing import approach
                    gTS[(i-offset)/3][s]=AlignmentUtils.getDiploidValue(convert[input.charAt(i+1)], convert[input.charAt(i)]);
                }
            }
        }
        txtL=null;
    }


    int getSiteNumber() {
        return siteN;
    }

    byte[][] getGenoTS() {
        return gTS;
    }

    ArrayList<Position> getBlkPosList() {
        return blkPosList;
    }
}
