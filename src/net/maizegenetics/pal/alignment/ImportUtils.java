/*
 * ImportUtils
 */
package net.maizegenetics.pal.alignment;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.*;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 * The class imports Alignment from various file formats.
 *
 * @author terry, Gabriel Rodrigues Alves Margarido
 */
public class ImportUtils {

    private static final Logger myLogger = Logger.getLogger(ImportUtils.class);
    private static final Pattern WHITESPACE_PATTERN = Pattern.compile("\\s");
    public static final int NUM_HAPMAP_NON_TAXA_HEADERS = 11;
    public static final int HAPMAP_SNPID_COLUMN_INDEX = 0;
    public static final int HAPMAP_CHROMOSOME_COLUMN_INDEX = 2;
    public static final int HAPMAP_POSITION_COLUMN_INDEX = 3;
    public static final int NUM_VCF_NON_TAXA_COLUMNS = 9;
    public static final int VCF_CHROMOSOME_COLUMN_INDEX = 0;
    public static final int VCF_POSITION_COLUMN_INDEX = 1;
    public static final int VCF_SNPID_COLUMN_INDEX = 2;
    public static final int VCF_REF_COLUMN_INDEX = 3;
    public static final int VCF_ALT_COLUMN_INDEX = 4;
    public static final int VCF_FORMAT_COLUMN_INDEX = 8;

    private ImportUtils() {
        // Utility Class - do not instantiate.
    }

    public static Alignment readGuessFormat(String fileName, boolean readSBit) {
        try {
            if (fileName.endsWith("hmp.h5") || fileName.endsWith("mhmp.h5")) {
                IHDF5Reader reader = HDF5Factory.openForReading(fileName);
                boolean geno = reader.exists(HapMapHDF5Constants.GENOTYPES);
                reader.close();
                if (geno) {
                    return MutableNucleotideAlignmentHDF5.getInstance(fileName);
                }
                return BitAlignmentHDF5.getInstance(fileName, readSBit);
            } else if (fileName.endsWith("hmp.txt.gz") || fileName.endsWith("hmp.txt")) {
                return readFromHapmap(fileName, readSBit, null);
            } else if (fileName.endsWith(".vcf") || fileName.endsWith(".vcf.gz")) {
                return readFromVCF(fileName, null);
            } else {
                return readFasta(fileName, readSBit);
            }
        } catch (Exception e) {
            System.err.println("Error reading:" + fileName);
            return null;
        }

    }

    /*
     * Counts number of Header rows in a VCF files
     * Use in conjunction with util.getNumberLines to count numSites
     */
    private static int getNumHeaderRowsVCF(String filename) {
        BufferedReader fileIn = null;
        try {
            int numHeader = 0;
            fileIn = Utils.getBufferedReader(filename, 1000000);

            String currLine = fileIn.readLine();
            while (currLine != null) {
                if (currLine.substring(0, 1).equals("#")) {
                    numHeader++;
                } else {
                    break;
                }
                currLine = fileIn.readLine();
            }

            return numHeader;
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Error in getNumSiteVCF, unable to read VCF file: " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                fileIn.close();
            } catch (Exception e) {
                // do nothing
            }
        }
    }

    // add support for SNP ID?
    public static Alignment readFromVCF(final String filename, ProgressListener listener, int maxKeptAlleles) {

        int minPosition = Integer.MAX_VALUE;
        String currLocus = null;

        Pattern colonPattern = Pattern.compile(":");
        Pattern commaPattern = Pattern.compile(",");

        long currentTime = System.currentTimeMillis();
        int numHeader = getNumHeaderRowsVCF(filename);
        int numSites = Utils.getNumberLines(filename) - numHeader;
        myLogger.info("readFromVCF: Number of Header Rows: " + numHeader);
        myLogger.info("readFromVCF: Number of Sites: " + numSites);

        long prevTime = currentTime;
        currentTime = System.currentTimeMillis();
        myLogger.info("readFromVCF: Time to count lines: " + ((currentTime - prevTime) / 1000));


        BufferedReader fileIn = null;
        try {
            fileIn = Utils.getBufferedReader(filename, 1000000);
            String currLine = null;
            for (int i = 0; i < numHeader; i++) {
                currLine = fileIn.readLine();
            }

            // get taxa names
            String[] header = WHITESPACE_PATTERN.split(currLine);
            int numTaxa = header.length - NUM_VCF_NON_TAXA_COLUMNS;
            String[] taxaNames = new String[numTaxa];
            System.arraycopy(header, NUM_VCF_NON_TAXA_COLUMNS, taxaNames, 0, numTaxa);
            IdGroup idGroup = new SimpleIdGroup(taxaNames);

            String[] snpID = new String[numSites];

            MutableVCFAlignment result = MutableVCFAlignment.getInstance(idGroup, numSites, numTaxa, numSites, maxKeptAlleles);


            int prevPos = -1;
            int currPos = -1;
            int locusStart = 0;

            for (int site = 0; site < numSites; site++) {
                currLine = fileIn.readLine();
                String[] currSite = WHITESPACE_PATTERN.split(currLine);
                // 1	6407	S1000_6407	.	A,C	20	PASS	NS=929;DP=2210;AF=0.01,0.01	GT:AD:DP:GQ:PL	./.	0/0:2,0,0:2:79:0,6,72	1/2:0,1,1:2:54:0,5,34
                String temp = currSite[VCF_CHROMOSOME_COLUMN_INDEX];
                currPos = Integer.parseInt(currSite[VCF_POSITION_COLUMN_INDEX]);
                snpID[site] = currSite[VCF_SNPID_COLUMN_INDEX];
                String ref = currSite[VCF_REF_COLUMN_INDEX];
                String alt = currSite[VCF_ALT_COLUMN_INDEX];
                String format = currSite[VCF_FORMAT_COLUMN_INDEX];
                String[] dataByTaxa = new String[numTaxa];
                System.arraycopy(currSite, NUM_VCF_NON_TAXA_COLUMNS, dataByTaxa, 0, numTaxa);

                //Currently this code can only support alleles that are encoded in single charactor, 
                //Indels need to be encoded as "+" or "-", multi-charactor indels are not supported 
                if (ref.length() > 1) {
                    throw new IllegalStateException("ImportUtils: readFromVCF: the reference allele must be in single charactor. Multi-character indels need to be converted to '+' or '-'. At position " + currPos + ", the ref allele '" + ref + "' is not supported.");
                }
                String[] altsplitted = alt.split(",");
                for (String t : altsplitted) {
                    if (t.length() != 1) {
                        throw new IllegalStateException("ImportUtils: readFromVCF: the alternative allele must be in single charactor. Multi-character indels need to be converted to '+' or '-'. At position " + currPos + ", the alternative allele '" + t + "' is not supported.");
                    }
                }
                alt = alt.replaceAll(",", "");

                // find alleles for current site, check to see if number of alleles is supported
                int numAlleles = alt.length() + 1;
                int oriNumAlleles = numAlleles;
                String alleleString = ref + alt;
                if (numAlleles > maxKeptAlleles) {
                    System.out.println("read VCF position " + currPos + ": extra allele(s) removed.");
                    numAlleles = maxKeptAlleles;
                    alleleString = alleleString.substring(0, numAlleles);
                    //throw new IllegalStateException("ImportUtils: readFromVCF: number of Alleles is larger than allowed currently in TASSEL: " + numAlleles + " alleles found, " + maxKeptAlleles + " alleles allowed, in line " + (numHeader + site + 1));
                }


                byte[] alleles = new byte[numAlleles];

                for (int allele = 0; allele < numAlleles; allele++) {
                    String currAllele = alleleString.substring(allele, allele + 1);
                    if (currAllele.equals(".")) {
                        alleles[allele] = (byte) -1;
                    } else if (currAllele.equalsIgnoreCase("A")) {
                        alleles[allele] = NucleotideAlignmentConstants.A_ALLELE;
                    } else if (currAllele.equalsIgnoreCase("C")) {
                        alleles[allele] = NucleotideAlignmentConstants.C_ALLELE;
                    } else if (currAllele.equalsIgnoreCase("G")) {
                        alleles[allele] = NucleotideAlignmentConstants.G_ALLELE;
                    } else if (currAllele.equalsIgnoreCase("T")) {
                        alleles[allele] = NucleotideAlignmentConstants.T_ALLELE;
                    } else if (currAllele.equals("+")) {
                        alleles[allele] = NucleotideAlignmentConstants.INSERT_ALLELE;
                    } else if (currAllele.equals("-")) {
                        alleles[allele] = NucleotideAlignmentConstants.GAP_ALLELE;
                    } else {
                        throw new IllegalStateException("ImportUtils: readFromVCF: a unsupported allele detected in this VCF file: " + currAllele + " in line " + (numHeader + site + 1));
                    }
                }


                result.setCommonAlleles(site, alleles);
                result.setPositionOfSite(site, currPos);
                if (!ref.equals(".")) {
                    result.setReferenceAllele(site, NucleotideAlignmentConstants.getNucleotideDiploidByte(ref));
                }

                // get the possible alleles for each site in to an byte array
                // result.setCommonAlleles(site, values);

                // figure out order of format
                int genoIndex = -1;
                int alleleDepthIndex = -1;

                String[] formatSplit = colonPattern.split(format);
                int numDataFields = formatSplit.length;
                for (int i = 0; i < formatSplit.length; i++) {
                    if (formatSplit[i].equalsIgnoreCase("GT")) {
                        genoIndex = i;
                    } else if (formatSplit[i].equalsIgnoreCase("AD")) {
                        alleleDepthIndex = i;
                    }
                }

                if (genoIndex == -1) {
                    throw new IllegalStateException("ImportUtils: readFromVCF: no genotype data found in this VCF file at line: " + (numHeader + site + 1));
                }

                //if (alleleDepthIndex == -1) {
                //throw new IllegalStateException("ImportUtils: readFromVCF: no allele depth data found in this VCF file at line: " + (numHeader + site + 1));
                //}


                for (int taxa = 0; taxa < numTaxa; taxa++) {
                    String[] dataSplit = colonPattern.split(dataByTaxa[taxa]);
                    if (dataSplit[0].startsWith("./.")) {
                        short[] depths = new short[numAlleles];
                        java.util.Arrays.fill(depths, (short) 0);
                        result.setDepthForAlleles(taxa, site, depths);
                        continue;
                    }
                    // for whatever reason if the actual data fields do not match up with the format column
                    // assume the data is unknown
                    byte calledGenotypeValue = (byte) 0xFF;
                    if (dataSplit.length != numDataFields) {
                        result.setBase(taxa, 0, calledGenotypeValue);
                    } else {
                        String[] genotypes = Pattern.compile("[/|]").split(dataSplit[genoIndex]);
                        // String genotypes = dataSplit[genoIndex].replaceAll("/", "").replaceAll("|", "");
                        if (genotypes.length > 2) {
                            throw new IllegalStateException("ImportUtils: readFromVCF: number of genotypes larger than supported by TASSEL at line: " + (numHeader + site + 1));
                        }
                        boolean recallgenotypeflag = false; //this flag  indicate whethere to recall genotypes, this happens when some alleles are removed
                        for (int i = 0; i < genotypes.length; i++) {
                            calledGenotypeValue <<= 4;
                            String currGenotype = genotypes[i];
                            if (currGenotype.equals(".")) {
                                calledGenotypeValue |= 0x0F;
                            } else {
                                int currGenoInt = Integer.parseInt(currGenotype);
                                if (currGenoInt == 14) {
                                    calledGenotypeValue |= 0x0E;
                                } else {
                                    if (currGenoInt > (numAlleles - 1)) {
                                        recallgenotypeflag = true;
                                        break;
                                    } else {
                                        currGenotype = alleleString.substring(currGenoInt, currGenoInt + 1);
                                        if (currGenotype.equalsIgnoreCase("A")) {
                                            calledGenotypeValue |= 0x00;
                                        } else if (currGenotype.equalsIgnoreCase("C")) {
                                            calledGenotypeValue |= 0x01;
                                        } else if (currGenotype.equalsIgnoreCase("G")) {
                                            calledGenotypeValue |= 0x02;
                                        } else if (currGenotype.equalsIgnoreCase("T")) {
                                            calledGenotypeValue |= 0x03;
                                        } else if (currGenotype.equals("+")) {
                                            calledGenotypeValue |= 0x04;
                                        } else if (currGenotype.equals("-")) {
                                            calledGenotypeValue |= 0x05;
                                        } else {
                                            throw new IllegalStateException("ImportUtils: readFromVCF: a unsupported allele detected in this VCF file: " + currGenotype + " in line " + (numHeader + site + 1));
                                        }
                                    }
                                }
                            }
                        }

                        //the called genotype will be set later, incase the genotype needs to be recalled
                        if (alleleDepthIndex >= 0) {
                            String alleleDepths = dataSplit[alleleDepthIndex];
                            String[] stringDepths = commaPattern.split(alleleDepths);
                            if (stringDepths.length != oriNumAlleles) {
                                throw new IllegalStateException("ImportUtils: readFromVCF: number of allele depth values does not match number of alleles in line: " + (numHeader + site + 1) + " taxa number: " + taxa);
                            }

                            short[] depths = new short[numAlleles];
                            int[] intDepths = new int[numAlleles];
                            for (int i = 0; i < numAlleles; i++) {
                                int depth = Integer.parseInt(stringDepths[i]);
                                intDepths[i] = depth;
                                if (depth > 32767) {
                                    myLogger.info("Depth value for genotype " + i + " had an original value of " + depth + ". Converted to the maximum of 32767. In line: " + (numHeader + site + 1) + " taxa number: " + taxa);
                                    depth = 32767;
                                }
                                depths[i] = (short) depth;
                            }
                            if (recallgenotypeflag == true) {
                                //recall genotype now
                                calledGenotypeValue = VCFUtil.resolveVCFGeno(alleles, intDepths);
                            }
                            result.setDepthForAlleles(taxa, site, depths);
                        }

                        result.setBase(taxa, site, calledGenotypeValue);

                    }
                }

                if (currLocus == null) {
                    currLocus = temp;
                    minPosition = currPos;
                } else if (!temp.equals(currLocus)) {
                    Locus newLocus = new Locus(currLocus, currLocus, minPosition, prevPos, null, null);

                    for (int i = locusStart; i < site; i++) {
                        result.setLocusOfSite(i, newLocus);
                    }
                    currLocus = temp;
                    minPosition = currPos;
                    locusStart = site;
                    prevPos = -1;
                }

                if (currPos < prevPos) {
                    throw new IllegalStateException("ImportUtils: readFromVCF: Sites are not properly sorted for chromosome: " + currLocus + " at " + currPos + " and " + prevPos);
                }

                prevPos = currPos;
            }

            if (currLocus != null) {
                Locus newLocus = new Locus(currLocus, currLocus, minPosition, prevPos, null, null);

                for (int i = locusStart; i < numSites; i++) {
                    result.setLocusOfSite(i, newLocus);
                }
            }
            result.clean();
            return result;
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("ImportUtils: readFromVCF: Problem creating Alignment: " + filename + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                fileIn.close();
            } catch (Exception ex) {
                // do nothing
            }
        }
    }

    public static Alignment readFromVCF(final String filename, ProgressListener listener) {
        return readFromVCF(filename, listener, VCFUtil.VCF_DEFAULT_MAX_NUM_ALLELES);
    }

    public static Alignment readFromHapmap(final String filename, ProgressListener listener) {
        return readFromHapmap(filename, true, listener);
    }

    public static Alignment readFromHapmap(final String filename, boolean isSBit, ProgressListener listener) {

        int minPosition = Integer.MAX_VALUE;
        String currLocus = null;
        List<Locus> loci = new ArrayList<Locus>();
        List<Integer> lociOffsets = new ArrayList<Integer>();

        long currentTime = System.currentTimeMillis();
        int numSites = Utils.getNumberLines(filename) - 1;
        myLogger.info("readFromHapmap: Number of Sites: " + numSites);

        long prevTime = currentTime;
        currentTime = System.currentTimeMillis();
        myLogger.info("readFromHapmap: Time to count lines: " + ((currentTime - prevTime) / 1000));


        BufferedReader fileIn = null;
        try {
            int numThreads = Runtime.getRuntime().availableProcessors();
            ExecutorService pool = Executors.newFixedThreadPool(numThreads);

            fileIn = Utils.getBufferedReader(filename, 1000000);
            String[] header = WHITESPACE_PATTERN.split(fileIn.readLine());
            int lineInFile = 1;
            int numTaxa = header.length - NUM_HAPMAP_NON_TAXA_HEADERS;
            String[] snpIDs = new String[numSites];
            int prevPosition = -1;

            OpenBitSet[][] theData;
            byte[][] alleles = new byte[numSites][TasselPrefs.getAlignmentMaxAllelesToRetain()];
            int numDataRows = TasselPrefs.getAlignmentMaxAllelesToRetain();
            if (TasselPrefs.getAlignmentRetainRareAlleles()) {
                numDataRows++;
            }
            int numSitesToProcess = 1;
            if (isSBit) {
                theData = new OpenBitSet[numDataRows][numSites];
                numSitesToProcess = 1;
            } else {
                theData = new OpenBitSet[numDataRows][numTaxa];
                for (int al = 0; al < numDataRows; al++) {
                    for (int t = 0; t < numTaxa; t++) {
                        theData[al][t] = new OpenBitSet(numSites);
                    }
                }
                numSitesToProcess = 64;
            }

            int[] physicalPositions = new int[numSites];
            int count = 0;
            String[][] tokens = new String[numSitesToProcess][];
            int currentSite = 0;
            for (int site = 0; site < numSites; site++) {

                lineInFile++;

                String input = fileIn.readLine();
                tokens[count] = WHITESPACE_PATTERN.split(input);

                snpIDs[site] = new String(tokens[count][HAPMAP_SNPID_COLUMN_INDEX]);
                int position = Integer.parseInt(tokens[count][HAPMAP_POSITION_COLUMN_INDEX]);
                String temp = new String(tokens[count][HAPMAP_CHROMOSOME_COLUMN_INDEX]);
                if (currLocus == null) {
                    lociOffsets.add(site);
                    currLocus = temp;
                    minPosition = position;
                    prevPosition = -1;
                } else if (!temp.equals(currLocus)) {
                    loci.add(new Locus(currLocus, currLocus, minPosition, prevPosition, null, null));
                    lociOffsets.add(site);
                    currLocus = temp;
                    minPosition = position;
                    prevPosition = -1;
                }

                if (position < prevPosition) {
                    throw new IllegalStateException("ImportUtils: readFromHapmap: Sites are not properly sorted for chromosome: " + currLocus + " at " + position + " and " + prevPosition);
                }

                count++;

                if (count == numSitesToProcess) {
                    pool.execute(ProcessLineOfHapmap.getInstance(alleles, theData, TasselPrefs.getAlignmentRetainRareAlleles(), tokens, count, currentSite, numTaxa, lineInFile, isSBit));
                    count = 0;
                    currentSite += numSitesToProcess;
                    tokens = new String[numSitesToProcess][];
                }

                physicalPositions[site] = position;
                prevPosition = position;

                if (listener != null) {
                    listener.progress((int) (((double) (site + 1) / (double) numSites) * 100.0), null);
                }
            }

            if (count != 0) {
                pool.execute(ProcessLineOfHapmap.getInstance(alleles, theData, TasselPrefs.getAlignmentRetainRareAlleles(), tokens, count, currentSite, numTaxa, lineInFile, isSBit));
            }


            pool.shutdown();
            if (!pool.awaitTermination(6000, TimeUnit.SECONDS)) {
                throw new IllegalStateException("ImportUtils: readFromHapmap: processing threads timed out.");
            }

            if (currLocus != null) {
                loci.add(new Locus(currLocus, currLocus, minPosition, prevPosition, null, null));
            }

            prevTime = currentTime;
            currentTime = System.currentTimeMillis();
            myLogger.info("readFromHapmap: Time to read file: " + ((currentTime - prevTime) / 1000));

            String[] taxaNames = new String[numTaxa];
            System.arraycopy(header, NUM_HAPMAP_NON_TAXA_HEADERS, taxaNames, 0, numTaxa);
            IdGroup idGroup = new SimpleIdGroup(taxaNames);

            Locus[] lociFinal = new Locus[loci.size()];
            loci.toArray(lociFinal);
            int[] offsetsFinal = new int[lociOffsets.size()];
            for (int i = 0; i < lociOffsets.size(); i++) {
                offsetsFinal[i] = ((Integer) lociOffsets.get(i)).intValue();
            }

            Alignment result = BitAlignment.getNucleotideInstance(idGroup, alleles, theData, null, null, physicalPositions, TasselPrefs.getAlignmentMaxAllelesToRetain(), lociFinal, offsetsFinal, snpIDs, TasselPrefs.getAlignmentRetainRareAlleles(), isSBit);

            prevTime = currentTime;
            currentTime = System.currentTimeMillis();
            myLogger.info("readFromHapmap: Time to create Alignment: " + ((currentTime - prevTime) / 1000));

            return result;
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("ImportUtils: readFromHapmap: Problem creating Alignment: " + filename + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                fileIn.close();
            } catch (Exception ex) {
                // do nothing
            }
        }

    }

    public static Alignment readFasta(String filename, boolean isSBit) throws FileNotFoundException, IOException {

        BufferedReader reader = Utils.getBufferedReader(filename);

        List taxa = new ArrayList();
        List sequences = new ArrayList();

        String line = null;
        line = reader.readLine();
        boolean sequence = false;
        int sequenceLength = -1;
        int count = 1;
        while (line != null) {

            line = line.trim();

            if (line.startsWith(";")) {
                line = reader.readLine();
            } else if (line.startsWith(">")) {
                StringTokenizer tokens = new StringTokenizer(line);
                String taxaName = tokens.nextToken();
                if (taxaName.length() == 1) {
                    taxaName = tokens.nextToken();
                } else {
                    taxaName = taxaName.substring(1).trim();
                }
                taxa.add(taxaName);
                sequence = true;
                line = reader.readLine();
            } else if (sequence) {
                StringBuilder builder = new StringBuilder();
                while ((line != null) && (!line.startsWith(">")) && (!line.startsWith(";"))) {
                    line = line.trim().toUpperCase();
                    builder.append(line);
                    line = reader.readLine();
                }
                String temp = builder.toString();
                if (sequenceLength == -1) {
                    sequenceLength = temp.length();
                } else if (sequenceLength != temp.length()) {
                    throw new IllegalStateException("ImportUtils: readFasta: Sequence: " + count + " Differs in Length.");
                }
                sequences.add(temp);
                sequence = false;
                count++;
            } else {
                myLogger.error("readFasta: file: " + filename + " invalid format.");
                throw new IllegalArgumentException("Import: readFasta: invalid format.");
            }

        }

        String[] taxaNames = new String[taxa.size()];
        taxa.toArray(taxaNames);
        IdGroup idGroup = new SimpleIdGroup(taxaNames);

        String[] sequenceArray = new String[sequences.size()];
        sequences.toArray(sequenceArray);

        Locus unknown = new Locus("Unknown", "0", 0, sequenceArray[0].length(), null, null);
        return BitAlignment.getNucleotideInstance(idGroup, sequenceArray, null, null, null, TasselPrefs.getAlignmentMaxAllelesToRetain(), new Locus[]{unknown}, new int[]{0}, null, TasselPrefs.getAlignmentRetainRareAlleles(), isSBit);

    }

    public static Alignment readAlignmentFromSerialGZ(String inFile) {

        Alignment alignment = null;
        long time = System.currentTimeMillis();
        FileInputStream fis = null;
        GZIPInputStream gs = null;
        ObjectInputStream ois = null;
        try {
            File theFile = new File(Utils.addSuffixIfNeeded(inFile, ".serial.gz"));
            myLogger.info("readAlignmentFromSerialGZ: Reading:" + theFile);
            fis = new FileInputStream(theFile);
            gs = new GZIPInputStream(fis);
            ois = new ObjectInputStream(gs);
            alignment = (Alignment) ois.readObject();

        } catch (Exception ee) {
            ee.printStackTrace();
        } finally {
            try {
                ois.close();
                gs.close();
                fis.close();
            } catch (Exception e) {
                // do nothing
            }
        }
        myLogger.info("readAlignmentFromSerialGZ: Time: " + (System.currentTimeMillis() - time) + "  Sites: " + alignment.getSiteCount() + "  Taxa: " + alignment.getSequenceCount());
        return alignment;
    }

    public static Alignment readFromPLink(final String pedFilename, final String mapFilename, ProgressListener listener) {
        return readFromPLink(pedFilename, mapFilename, true, listener);
    }

    public static Alignment readFromPLink(final String pedFilename, final String mapFilename, boolean isSBit, ProgressListener listener) {

        int minPosition = Integer.MAX_VALUE;
        String currLocus = null;
        List<Locus> loci = new ArrayList<Locus>();
        List<Integer> lociOffsets = new ArrayList<Integer>();

        long currentTime = System.currentTimeMillis();
        int numSites = Utils.getNumberLines(pedFilename) - 1;
        myLogger.info("readFromHapmap: Number of Sites: " + numSites);

        long prevTime = currentTime;
        currentTime = System.currentTimeMillis();
        myLogger.info("readFromHapmap: Time to count lines: " + ((currentTime - prevTime) / 1000));


        BufferedReader fileIn = null;
        try {
            int numThreads = Runtime.getRuntime().availableProcessors();
            ExecutorService pool = Executors.newFixedThreadPool(numThreads);

            fileIn = Utils.getBufferedReader(pedFilename, 1000000);
            String[] header = WHITESPACE_PATTERN.split(fileIn.readLine());
            int lineInFile = 1;
            int numTaxa = header.length - NUM_HAPMAP_NON_TAXA_HEADERS;
            String[] snpIDs = new String[numSites];
            int prevPosition = -1;

            OpenBitSet[][] theData;
            byte[][] alleles = new byte[numSites][TasselPrefs.getAlignmentMaxAllelesToRetain()];
            int numDataRows = TasselPrefs.getAlignmentMaxAllelesToRetain();
            if (TasselPrefs.getAlignmentRetainRareAlleles()) {
                numDataRows++;
            }
            int numSitesToProcess = 1;
            if (isSBit) {
                theData = new OpenBitSet[numDataRows][numSites];
                numSitesToProcess = 1;
            } else {
                theData = new OpenBitSet[numDataRows][numTaxa];
                for (int al = 0; al < numDataRows; al++) {
                    for (int t = 0; t < numTaxa; t++) {
                        theData[al][t] = new OpenBitSet(numSites);
                    }
                }
                numSitesToProcess = 64;
            }

            int[] physicalPositions = new int[numSites];
            int count = 0;
            String[][] tokens = new String[numSitesToProcess][];
            int currentSite = 0;
            for (int site = 0; site < numSites; site++) {

                lineInFile++;

                String input = fileIn.readLine();
                tokens[count] = WHITESPACE_PATTERN.split(input);

                snpIDs[site] = new String(tokens[count][HAPMAP_SNPID_COLUMN_INDEX]);
                int position = Integer.parseInt(tokens[count][HAPMAP_POSITION_COLUMN_INDEX]);
                String temp = new String(tokens[count][HAPMAP_CHROMOSOME_COLUMN_INDEX]);
                if (currLocus == null) {
                    lociOffsets.add(site);
                    currLocus = temp;
                    minPosition = position;
                    prevPosition = -1;
                } else if (!temp.equals(currLocus)) {
                    loci.add(new Locus(currLocus, currLocus, minPosition, prevPosition, null, null));
                    lociOffsets.add(site);
                    currLocus = temp;
                    minPosition = position;
                    prevPosition = -1;
                }

                if (position < prevPosition) {
                    throw new IllegalStateException("ImportUtils: readFromHapmap: Sites are not properly sorted for chromosome: " + currLocus + " at " + position + " and " + prevPosition);
                }

                count++;

                if (count == numSitesToProcess) {
                    pool.execute(ProcessLineOfHapmap.getInstance(alleles, theData, TasselPrefs.getAlignmentRetainRareAlleles(), tokens, count, currentSite, numTaxa, lineInFile, isSBit));
                    count = 0;
                    currentSite += numSitesToProcess;
                    tokens = new String[numSitesToProcess][];
                }

                physicalPositions[site] = position;
                prevPosition = position;

                if (listener != null) {
                    listener.progress((int) (((double) (site + 1) / (double) numSites) * 100.0), null);
                }
            }

            if (count != 0) {
                pool.execute(ProcessLineOfHapmap.getInstance(alleles, theData, TasselPrefs.getAlignmentRetainRareAlleles(), tokens, count, currentSite, numTaxa, lineInFile, isSBit));
            }


            pool.shutdown();
            if (!pool.awaitTermination(6000, TimeUnit.SECONDS)) {
                throw new IllegalStateException("ImportUtils: readFromHapmap: processing threads timed out.");
            }

            if (currLocus != null) {
                loci.add(new Locus(currLocus, currLocus, minPosition, prevPosition, null, null));
            }

            prevTime = currentTime;
            currentTime = System.currentTimeMillis();
            myLogger.info("readFromHapmap: Time to read file: " + ((currentTime - prevTime) / 1000));

            String[] taxaNames = new String[numTaxa];
            System.arraycopy(header, NUM_HAPMAP_NON_TAXA_HEADERS, taxaNames, 0, numTaxa);
            IdGroup idGroup = new SimpleIdGroup(taxaNames);

            Locus[] lociFinal = new Locus[loci.size()];
            loci.toArray(lociFinal);
            int[] offsetsFinal = new int[lociOffsets.size()];
            for (int i = 0; i < lociOffsets.size(); i++) {
                offsetsFinal[i] = ((Integer) lociOffsets.get(i)).intValue();
            }

            Alignment result = BitAlignment.getNucleotideInstance(idGroup, alleles, theData, null, null, physicalPositions, TasselPrefs.getAlignmentMaxAllelesToRetain(), lociFinal, offsetsFinal, snpIDs, TasselPrefs.getAlignmentRetainRareAlleles(), isSBit);

            prevTime = currentTime;
            currentTime = System.currentTimeMillis();
            myLogger.info("readFromHapmap: Time to create Alignment: " + ((currentTime - prevTime) / 1000));

            return result;
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("ImportUtils: readFromHapmap: Problem creating Alignment: " + pedFilename + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                fileIn.close();
            } catch (Exception ex) {
                // do nothing
            }
        }

    }
}
