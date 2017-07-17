package net.maizegenetics.gbs.tagdist;

import java.util.TreeSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;

import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.OpenBitSet;

/**
 * This class contains methods which process TagsByTaxa files one line at a time, as opposed
 * to the methods in {@link AbstractTagsByTaxa}, which hold an entire file in memory.
 * @author jvh39 */
public class TagsByTaxaUtils {

    public static void main(String[] args) {
        printCoverage("/media/jvh053111/all_sorghum_20120326/tbt/70MU0AAXX_s_1.tbt.byte", TagsByTaxa.FilePacking.Byte, false);
//       printSumCounts("/media/SSD/all_sorghum_merged_032812.tbt.byte", TagsByTaxa.FilePacking.Byte, true);
//       printSumCountsOfAll("/media/jvh053111/all_sorghum_20120326/tbt/", TagsByTaxa.FilePacking.Byte);
//       printTotalTagsAndTaxa("/media/jvh053111/all_sorghum_20120326/tbt/", TagsByTaxa.FilePacking.Byte);
//       replaceInNames("\\s", "_", "/media/Data/Sorghum/merged_all_sorghum_081611.tbt.bin", TagsByTaxa.FilePacking.Bit, "/media/Cache/all_sorghum_spacesremoved_110819.tbt.bin");
    }

    private static class TBTRecord {

        long[] sequence;
        byte tagLength;
        long[] bitDistribution;
        byte[] byteDistribution;

        /** Returns a byte[] of taxon distribution for the current tag, regardless of file format.*/
        public byte[] taxonDist() {
            byte[] result = new byte[taxaNames.length];
            if (inputFormat == TagsByTaxa.FilePacking.Bit) {
                OpenBitSet currBitSet = new OpenBitSet(currRecord.bitDistribution, longsInBitset);
                for (int i = 0; i < currBitSet.capacity(); i++) {
                    if (currBitSet.fastGet(i)) {
                        result[i]++;
                    }
                }
            } else if (inputFormat == TagsByTaxa.FilePacking.Byte) {
                return byteDistribution;
            } else {
                System.out.println("Taxon distribution of the current tag is blank.");
                return null;
            }
            return result;
        }

        public float taxonCoverage() {
            long zeroValues = 0, nonZeroValues = 0;
            for (byte value : taxonDist()) {
                if (value > 0) {
                    zeroValues++;
                } else {
                    nonZeroValues++;
                }
            }
            return ((float) nonZeroValues) / ((float) zeroValues + nonZeroValues);
        }

        public TBTRecord() {
        }
    ;
    }

    /**These fields are declared static here so they aren't re-instantiated for every new record that's read.*/
    private static TBTRecord currRecord = new TBTRecord();
    private static int numTags;
    private static int tagLengthInLong;
    private static int numTaxa;
    private static int longsInBitset;
    private static String[] taxaNames;
    private static DataOutputStream outputStream;
    private static DataInputStream inputStream;
    private static BufferedReader reader;
    private static BufferedWriter writer;
    private static String inputFileName;
    private static String outputFileName;
    private static TagsByTaxa.FilePacking inputFormat, outputFormat;

    /**Prints bitsets in human-readable format for debugging purposes.
     *@param bitset An OpenBitSet object.
     *@return result A string of 1s and 0s representing the contents of the bitset.     */
    public static String bitsetToString(OpenBitSet bitset) {
        String result = "";
        for (int i = 0; i < bitset.size(); i++) {
            if (bitset.fastGet(i)) {
                result = result + "1";
            } else {
                result = result + "0";
            }
        }
        return result;
    }

    /**Stores bitsets in human-readable format.
     *@param bitset An OpenBitSet object.
     *@return result A string array of 1s and 0s representing the contents of the bitset.     */
    public static String[] bitsetToStringArray(OpenBitSet bitset) {
        String[] result = new String[(int) bitset.size()];
        for (int i = 0; i < bitset.size(); i++) {
            if (bitset.fastGet(i)) {
                result[i] = "1";
            } else {
                result[i] = "0";
            }
        }
        return result;
    }

    /**Converts a text .tbt file to binary, using only the memory required for the file buffer.*/
    public static void streamTextToBinary(String inputFileName, String outputFileName) {
        int tagNum, tagLengthInLong, taxaNum;
        try {
            System.out.println("Reading " + inputFileName + ".");
            BufferedReader br = new BufferedReader(new FileReader(inputFileName), 65536);
            int hapsOutput = 0; //Number of tags written to output file

            //Read header
            String[] inputLine = br.readLine().split("\t");
            tagNum = Integer.parseInt(inputLine[0]);
            tagLengthInLong = Integer.parseInt(inputLine[1]);
            taxaNum = Integer.parseInt(inputLine[2]);
            String[] taxaNames = new String[taxaNum];   //Taxa name array
            OpenBitSet obs = new OpenBitSet(taxaNum);   //Distribution bitset

            //Read taxon names
            inputLine = br.readLine().trim().split("\t");
            for (int t = 0; t < taxaNum; t++) {
                taxaNames[t] = inputLine[t];
            }

            //Write header
            System.out.println("Writing " + outputFileName + ".");
            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFileName), 65536));
            fw.writeInt(tagNum);                                                //Number of tags
            fw.writeInt(tagLengthInLong);                                       //Tag length
            fw.writeInt(taxaNum);                                               //Number of taxa
            for (int t = 0; t < taxaNum; t++) {                                 //Taxon names
                fw.writeUTF(taxaNames[t]);
            }

            //Loop through tags
            for (int i = 0; i < tagNum; i++) {
                inputLine = br.readLine().split("\t");                              //Read tag...
                long[] tagSequence = BaseEncoder.getLongArrayFromSeq(inputLine[0]); //Sequence
                byte tagLength = Byte.parseByte(inputLine[1]);                      //Length
                for (int t = 0; t < taxaNum; t++) {
                    obs.fastClear(t);                                               //Clear bit, in case it is still set
                    if (inputLine[t + 2].matches("1")) {
                        obs.set(t);
                    }                    //Distribution
                }

                for (int j = 0; j < tagLengthInLong; j++) {
                    fw.writeLong(tagSequence[j]); //Write sequence
                }
                fw.writeByte(tagLength);                                                //Write length
                long[] obsInLong = obs.getBits();                                         //Put bitset in array
                for (int t = 0; t < obsInLong.length; t++) {
                    fw.writeLong(obsInLong[t]);                                         //Write bitDistribution
                }
                hapsOutput++;
                if (hapsOutput % 10000 == 0) {
                    System.out.println("Wrote " + hapsOutput + " tags.");
                }
            }
            fw.close();
            br.close();

            System.out.println("Number of Taxa in file:" + taxaNum);
            System.out.println("Number of Haplotypes in file:" + hapsOutput);
        } catch (Exception e) {
            System.out.println("Catch in writing output file e=" + e);
            e.printStackTrace();
        }
    }

    /**Writes a single TBT record to disk.  It  uses the value of <b>format</b>
     * to determine whether to write text, binary, or another format.
     * @param currRecord  A TBTRecord object.*/
    private static void writeRecord(TBTRecord currRecord) {
        switch (outputFormat) {
            case Bit:
                try {
                    for (Long i : currRecord.sequence) {
                        outputStream.writeLong(i);
                    }
                    outputStream.write(currRecord.tagLength);
                    for (Long i : currRecord.bitDistribution) {
                        outputStream.writeLong(i);
                    }
                } catch (Exception e) {
                    System.out.println("Caught exception while writing binary TBT record: " + e);
                }
                break;

            case Text:
                //Records are output to text two different ways, depending on which format they are stored in.
                String outputLine = null;
                if (inputFormat.equals(TagsByTaxa.FilePacking.Bit)) {
                    outputLine =
                            BaseEncoder.getSequenceFromLong(currRecord.sequence) + "\t"
                            + Byte.toString(currRecord.tagLength) + "\t"
                            + bitsetToString(new OpenBitSet(currRecord.bitDistribution, currRecord.bitDistribution.length)) + "\n";
                } else if (inputFormat.equals(TagsByTaxa.FilePacking.Byte)) {
                    outputLine =
                            BaseEncoder.getSequenceFromLong(currRecord.sequence) + "\t"
                            + Byte.toString(currRecord.tagLength) + "\t";
                    for (int i = 0; i < currRecord.byteDistribution.length; i++) {
                        outputLine = outputLine + currRecord.byteDistribution[i] + "\t";
                    }
                    outputLine = outputLine + "\n";
                }
                try {
                    writer.write(outputLine);
                } catch (Exception e) {
                    System.out.println("Caught Exception while writing TBT text record: " + e);
                }
                break;

            case Byte:
                try {
                    for (Long i : currRecord.sequence) {
                        outputStream.writeLong(i);
                    }
                    outputStream.write(currRecord.tagLength);
                    outputStream.write(currRecord.byteDistribution);
                } catch (Exception e) {
                    System.out.println("Caught exception while writing binary TBT record: " + e);
                }
                break;
        }
    }

    /**Reads a single TBT record to disk.  It  uses the value of <b>format</b>
     * to determine whether to read text, binary, or another format.
     * @param currRecord  A TBTRecord object.*/
    private static void readRecord() {
        switch (inputFormat) {
            case Bit:
                try {
                    for (int i = 0; i < tagLengthInLong; i++) {
                        currRecord.sequence[i] = inputStream.readLong();
                    }
                    currRecord.tagLength = inputStream.readByte();
                    for (int i = 0; i < longsInBitset; i++) {
                        currRecord.bitDistribution[i] = inputStream.readLong();
                    }
                } catch (Exception e) {
                    System.out.println("Caught exception while reading binary TBT record: " + e);
                }
                break;

            case Byte:
                currRecord.byteDistribution = new byte[numTaxa];
                try {
                    for (int i = 0; i < tagLengthInLong; i++) {
                        currRecord.sequence[i] = inputStream.readLong();
                    }
                    currRecord.tagLength = inputStream.readByte();
                    for (int i = 0; i < numTaxa; i++) {
                        currRecord.byteDistribution[i] = inputStream.readByte();
                    }
                } catch (Exception e) {
                    System.out.println("Caught exception while reading byte-encoded TBT record: " + e);
                }
                break;

            case Text:
                OpenBitSet bitset = new OpenBitSet(numTaxa);
                try {
                    String[] currLine = reader.readLine().split("\\s");
                    currRecord.sequence = BaseEncoder.getLongArrayFromSeq(currLine[0]);
                    currRecord.tagLength = Byte.parseByte(currLine[1]);

                    for (int i = 0; i < numTaxa; i++) {
                        if (currLine[i + 2].equals("1")) {
                            bitset.fastSet(i);
                        }
                    }
                    currRecord.bitDistribution = bitset.getBits();
                } catch (Exception e) {
                    System.out.println("Caught exception while reading text TBT record: " + e);
                }
                break;

            default:
                System.out.println("Couldn't read TBT record: file format is unknown.");
        }
    }

    /**Reads a TBT file header and initializes class variables that depend on
     * the information in the header.  It  uses the value of <b>format</b>
     * to determine whether to read text, binary, or another format.*/
    private static void readHeader() {
        switch (inputFormat) {
            case Bit:   //Fall through
            case Byte:
                try {
                    numTags = inputStream.readInt();
                    tagLengthInLong = inputStream.readInt();
                    numTaxa = inputStream.readInt();

                    taxaNames = new String[numTaxa];
                    for (int i = 0; i < numTaxa; i++) {
                        taxaNames[i] = inputStream.readUTF();
                    }
                    break;
                } catch (Exception e) {
                    System.out.println("Caught exception while reading binary TBT file header: " + e);
                    System.exit(0);
                }

            case Text:
                try {
                    String[] headerLine = reader.readLine().split("\\s");
                    numTags = Integer.parseInt(headerLine[0]);
                    tagLengthInLong = Integer.parseInt(headerLine[2]);
                    numTaxa = Integer.parseInt(headerLine[2]);

                    taxaNames = reader.readLine().split("\\s");
                    break;
                } catch (Exception e) {
                    System.out.println("Caught exception while reading TBT text file header: " + e);
                }
            default:
                System.out.println("Couldn't read header: the file format isn't recognized.");
                System.exit(0);
                break;
        }
        longsInBitset = BitUtil.bits2words(numTaxa);
        currRecord.sequence = new long[tagLengthInLong];
        currRecord.bitDistribution = new long[longsInBitset];
    }

    /**Writes a TBT file header, using the values of class variables that
     * pertain to the header.  It  uses the value of <b>format</b>
     * to determine whether to write text, binary, or another format.*/
    private static void writeHeader() {
        switch (outputFormat) {
            case Bit:   //Fall through
            case Byte:
                try {
                    outputStream.writeInt(numTags);
                    outputStream.writeInt(tagLengthInLong);
                    outputStream.writeInt(numTaxa);
                    for (String name : taxaNames) {
                        outputStream.writeUTF(name);
                    }
                } catch (Exception e) {
                    System.out.println("Caught exception while writing binary TBT file header: " + e);
                }
                break;

            case Text:
                try {
                    writer.write(numTags + "\t" + tagLengthInLong + "\t" + numTaxa + "\n");
                    for (String name : taxaNames) {
                        writer.write(name + "\t");
                    }
                    writer.write("\n");
                } catch (Exception e) {
                    System.out.println("Caught exception while writing text TBT file header: " + e);
                }
                break;

            default:
                System.out.println("Couldn't print header: the file format isn't recognized.");
                System.exit(0);
        }
    }

    public static void slice(String inputFileName, String[] sequences) {
        throw new UnsupportedOperationException("This function isn't finished yet.\n");
//        long[][] binarySequences= new long[tagLengthInLong][sequences.length];
//        for (int i= 0; i < sequences.length; i++) {
//            long[] currSeq=BaseEncoder.getLongArrayFromSeq(sequences[i]);
//            for (int j= 0; j < tagLengthInLong; j++) {
//                binarySequences[j][i]=currSeq[j];
//            }
//        }
//
//        for (int i= 0; i < numTags; i++) {
//            readRecord();
//            for (int j= 0; j < binarySequences[0].length; j++) {
//
//            }
//        }
    }

    /**Outputs taxon bitDistribution of all tags in the supplied TOPM file.
     *
     * @param inputFileName
     * @param topm
     */
    public static void filterByTOPM(String inputFileName, TagsOnPhysicalMap topm) {
        outputFileName = inputFileName.replaceFirst("tbt", "tbt.filteredbytopm");

        openIOStreams(inputFileName, TagsByTaxa.FilePacking.Byte, outputFileName, TagsByTaxa.FilePacking.Byte);
        readHeader();

        //Count tags present in topm
        int filteredLines = 0, unfilteredLines = numTags;
        for (int i = 0; i < unfilteredLines; i++) {
            readRecord();
            if (topm.getTagIndex(currRecord.sequence) > 0) {
                filteredLines++;
            }
            if (unfilteredLines % 1000000 == 0) {
                System.out.println("Counted " + (float) filteredLines / 1000000 + " million records.");
            }
        }
        closeIOStreams();

        openIOStreams(inputFileName, TagsByTaxa.FilePacking.Byte, outputFileName, TagsByTaxa.FilePacking.Byte);
        readHeader();
        //Update header with new (reduced) number of tags, then write file
        numTags = filteredLines;
        int writtenLines = 0;
        writeHeader();
        for (int i = 0; i < unfilteredLines; i++) {
            readRecord();
            if (topm.getTagIndex(currRecord.sequence) > 0) {
                writeRecord(currRecord);
                writtenLines++;
                if (writtenLines % 1000000 == 0) {
                    System.out.println("Wrote " + (float) writtenLines / 1000000 + " million records.");
                }
            }
        }
        System.out.println("Wrote " + writtenLines + " records.");
        closeIOStreams();
    }

    /** Replaces the supplied regex with the supplied string in taxon names.
     *
     * @param regex
     * @param replacement
     * @param inputFileName
     * @param format    Format of the input .tbt file (a TagsByTaxa.FilePacking enumerated value)
     * @param outputFileName
     */
    public static void replaceInNames(String regex, String replacement, String inputFileName, TagsByTaxa.FilePacking format, String outputFileName) {
        openIOStreams(inputFileName, format, outputFileName, format);
        readHeader();
        for (int i = 0; i < taxaNames.length; i++) {
            System.out.println("Replacing " + taxaNames[i] + " with ");
            taxaNames[i] = taxaNames[i].replaceAll(regex, replacement);
            System.out.println(taxaNames[i] + ".");
        }
        writeHeader();
        for (int i = 0; i < numTags; i++) {
            readRecord();
            writeRecord(currRecord);
        }
    }

    /**Prints the physical coordinates of the supplied taxa, as stored in the supplied TOPM file.
     *
     * @param inputFileName A reference to a TagsByTaxa file.
     * @param format    A TagsByTaxa.FilePacking enumerated value.
     * @param topm  A TagsOnPhysicalMap object containing coordinate info for the supplied taxa.
     * @param taxa  A list of taxon names.
     */
    public static void positions(String inputFileName, TagsByTaxa.FilePacking format, TagsOnPhysicalMap topm, String[] taxa) {
        outputFileName = inputFileName.replaceFirst("tbt.bin", "positions.tbt.bin");
        OpenBitSet currBitset;
        int linesPrinted = 0;
        ArrayList<Integer> columnsTracked = new ArrayList<Integer>();
        openIOStreams(inputFileName, format, outputFileName, format);
        readHeader();

        //Compare provided names to names in TBT.  Track the columns of ones that match.
        for (int i = 0; i < taxaNames.length; i++) {
            for (String taxon : taxa) {
                if (taxaNames[i].equals(taxon)) {
                    columnsTracked.add(i);
                }
            }
        }

        //Print header of output file
        try {
            outputStream.writeBytes("Printing positions of following taxa:\n");
            for (Integer integer : columnsTracked) {
                String taxon = taxaNames[integer];
                outputStream.writeBytes(taxon + "\n");
            }
        } catch (Exception e) {
            System.out.println("Caught exception while printing positions of taxa: " + e);
        }

        //Loop over tags and put bitDistribution into bitset
        for (int i = 0; i < numTags; i++) {
            readRecord();
            currBitset = new OpenBitSet(currRecord.bitDistribution, longsInBitset);

            //Loop over tracked columns.  If the current sequence is present in the column,
            //output its position as found in the TOPM file.
            for (Integer integer : columnsTracked) {
                if (currBitset.fastGet(integer)) {

                    //Write coordinates from TOPM record
                    int topmIndex = topm.getTagIndex(currRecord.sequence);
                    try {
                        outputStream.writeBytes(topm.printRow(topmIndex) + "\n");
                    } catch (Exception e) {
                        System.out.println("Caught exception while printing positions of taxa: " + e);
                    }
                }
            }
            linesPrinted++;
        }
        closeIOStreams();
    }

    /**Creates new reader & writer class variables.  This function is just here to abstract away a repetitive task.*/
    private static void openIOStreams(String inputFileName, TagsByTaxa.FilePacking inputFormat, String outputFileName, TagsByTaxa.FilePacking outputFormat) {
        TagsByTaxaUtils.inputFormat = inputFormat;
        TagsByTaxaUtils.outputFormat = outputFormat;
        TagsByTaxaUtils.inputFileName = inputFileName;
        if (outputFileName == null) {
            outputFileName = inputFileName + ".output";
        } else {
            TagsByTaxaUtils.outputFileName = outputFileName;
        }

        try {
            switch (inputFormat) {
                case Bit:   //Fall through
                case Byte:
                    inputStream = new DataInputStream(new BufferedInputStream(new FileInputStream(inputFileName), 65536));
                    break;
                case Text:
                    reader = new BufferedReader(new FileReader(inputFileName), 65536);
                    break;
                default:
                    System.out.println("Couldn't determine input file format.");
                    System.exit(0);
                    break;
            }

            switch (outputFormat) {
                case Bit:   //Fall through
                case Byte:
                    outputStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFileName), 65536));
                    break;
                case Text:
                    writer = new BufferedWriter(new FileWriter(outputFileName), 65536);
                    break;
                default:
                    System.out.println("Couldn't determine input file format.");
                    System.exit(0);
                    break;
            }
        } catch (Exception e) {
            System.out.println("Caught exception while opening input & output files in TagsByTaxaUtils: " + e);
        }
    }

    /**Shuts down open reader & writer class variables.  This function is just here to abstract away a repetitive task.*/
    private static void closeIOStreams() {
        try {
            switch (inputFormat) {
                case Byte: //Fall through
                case Bit:
                    inputStream.close();
                    break;
                case Text:
                    reader.close();
                    break;
            }

            switch (outputFormat) {
                case Byte: //Fall through
                case Bit:
                    outputStream.close();
                    break;
                case Text:
                    writer.close();
                    break;
            }

        } catch (Exception e) {
            System.out.println("Caught exception while opening input & output files in TagsByTaxaUtils: " + e);
        }
    }

    /**Returns a list of taxa in the specified TagsByTaxa file, along with the total count of all tags found in that taxon.
     *
     * @param inputFileName
     * @param format  A TagsByTaxa.FilePacking enumerated value.
     * @param progressIndication    Whether or not to provide feeback on number of tags read.
     * @return  A hash map of tag counts indexed by taxon names.
     */
    public static HashMap<String, Integer> sumCounts(String inputFileName, TagsByTaxa.FilePacking format, boolean progressIndication) {
        HashMap<String, Integer> results = new HashMap<String, Integer>();
        outputFileName = inputFileName + ".tmp";
        openIOStreams(inputFileName, format, outputFileName, format);
        readHeader();
        int[] tagCount = new int[numTaxa];

        //Loop over tags and sum counts
        for (int i = 0; i < numTags; i++) {
            if (progressIndication && (i % 100000 == 1)) {
                System.out.println("Finished reading " + i + "tags...");
            }
            readRecord();

            if (format == TagsByTaxa.FilePacking.Bit) {
                OpenBitSet bitset = new OpenBitSet(currRecord.bitDistribution, numTaxa); //(currRecord.bitDistribution, numTaxa);
                for (int j = 0; j < numTaxa; j++) {
                    if (bitset.fastGet(j)) {
                        tagCount[j]++;
                    }
                }

            } else if (format == TagsByTaxa.FilePacking.Byte) {
                for (int j = 0; j < numTaxa; j++) {
                    tagCount[j] += currRecord.byteDistribution[j];
                }

            } else {
                System.out.println("TBT file format not recognized.");
            }
        }

        //Setup hashmap
        for (int i = 0; i < numTaxa; i++) {
            results.put(taxaNames[i], tagCount[i]);
        }
        return results;
    }

    /**Prints a count of the total number of tags and taxa in the files contained in the specified directory.*/
    public static void printTotalTagsAndTaxa(String directoryName, TagsByTaxa.FilePacking format) {
        int totalTags = 0, totalTaxa = 0;
        String[] filenames = DirectoryCrawler.listFileNames(".*.tbt.bin|.*.tbt.txt|.*.tbt.byte", directoryName);
        for (String filename : filenames) {
            openIOStreams(filename, format, null, format);
            readHeader();
            totalTags += numTags;
            totalTaxa += numTaxa;
            System.out.println("Reading " + inputFileName + ".  Counted " + numTags + " tags over " + numTaxa + " taxa.");
        }
        System.out.println("Total: " + totalTags + " non-unique tags over " + totalTaxa + " non-unique taxa.");
    }

    /**Prints a list of taxa in the specified TagsByTaxa file, along with the total count of all tags found in that taxon.
     *
     * @param inputFileName
     * @param format  A TagsByTaxa.FilePacking enumerated value.
     * @param progressIndication    Whether or not to provide feedback on number of tags read.
     */
    public static void printSumCounts(String inputFileName, TagsByTaxa.FilePacking format, boolean progressIndication) {
        HashMap<String, Integer> results = sumCounts(inputFileName, format, progressIndication);
        for (String index : results.keySet()) {
            System.out.println(index + "\t" + results.get(index));
        }
    }

    /**Calls <b>printSumCounts</b> once for every file in the specified directory.
     *
     * @param inputFileName
     * @param format  A TagsByTaxa.FilePacking enumerated value.
     * @param progressIndication    Whether or not to provide feeback on number of tags read.   */
    public static void printSumCountsOfAll(String directoryName, TagsByTaxa.FilePacking format) {
        for (String filename : DirectoryCrawler.listFileNames(".*.tbt.bin|.*.tbt.txt|.*.tbt.byte", directoryName)) {
            System.out.println(filename + ":");
            printSumCounts(filename, format, false);
        }
    }

    /**Write out only the taxa that match a list of names*/
//           public static void selectByTaxa(String[] regexes, String inputFileName, TagsByTaxa.FilePacking format){
//               OpenBitSet outputThisTaxon = new OpenBitSet(numTaxa);
//               openIOStreams(inputFileName, format);
//               readHeader();
//
//               //Loop over taxon names and regexes.  If a given taxon name matches any of the regexes,
//               //it will be marked for output later.
//               for (int i= 0; i < numTaxa; i++) {
//                    for (String regex: regexes) {
//                       if(taxaNames[i].matches(regex)){
//                           outputThisTaxon.fastSet(i);
//                           break;
//                       }
//                   }
//               }
//           }

    /*Prints out statistics on the sparsity of a TBT matrix (i.e., how many tag counts are zero).*/
    public static void sparsity(String inputFileName, TagsByTaxa.FilePacking format) {
        long zeroValues = 0;
        long nonZeroValues = 0;
        OpenBitSet currBitSet;

        openIOStreams(inputFileName, format, inputFileName + ".tmp", format);
        readHeader();
        System.out.println("Measuring TBT matrix sparsity.");
        for (int i = 0; i < numTags; i++) {
            readRecord();
            if (i % 1000000 == 0) {
                System.out.println("Read " + i / 1000000 + " million tags.");
            }
            if (format.equals(TagsByTaxa.FilePacking.Bit)) {
                currBitSet = new OpenBitSet(currRecord.bitDistribution, longsInBitset);
                nonZeroValues += currBitSet.cardinality();
                zeroValues += currBitSet.size() - currBitSet.cardinality();
            } else if (format.equals(TagsByTaxa.FilePacking.Byte)) {
                for (byte value : currRecord.byteDistribution) {
                    if (value == 0) {
                        zeroValues++;
                    } else {
                        nonZeroValues++;
                    }
                }
            } else {
                System.out.println("This format isn't recognized by the sparsity function.");
                System.exit(0);
            }
        }
        long totalValues = zeroValues + nonZeroValues;
        System.out.println(
                "\nTags: " + numTags + "\n"
                + "Taxa: " + numTaxa + "\n"
                + "Size: " + (numTags * numTaxa) + " values.\n"
                + zeroValues + " zero values.\n"
                + nonZeroValues + " non-zero values.\n"
                + "Sparsity is " + (float) zeroValues / ((float) totalValues) + ".");
    }

    /** Prints tag coverage of each taxon in the file, and taxon coverage of each tag in the file.
     * If itemized=true, prints one record for each tag and for each taxon.
     *  Otherwise prints a summary of how many tags/taxa are in each coverage "bin" (0-10% coverage, 10-20% coverage, etc.)*/
    public static void printCoverage(String inputFileName, TagsByTaxa.FilePacking format, boolean itemized) {
        openIOStreams(inputFileName, format, inputFileName + ".tmp", format);
        readHeader();
        int[] taxonCoverage = new int[10], tagCoverage = new int[10],
                tagsCovered = new int[taxaNames.length], tagsUncovered = new int[taxaNames.length];

        //Print taxon coverage line by line
        if (itemized) {
            System.out.println("Total " + taxaNames.length + " taxa." + '\n' + "sequence" + '\t' + "taxon coverage");
        }
        for (int i = 0; i < numTags; i++) {
            readRecord();

            //Determine taxon coverage
            float coverage = currRecord.taxonCoverage();
            int bin = ((int) (coverage * 10));
            taxonCoverage[bin]++;
            if (itemized) {
                System.out.println(BaseEncoder.getSequenceFromLong(currRecord.sequence) + '\t' + coverage);
            }

            //Store counts/taxon now, determine tag coverage later
            byte[] currCounts = currRecord.taxonDist();
            for (int j = 0; j < tagsCovered.length; j++) {
                if (currCounts[j] > 0) {
                    tagsCovered[j]++;
                } else {
                    tagsUncovered[j]++;
                }
            }
        }

        //Summarize taxon coverage
        float total = 0;
        if (!itemized) {
            System.out.println("Taxon Coverage:");
        }
        for (int count : taxonCoverage) {
            total += ((float) count);
        }
        if (!itemized) {
            System.out.println("percent of taxa covered" + '\t' + "tags at coverage level" + '\t' + "fraction of total tags");
            for (int bin = 0; bin < 10; bin++) {
                int min = bin * 10, max = (bin + 1) * 10;
                float pct = ((float) taxonCoverage[bin]) / total;
                System.out.println(min + "-" + max + "%:" + '\t' + taxonCoverage[bin] + '\t' + pct);
            }
        }


        //Determine tag coverage
        if (itemized) {
            System.out.println("Total " + numTags + " tags." + '\n' + "taxon" + '\t' + "tag coverage" + '\t' + "fraction of total");
        }
        for (int i = 0; i < taxaNames.length; i++) {
            float coverage = ((float) tagsCovered[i]) / ((float) numTags);
            int bin = ((int) (coverage * 10));
            tagCoverage[bin]++;
            if (itemized) {
                System.out.println(taxaNames[i] + '\t' + tagsCovered[i] + '\t' + coverage);
            }
        }

        //Summarize tag coverage
        if (!itemized) {
            System.out.println("Tag Coverage:");
        }
        if (!itemized) {
            System.out.println("percent of tags covered" + '\t' + "taxa at coverage level" + '\t' + "fraction of total taxa");
            for (int bin = 0; bin < 10; bin++) {
                int min = bin * 10, max = (bin + 1) * 10;
                float pct = ((float) tagCoverage[bin]) / ((float) taxaNames.length);
                System.out.println(min + "-" + max + "%:" + '\t' + tagCoverage[bin] + '\t' + pct);
            }
        }
    }

    public static TagsByTaxa.FilePacking format(String filename) {
        TagsByTaxa.FilePacking format = null;
        String fileExtension = filename.substring(filename.lastIndexOf("."));
        if (fileExtension.equals(".bin")) {
            format = TagsByTaxa.FilePacking.Bit;
        }
        if (fileExtension.equals(".byte")) {
            format = TagsByTaxa.FilePacking.Byte;
        }
        if (fileExtension.equals(".shrt")) {
            format = TagsByTaxa.FilePacking.Short;
        }
        if (format == null) {
            System.out.println("Couldn't identify format of input file: " + filename);
            System.exit(0);
        }
        return format;
    }

    /** Merge taxa with identical names (merge their (binary) tag counts into a single column).
     * 
     * The input must be a TagsByTaxaBit file in either binary or text format. The output is written in 
     * binary TagsByTaxaBit format. Other TagsByTaxa formats (Byte, Short, etc) are not currently supported. */
    public static void mergeTaxaByName(String inputFileName, String outputFileName, TagsByTaxa.FilePacking format, boolean caseSensitive) {
        openIOStreams(inputFileName, format, outputFileName, format);
        readHeader();

        //Get non-unique portion of name; merge into a hash table
        TreeSet<String> uniqueNames = new TreeSet<String>();
        OpenBitSet oldBitset, newBitset;
        for (int i = 0; i < taxaNames.length; i++) {
            if (caseSensitive) {
                taxaNames[i] = taxaNames[i].substring(0, taxaNames[i].indexOf(":"));
            } else {
                taxaNames[i] = taxaNames[i].substring(0, taxaNames[i].indexOf(":")).toUpperCase();
            }
            uniqueNames.add(taxaNames[i]);
        }

        //Create sorted array of new unique names
        String[] nameArray = uniqueNames.toArray(new String[uniqueNames.size()]);
        Arrays.sort(nameArray);

        // Loop over old columns. For each element in taxaNames, the corresponding 
        // element in newColumn indicates the new column where it belongs.
        int[] newColumn = new int[taxaNames.length];
        for (int i = 0; i < taxaNames.length; i++) {
            for (int j = 0; j < nameArray.length; j++) {
                if (taxaNames[i].equals(nameArray[j])) {
                    newColumn[i] = j;
                }
            }
        }

        //Update file properties to reflect new file size; write header.
        taxaNames = Arrays.copyOf(nameArray, nameArray.length);
        numTaxa = taxaNames.length;
        writeHeader();

        //Loop over tags; if tag is present in a non-unique taxon, set its presence 
        //in the corresponding unique taxon, based on the value in newColumn[].
        for (int i = 0; i < numTags; i++) {
            readRecord();
            oldBitset = new OpenBitSet(currRecord.bitDistribution, longsInBitset);
            newBitset = new OpenBitSet(uniqueNames.size());
            for (int j = 0; j < oldBitset.capacity(); j++) {
                if (oldBitset.fastGet(j)) {
                    newBitset.fastSet(newColumn[j]);
                }
            }

            //Update current record to reflect new tag bitDistribution, write record, and then
            //change bitDistribution back to the previous size to read another old record.
            currRecord.bitDistribution = newBitset.getBits();
            writeRecord(currRecord);
            currRecord.bitDistribution = new long[longsInBitset];
        }
        closeIOStreams();
    }

    /**Converts a binary .tbt file to text, using only the memory required for the file buffer.*/
    public static void streamBinaryToText(String inputFileName, int maxRecords) {
        openIOStreams(inputFileName, format(inputFileName), inputFileName.replaceFirst("\\.bin$|\\.byte$", ".txt"), TagsByTaxa.FilePacking.Text);
        readHeader();
        writeHeader();
        int writtenLines = 0;
        for (int i = 0; i < numTags && i < maxRecords; i++) {
            readRecord();
            writeRecord(currRecord);
            writtenLines++;
        }
        System.out.println("Wrote " + writtenLines + " records.");
        closeIOStreams();
    }
}
