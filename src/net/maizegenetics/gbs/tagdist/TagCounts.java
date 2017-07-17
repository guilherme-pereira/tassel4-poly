package net.maizegenetics.gbs.tagdist;

import cern.colt.GenericSorting;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;

import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.util.BaseEncoder;

/**
 * Holds tags counts.  Tags sequences are compressed in long, tags lengths are tracked,
 * and read counts are stored.  Has basic filtering methods.
 * Read counts can be modified.
 *
 * This class is not resizable in terms of new tags - for growing tag count lists use TagCountMutable.java
 *
 * User: ed
 */
public class TagCounts extends AbstractTags {

    int[] readCount;

    public TagCounts() {
        this.tagLengthInLong = 2;
        initMatrices(10);
    }

    public TagCounts(int tagLengthInLong, int maxSize) {
        this.tagLengthInLong = tagLengthInLong;
        initMatrices(maxSize);
    }

    public TagCounts(String inFile, FilePacking binary) {
        readDistFile(inFile, binary);
    }

    public void readDistFile(String infile, FilePacking binary) {
        File in = new File(infile);
        System.out.println("Reading Haplotypes distribution from:" + in.toString());
        switch (binary) {
            case Text:
                readTextTagCountFile(in);
                break;
            default:
                readBinaryTagCountFile(in);
                break;
        }
    }

    protected void initMatrices(int tagNum) {
        tags = new long[tagLengthInLong][tagNum];
        tagLength = new byte[tagNum];
        readCount = new int[tagNum];
    }
    
    /**
     * Convert TagCounts to FASTA file for alignment using Blast
     * @param outfile 
     */
    public void toFASTA (String outfile) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfile), 65536);
            for (int i = 0; i < this.getTagCount(); i++) {
                bw.write(">"+String.valueOf(i));
                bw.newLine();
                bw.write(BaseEncoder.getSequenceFromLong(this.getTag(i)).substring(0, this.getTagLength(i)));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /** Reads in a binary tag count file line-by-line and constructs a FASTQ file for alignment. */
    public void toFASTQ(String infile, String outfile) {
        int tagsRead = 0;
        long[] currSequence = new long[tagLengthInLong];
        String textSequence;
        byte currLength;
        int currCount;

        try {
            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfile), 65536));
            DataInputStream rw = new DataInputStream(
                    new BufferedInputStream(
                    new FileInputStream(infile),
                    4000000));

            //Read header
            int tagNum = rw.readInt();
            tagLengthInLong = rw.readInt();

            while (rw.available() != 0) {

                //Read single record
                for (int j = 0; j < tagLengthInLong; j++) {
                    currSequence[j] = rw.readLong();
                }

                currLength = rw.readByte();
                currCount = rw.readInt();
                tagsRead++;
                textSequence = BaseEncoder.getSequenceFromLong(currSequence);
                textSequence = textSequence.substring(0, currLength);  //Remove any poly-A padding
                //Write to FASTQ file with bogus quality score
                fw.writeBytes(
                        "@length=" + currLength + "count=" + currCount + "\n"
                        + textSequence + "\n"
                        + "+\n"
                        + "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff\n");
                
            }
            rw.close();
            fw.close();


        } catch (Exception e) {
            System.out.println("Catch in reading TagCount file e=" + e);


        }
        System.out.println("Number of Tags in file:" + tagsRead);

    }

    void readTextTagCountFile(File inFile) {
        int tagsRead = 0;
        String[] inputLine;


        try {
            BufferedReader br = new BufferedReader(new FileReader(inFile), 65536);
            inputLine = br.readLine().trim().split("\\s");


            int tagNum = Integer.parseInt(inputLine[0]);
            tagLengthInLong = Integer.parseInt(inputLine[1]);
            initMatrices(
                    tagNum);


            for (int i = 0; i < tagNum; i++) {
                inputLine = br.readLine().split("\\s");


                long[] tt = BaseEncoder.getLongArrayFromSeq(inputLine[0]);


                for (int j = 0; j
                        < tt.length; j++) {
                    tags[j][i] = tt[j];


                }
                tagLength[i] = Byte.valueOf(inputLine[1]);
                readCount[i] = Integer.valueOf(inputLine[2]);
                tagsRead++;

            }


            br.close();


        } catch (Exception e) {
            System.out.println("Catch in reading TagCount file e=" + e);
            e.printStackTrace();


        }
        System.out.println("Number of Tags in file:" + tagsRead);
        System.out.println("Number of Tags in memory:" + tags[0].length);


    }

    public void readBinaryTagCountFile(File inFile) {
        int tagsRead = 0;
        try {
            DataInputStream rw = new DataInputStream(new BufferedInputStream(new FileInputStream(inFile), 4000000));
            int tagNum = rw.readInt();
            tagLengthInLong = rw.readInt();
            initMatrices(tagNum);
            for (int i = 0; i < tagNum; i++) {
                for (int j = 0; j < tagLengthInLong; j++) {
                    tags[j][i] = rw.readLong();
                }
                tagLength[i] = rw.readByte();
                readCount[i] = rw.readInt();
                tagsRead++;
            }
            rw.close();
        } catch (Exception e) {
            System.out.println("Catch in reading TagCount file e=" + e);


        }
        System.out.println("Number of Tags in file:" + tagsRead);


    }

    public void writeTagCountFile(String outFile, FilePacking binary, int minCount) {
        int hapsOutput = 0;
        int[] outTagsAndReadsKept = tagsWCountsGreaterThanMin(minCount);
        System.out.println(outTagsAndReadsKept[0] + " tags will be output to " + outFile);
        System.out.println("These " + outTagsAndReadsKept[0] + " tags were covered by " + outTagsAndReadsKept[1] + " matching reads");
        switch (binary) {
            case Text:
                hapsOutput = writeTextTagCountFile(outFile, outTagsAndReadsKept[0], minCount);
                break;
            default:
                hapsOutput = writeByteTagCountFile(outFile, binary, outTagsAndReadsKept[0], minCount);
                break;
        }
        System.out.println("Tags written to:" + outFile.toString());
        System.out.println("Number of Tags in file:" + hapsOutput);
    }

    private int writeTextTagCountFile(String outFile, int outReads, int minCount) {
        int hapsOutput = 0;
        try {
            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile), 65536));
            fw.writeBytes(outReads + "\t" + tagLengthInLong + "\n");

            for (int i = 0; i < tags[0].length; i++) {
                if (getReadCount(i) < minCount) {
                    continue;
                }
                fw.writeBytes(BaseEncoder.getSequenceFromLong(getTag(i)) + "\t");
                fw.writeBytes(getTagLength(i) + "\t");
                fw.writeBytes(getReadCount(i) + "\t");
                fw.writeBytes("\n");
                hapsOutput++;
            }
            fw.close();
        } catch (Exception e) {
            System.out.println("Catch in writeTextDistFile writing output file e=" + e);
            e.printStackTrace();
        }
        return hapsOutput;
    }

    private int writeByteTagCountFile(String outFile, FilePacking binary, int outReads, int minCount) {
        int hapsOutput = 0;
        try {
            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile), 65536));
            fw.writeInt(outReads);
            fw.writeInt(tagLengthInLong);
            for (int i = 0; i < tags[0].length; i++) {
                if (getReadCount(i) < minCount) {
                    continue;
                }
                for (int j = 0; j < tagLengthInLong; j++) {
                    fw.writeLong(tags[j][i]);
                }
                fw.writeByte(tagLength[i]);
                fw.writeInt(getReadCount(i));
                hapsOutput++;
            }
            fw.close();
        } catch (Exception e) {
            System.out.println("Catch in writeTextDistFile writing output file e=" + e);
            e.printStackTrace();
        }
        return hapsOutput;
    }

    /** @param index Array index of the read whose count this function should return. */
    public int getReadCount(int index) {
        if (index >= readCount.length) {
            return -1;


        }
        return readCount[index];


    }

    protected int[] tagsWCountsGreaterThanMin(int minCount) {
        int nPassingTags = 0;
        int nKeptReads = 0;
        for (int i = 0; i < getTagCount(); i++) {
            int readCount = getReadCount(i);
            if (readCount >= minCount) {
                nPassingTags++;
                nKeptReads += readCount;
            }
        }
        int[] retArr = {nPassingTags, nKeptReads};
        return retArr;
    }

    public void setTag(long[] sequence, byte length, int count, int index) {
        for (int i = 0; i < tagLengthInLong; i++) {
            tags[i][index] = sequence[i];
            tagLength[index] = length;
            readCount[index] = count;
        }
    }

    public void sort() {
        System.out.println("Position index sort begin.");
        GenericSorting.quickSort(0, this.getSize(), this, this);
        System.out.println("Position index sort end.");


    }

    protected void printRows(int numRows) {
        for (int i = 0; i
                < numRows; i++) {
            System.out.println(BaseEncoder.getSequenceFromLong(tags[0][i])
                    + BaseEncoder.getSequenceFromLong(tags[1][i])
                    + " " + this.getTagLength(i) + " "
                    + getReadCount(i));


        }
    }

    protected void collapseCounts() {
        sort();//requires that the reads are sorted
        int collapsedRows = 0;


        for (int i = 1; i
                < this.getSize(); i++) {
            if ((tags[0][i - 1] == tags[0][i]) && (tags[1][i - 1] == tags[1][i])) {
                readCount[i] += readCount[i - 1];
                readCount[i - 1] = 0;
                collapsedRows++;

            }


        }
        System.out.println("Rows collapsed:" + collapsedRows);


    }

    public int getSize() {
        return tags[0].length;


    }

    public int getTotalCount() {
        int totalCount = 0;


        for (int i = 0; i
                < readCount.length; i++) {
            totalCount += readCount[i];


        }
        return totalCount;


    }

    @Override
    public void swap(int index1, int index2) {
        long temp;
        for (int i = 0; i < tagLengthInLong; i++) {
            temp = tags[i][index1];
            tags[i][index1] = tags[i][index2];
            tags[i][index2] = temp;
        }
        byte tl;
        tl = tagLength[index1];
        tagLength[index1] = tagLength[index2];
        tagLength[index2] = tl;
        int t3;
        t3 = readCount[index1];
        readCount[index1] = readCount[index2];
        readCount[index2] = t3;
    }

    @Override
    public int compare(int index1, int index2) {
        for (int i = 0; i
                < tagLengthInLong; i++) {
            if (tags[i][index1] < tags[i][index2]) {
                return -1;


            }
            if (tags[i][index1] > tags[i][index2]) {
                return 1;


            }
        }
        if (readCount[index1] < readCount[index2]) {
            return -1;


        }
        if (readCount[index1] > readCount[index2]) {
            return 1;


        }
        return 0;


    }

    public byte compare(long[] sequence1, long[] sequence2) {
        for (int i = 0; i < tagLengthInLong; i++) {
            if (sequence1[i] < sequence2[i]) {
                return -1;
            }
            if (sequence1[i] > sequence2[i]) {
                return 1;
            }
        }
        return 0;
    }

    @Override
    public int[] getTagIndexSet(long[] read) {
        throw new UnsupportedOperationException("Not supported yet.");


    }

    @Override
    public boolean areTagsUnique() {
        throw new UnsupportedOperationException("Not supported yet.");


    }

    @Override
    public int getTagSizeInLong() {
        return tagLengthInLong;


    }

    @Override
    public int getTagLength(int index) {
        return (int) tagLength[index];


    }

    /** I don't know why we have 4 or 5 methods with different
     * names that do exactly the same thing.  I'm just putting this here
     * because it's defined in the abstract class and its absence is causing
     * problems.
     * @return The number of tags in the current object's tags[][] matrix.
     */
    @Override
    public int getTagCount() {
        return tags[0].length;
    }

    public long[][] getTags() {
        return tags;

    }
}
