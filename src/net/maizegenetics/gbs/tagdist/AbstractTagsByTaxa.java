/*
 * AbstractTagsByTaxa
 */
package net.maizegenetics.gbs.tagdist;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.nio.channels.FileChannel;
import java.nio.ByteBuffer;
import java.nio.CharBuffer;

import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.util.BitUtil;

/**
 * Basic method implementation of Tags by Taxa, including methods for reading and writing
 * files
 * @author ed, Gabriel Rodrigues Alves Margarido
 */
public abstract class AbstractTagsByTaxa extends AbstractTags implements TagsByTaxa {

    protected int taxaNum = 0;
    protected String[] taxaNames;

    @Override
    public int getTaxaCount() {
        return taxaNames.length;
    }

    @Override
    public String getTaxaName(int taxaIndex) {
        return taxaNames[taxaIndex];
    }

    @Override
    public String[] getTaxaNames() {
        return taxaNames;
    }

    @Override
    public int getIndexOfTaxaName(String taxon) {
        for (int i = 0; i < taxaNames.length; i++) {
            if (taxon.equals(taxaNames[i])) {
                return i;
            }
        }
        return -1;
    }

    /**Truncates each taxon name to everything up to the first colon, and converts
     * all letters to lower case.  This causes tag counts from taxa with identical
     * names to be merged. */
    public void truncateTaxonNames() {
        for (int i = 0; i < taxaNames.length; i++) {
            taxaNames[i] = taxaNames[i].substring(0, taxaNames[i].indexOf(":")).toLowerCase();
        }
    }

    synchronized public void addReadsToTagTaxon(int tagIndex, int taxaIndex, int addValue) {
        setReadCountForTagTaxon(tagIndex, taxaIndex, addValue + getReadCountForTagTaxon(tagIndex, taxaIndex));
    }

    public short[] getTaxaReadCountsForTag(int tagIndex) {
        short[] r = new short[getTaxaCount()];
        for (int i = 0; i < getTaxaCount(); i++) {
            r[i] = (short) getReadCountForTagTaxon(tagIndex, i);
        }
        return r;
    }

    /**@param tagIndex Number of organisms/bits that must be stored in the bitset.
     *@return r A new OpenBitSet object containing the tag distribution of its parent TagsByTaxa object.*/
    public OpenBitSet getTaxaReadBitsForTag(int tagIndex) {
        OpenBitSet r = new OpenBitSet(getTaxaCount());
        for (int i = 0; i < getTaxaCount(); i++) {
            if (getReadCountForTagTaxon(tagIndex, i) > 0) {
                r.fastSet(i);
            }
        }
        return r;
    }

    public int getNumberOfTaxaWithTag(int readIndex) {   // how many taxa was a given read seen in?
        int nTaxaWData = 0;
        for (int i = 0; i < getTaxaCount(); i++) {
            if (getReadCountForTagTaxon(readIndex, i) > 0) {
                nTaxaWData++;
            }
        }
        return nTaxaWData;
    }

    public void readDistFile(File inFile, FilePacking numberType) {
        System.out.println("Reading Haplotypes distribution from:" + inFile.toString());
        switch (numberType) {
            case Text:
                readTextDistFile(inFile);
                break;
            case Bit:
                readBitDistFile(inFile);
                break;
            default:
                readByteShortDistFile(inFile, numberType);
                break;
        }
    }

    public void readChannel(File inFile) {
        try {
            int hapsOutput = 0;
            ByteBuffer bb = ByteBuffer.allocate(1024);
            FileChannel inputChannel = new FileInputStream(inFile).getChannel();
            while (inputChannel.read(bb) != -1) {
                CharBuffer inputBuffer = bb.asCharBuffer();
                inputBuffer.flip();
                while (inputBuffer.hasRemaining()) {
                    System.out.println(inputBuffer.get());
                }
                inputBuffer.clear();
            }
            inputChannel.close();
        } catch (Exception e) {
            System.out.println("Catch in reading input channel: " + e);
            e.printStackTrace();
        }
    }

    void readTextDistFile(File inFile) {
        int hapsOutput = 0;
        String[] inputLine;
        try {
            BufferedReader br = new BufferedReader(new FileReader(inFile), 65536);
            inputLine = br.readLine().split("\t");
            int tagNum = Integer.parseInt(inputLine[0]);
            tagLengthInLong = Integer.parseInt(inputLine[1]);
            taxaNum = Integer.parseInt(inputLine[2]);
            initMatrices(taxaNum, tagNum);
            inputLine = br.readLine().trim().split("\t");
            for (int t = 0; t < taxaNum; t++) {
                taxaNames[t] = inputLine[t];
            }
            for (int i = 0; i < tagNum; i++) {
                inputLine = br.readLine().split("\t");
                long[] tt = BaseEncoder.getLongArrayFromSeq(inputLine[0]);
                for (int j = 0; j < tt.length; j++) {
                    tags[j][i] = tt[j];
                }
                tagLength[i] = Byte.valueOf(inputLine[1]);
                for (int t = 0; t < taxaNum; t++) {
                    setReadCountForTagTaxon(i, t, Byte.valueOf(inputLine[t + 2]));
                }
                hapsOutput++;
            }
        } catch (Exception e) {
            System.out.println("Catch in writing output file e=" + e);
            e.printStackTrace();
        }
        System.out.println("Number of Taxa in file:" + taxaNum);
        System.out.println("Number of Haplotypes in file:" + hapsOutput);
    }

    void readByteShortDistFile(File inFile, FilePacking countType) {
        int hapsOutput = 0;
        try {
            DataInputStream rw = new DataInputStream(new BufferedInputStream(new FileInputStream(inFile), 4000000));
            int tagNum = rw.readInt();
            tagLengthInLong = rw.readInt();
            taxaNum = rw.readInt();
            initMatrices(taxaNum, tagNum);
            for (int t = 0; t < taxaNum; t++) {
                taxaNames[t] = rw.readUTF();
            }
            for (int i = 0; i < tagNum; i++) {
                for (int j = 0; j < tagLengthInLong; j++) {
                    tags[j][i] = rw.readLong();
                }
                tagLength[i] = rw.readByte();
                for (int t = 0; t < taxaNum; t++) {
                    if (countType == FilePacking.Short) {
                        setReadCountForTagTaxon(i, t, rw.readShort());
                    } else {
                        setReadCountForTagTaxon(i, t, rw.readByte());
                    }
                }
                hapsOutput++;
            }
            rw.close();
        } catch (Exception e) {
            System.out.println("Catch in writing output file e=" + e);
        }
        System.out.println("Number of Taxa in file:" + taxaNum);
        System.out.println("Number of Haplotypes in file:" + hapsOutput);
    }

    void readBitDistFile(File inFile) {
        int hapsOutput = 0;
        try {
            DataInputStream rw = new DataInputStream(new BufferedInputStream(new FileInputStream(inFile), 4000000));
            int tagNum = rw.readInt();
//                    System.out.println(tagNum);
            tagLengthInLong = rw.readInt();
//                    System.out.println(tagLengthInLong);
            taxaNum = rw.readInt();
//                    System.out.println(taxaNum);
            initMatrices(taxaNum, tagNum);
            for (int t = 0; t < taxaNum; t++) {
                taxaNames[t] = rw.readUTF();
//                    System.out.println(taxaNames[t]);
            }
            int numberOfLongs = BitUtil.bits2words(taxaNum);
            long[] distInLong = new long[numberOfLongs];
            OpenBitSet obs;
            for (int i = 0; i < tagNum; i++) {
                for (int j = 0; j < tagLengthInLong; j++) {
                    tags[j][i] = rw.readLong();
                }
                tagLength[i] = rw.readByte();
                for (int j = 0; j < numberOfLongs; j++) {
                    distInLong[j] = rw.readLong();
                }
                obs = new OpenBitSet(distInLong, taxaNum);
                for (int t = 0; t < taxaNum; t++) {
                    if (obs.fastGet(t)) {
                        setReadCountForTagTaxon(i, t, 1);
//                        System.out.print("1\t");
                    } else {
//                        System.out.print("0\t");
                    }
                }

                hapsOutput++;
            }
            rw.close();
        } catch (Exception e) {
            System.out.println("Catch in reading input line " + hapsOutput + ": " + e);
            e.printStackTrace();
        }
        System.out.println("Number of Taxa in file:" + taxaNum);
        System.out.println("Number of Haplotypes in file:" + hapsOutput);
    }

    public void writeDistFile(File outFile, FilePacking numberType, int minCount) {
        int hapsOutput = 0;
        int outReads = minCount > 0 ? readsWCountsGreaterThanMin(minCount) : getTagCount();
        System.out.println(outReads + " tags will be output to " + outFile.getName());
        try {
            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile), 65536));
            switch (numberType) {
                case Text:
                    hapsOutput = writeTextDistFile(fw, outReads, minCount);
                    break;
                case Bit:
                    hapsOutput = writeBitDistFile(fw, outReads, minCount);
                    break;
                default:
                    hapsOutput = writeByteShortDistFile(fw, numberType, outReads, minCount);
                    break;
            }
            fw.flush();
            fw.close();
            System.out.println("Tags written to:" + outFile.toString());
            System.out.println("Number of tags in file:" + hapsOutput);
        } catch (Exception e) {
            System.out.println("Catch in writing output file e=" + e);
            e.printStackTrace();
        }
    }

    private int writeTextDistFile(DataOutputStream fw, int outReads, int minCount) {
        int hapsOutput = 0;
        try {
            fw.writeBytes(outReads + "\t" + tagLengthInLong + "\t" + taxaNum + "\n");
            for (int t = 0; t < taxaNum; t++) {
                fw.writeBytes(taxaNames[t] + "\t");
            }
            fw.writeBytes("\n");

            for (int i = 0; i < tags[0].length; i++) {
                if (minCount > 0 && getReadCount(i) < minCount) {
                    continue;
                }
                fw.writeBytes(BaseEncoder.getSequenceFromLong(getTag(i)) + "\t");
                fw.writeBytes(getTagLength(i) + "\t");
                for (int t = 0; t < taxaNum; t++) {
                    fw.writeBytes(getReadCountForTagTaxon(i, t) + "\t");
                }
                fw.writeBytes("\n");
                hapsOutput++;
            }

        } catch (Exception e) {
            System.out.println("Catch in writeTextDistFile writing output file e=" + e);
            e.printStackTrace();
        }
        return hapsOutput;
    }

    private int writeByteShortDistFile(DataOutputStream fw, FilePacking countType, int outReads, int minCount) {
        int hapsOutput = 0;
        try {
            fw.writeInt(outReads);
            fw.writeInt(tagLengthInLong);
            fw.writeInt(taxaNum);
            for (int t = 0; t < taxaNum; t++) {
                fw.writeUTF(taxaNames[t]);
            }

            for (int i = 0; i < tags[0].length; i++) {
                if (minCount > 0 && getReadCount(i) < minCount) {
                    continue;
                }
                for (int j = 0; j < tagLengthInLong; j++) {
                    fw.writeLong(tags[j][i]);
                }
                fw.writeByte(tagLength[i]);
                for (int t = 0; t < taxaNum; t++) {
                    if (countType == FilePacking.Short) {
                        fw.writeShort(getReadCountForTagTaxon(i, t));
                    } else {
                        fw.writeByte(getReadCountForTagTaxon(i, t));
                    }
                }
                hapsOutput++;
            }

        } catch (Exception e) {
            System.out.println("Catch in writeTextDistFile writing output file e=" + e);
            e.printStackTrace();
        }
        return hapsOutput;
    }

    private int writeBitDistFile(DataOutputStream fw, int outReads, int minCount) {
        int hapsOutput = 0;
        OpenBitSet obs = new OpenBitSet(taxaNum);
        try {
            fw.writeInt(outReads);
            fw.writeInt(tagLengthInLong);
            fw.writeInt(taxaNum);
            for (int t = 0; t < taxaNum; t++) {
                fw.writeUTF(taxaNames[t]);
            }

            for (int i = 0; i < tags[0].length; i++) {
                if (minCount > 0 && getReadCount(i) < minCount) {
                    continue;  //Skip tags whose count is lower than a cutoff
                }
                for (int j = 0; j < tagLengthInLong; j++) {
                    fw.writeLong(tags[j][i]);
                }
                fw.writeByte(tagLength[i]);
                obs = new OpenBitSet(taxaNum);
                for (int t = 0; t < taxaNum; t++) {
                    if (getReadCountForTagTaxon(i, t) > 0) {
                        obs.set(t);
                    }
                }
                long[] obsInLong = obs.getBits();
                for (int t = 0; t < obsInLong.length; t++) {
                    fw.writeLong(obsInLong[t]);
                }
                hapsOutput++;
            }
        } catch (Exception e) {
            System.out.println("Catch in writeTextDistFile writing output file e=" + e);
            e.printStackTrace();
        }
        return hapsOutput;
    }

    /** Writes a tag count file in binary (true) or text (false) format, containing the tag sequences from this TBT file.*/
    public void writeReadCountFile(File outFile, boolean binary) {
        int hapsOutput = 0;
        try {
            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile), 4000000));
            if (binary) {
                fw.writeInt(tags[0].length);
                fw.writeInt(tagLengthInLong);
            } else {
                fw.writeBytes(getTagCount() + "\t" + tagLengthInLong + "\n");
            }
            for (int i = 0; i < tags[0].length; i++) {
                if (!binary) {
                    fw.writeBytes(BaseEncoder.getSequenceFromLong(getTag(i)) + '\t' + tags[0].length + '\t' + getReadCount(i) + "\n");
                } else {
                    for (int j = 0; j < tagLengthInLong; j++) {
                        fw.writeLong(tags[j][i]);
                    }
                    fw.writeByte(getTagLength(i));
                    fw.writeInt(getReadCount(i));
                }
                hapsOutput++;
            }
            fw.flush();
            fw.close();
            System.out.println("Reads written to:" + outFile.toString());
            System.out.println("Number of Reads in file:" + hapsOutput);
        } catch (Exception e) {
            System.out.println("Catch in writing output file e=" + e);
        }
    }

    private int readsWCountsGreaterThanMin(int minCount) {
        if (minCount < 1) {
            return getTagCount();
        }
        int sum = 0;
        for (int i = 0; i < getTagCount(); i++) {
            if (getReadCount(i) >= minCount) {
                sum++;
            }
        }
        return sum;
    }

    public int getReadCount(int tagIndex) {
        int sum = 0;
        for (int i = 0; i < getTaxaCount(); i++) {
            sum += getReadCountForTagTaxon(tagIndex, i);
        }
        return sum;
    }

    @Override
    public boolean areTagsUnique() {
        return true;
    }
}
