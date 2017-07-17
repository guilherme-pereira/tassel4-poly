/*
 * AbstractTags
 */
package net.maizegenetics.gbs.tagdist;

import cern.colt.GenericSorting;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Arrays;
import net.maizegenetics.gbs.util.BaseEncoder;
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.NucleotideCompound;

/**
 * Basic methods for working with PE Tags, including sorting and search.
 * @author Fei Lu
 */
public abstract class AbstractPETags implements PETags {
    /**Tag length in Long primitive data type*/
    protected int tagLengthInLong; 
    /**
     * Forward tags in Long array
     * Columns are sequence, rows are index of sequence. This is faster than the other way
     */
    protected long[][] tagsF;
    /**
     * Backward tags in Long array
     * Columns are sequence, rows are index of sequence. This is faster than the other way
     */
    protected long[][] tagsB;
    /**Length of forward tags*/
    protected short[] tagFLength;  // length of tag (number of bases)
    /**Length of backward tags*/
    protected short[] tagBLength;
    /**
     * Length of contigs in Long primitive data type, 0 means contig doesn't exist
     */
    protected byte[] contigLengthInLong;
    /** Contigs in Long array*/
    protected long[][] contig;
    /**Length of contigs, 0 means contig doesn't exist*/
    protected short[] contigLength;  //entire length of PE tags (cutsite to cutsite), 0 means contig doesn't exist
     
    @Override
    public long[] getTagF(int index) {
        return tagsF[index];
    }
    
    @Override
    public long[] getTagB(int index) {
        return tagsB[index];
    }

    @Override
    public int getTagCount() {
        return tagsF.length;
    }
    
    @Override
    public int getTagIndex(long[] tagF, long[] tagB) {
        //code inspired by COLT lower bound function
        int first = 0;
        int len = this.getTagCount() - first;
        int comp = 0;
        while (len > 0) {
            int half = len / 2;
            int middle = first + half;
            if ((comp = compareTags(middle, tagF, tagB)) < 0) {
                first = middle + 1;
                len -= half + 1;
            } else {
                len = half;
            }
        }
        if ((first < this.getTagCount()) && (compareTags(first, tagF, tagB) == 0)) {
            return first;
        }
        return -(first + 1);
    }

    @Override
    public long[] getContig (int index) {
        if (this.getContigLengthInLong(index) == 0) return null;
        return contig[index];
    }
    
    @Override
    public int getContigCount () {
        int sum = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.getContigLengthInLong(i) == 0) continue;
            sum++;
        }
        return sum;
    }
    
    @Override
    public short getContigLength (int index) {
        return contigLength[index];
    }
    
    @Override
    public byte getContigLengthInLong (int index) {
        return contigLengthInLong[index];
    }
    
    /**
     * Write forword and backword fasta files of PE tag for Bowtie2 alignment, the sequence is in format of index_f/b_c/n, where f means forward, b means backword, c means contig, n means no contig
     * When contig exist, the forward and backword tags are the contig and reverse complement of the contig, respectively.
     * @param fFastaFileS
     * @param bFastaFileS 
     */
    public void mkFastaFile (String fFastaFileS, String bFastaFileS) {
        try {
            BufferedWriter bwf = new BufferedWriter (new FileWriter (fFastaFileS), 65536);
            BufferedWriter bwb = new BufferedWriter (new FileWriter (bFastaFileS), 65536);
            for (int i = 0; i < this.getTagCount(); i++) {
                String forwardName = String.valueOf(i)+"_f_";
                String backwardName = String.valueOf(i)+"_b_";
                String forwardTag;
                String backwardTag;
                if (this.getContigLength(i) == 0) {
                    forwardName = forwardName+"n";
                    backwardName = backwardName+"n";
                    forwardTag = BaseEncoder.getSequenceFromLong(this.getTagF(i)).substring(0, this.getTagFLength(i));
                    backwardTag = BaseEncoder.getSequenceFromLong(this.getTagB(i)).substring(0, this.getTagBLength(i));
                }
                else {
                    forwardName = forwardName+"c";
                    backwardName = backwardName+"c";
                    forwardTag = BaseEncoder.getSequenceFromLong(this.getContig(i)).substring(0, this.getContigLength(i));
                    backwardTag = BaseEncoder.getReverseComplement(forwardTag);
                }
                
                bwf.write(">"+forwardName);
                bwf.newLine();
                bwf.write(forwardTag);
                bwf.newLine();
                bwb.write(">"+backwardName);
                bwb.newLine();
                bwb.write(backwardTag);
                bwb.newLine();
            }
            bwf.flush();
            bwf.close();
            bwb.flush();
            bwb.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        } 
    }
    
    /**
     * Write fasta files of PE tag for alignment
     * @param tagsFFastaFileS
     * @param tagsBFastaFileS
     * @param contigFastaFileS
     * @deprecated 
     */
    public void mkFastaFile (String tagsFFastaFileS, String tagsBFastaFileS, String contigFastaFileS) {
        try {
            BufferedWriter bwf = new BufferedWriter (new FileWriter (tagsFFastaFileS), 65536);
            BufferedWriter bwb = new BufferedWriter (new FileWriter (tagsBFastaFileS), 65536);
            BufferedWriter bwc = new BufferedWriter (new FileWriter (contigFastaFileS), 65536);
            for (int i = 0; i < this.getTagCount(); i++) {
                bwf.write(">"+String.valueOf(i));
                bwf.newLine();
                bwf.write(BaseEncoder.getSequenceFromLong(this.getTagF(i)).substring(0, this.getTagFLength(i)));
                bwf.newLine();
                bwb.write(">"+String.valueOf(i));
                bwb.newLine();
                bwb.write(BaseEncoder.getSequenceFromLong(this.getTagB(i)).substring(0, this.getTagBLength(i)));
                bwb.newLine();
                if (this.getContigLength(i) == 0) continue;
                bwc.write(">"+String.valueOf(i));
                bwc.newLine();
                bwc.write(BaseEncoder.getSequenceFromLong(this.getContig(i)).substring(0, this.getContigLength(i)));
                bwc.newLine();
            }
            bwf.flush();
            bwf.close();
            bwb.flush();
            bwb.close();
            bwc.flush();
            bwc.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }  
    }
    
    /**
     * Order forward and backward tags
     */
    public void orderTagsFTagsB () {
        for (int i = 0; i < this.getTagCount(); i++) {
            this.orderTagFTagB(i);
        }
    }
    
    /**
     * Order forward and backward tags for a particular PE tag
     * @param index
     * @return boolean value of that if the forward and backward is switched.
     */
    private boolean orderTagFTagB (int index) {
        if (this.compareTagFTagB(index) == 1) {
            this.switchTagFTagB(index);
            return true;
        }
        return false;
    }
    
    /**
     * Switch forward and backward tags
     * @param index 
     */
    private void switchTagFTagB (int index) {
        long[] temp = tagsF[index];
        tagsF[index] = tagsB[index];
        tagsB[index] = temp;
        short tl = tagFLength[index];
        tagFLength[index] = tagBLength[index];
        tagBLength[index] = tl;
    }
    
    /**
     * Make contigs from forward and backward tags
     */
    public void contigPETags () {
        SimpleGapPenalty gapPen = new SimpleGapPenalty((short) 10, (short) 10);
        SubstitutionMatrix<NucleotideCompound> subMatrix = SubstitutionMatrixHelper.getNuc4_4();
        int minOverlap = 30;
        double minIden = 0.97;
        int halfLength = 128;
        for (int i = 0; i < this.getTagCount(); i++) {
            String queryS = BaseEncoder.getSequenceFromLong(this.getTagF(i)).substring(0, this.getTagFLength(i));
            DNASequence query = new DNASequence(queryS);
            String hitS = BaseEncoder.getReverseComplement(BaseEncoder.getSequenceFromLong(this.getTagB(i)).substring(0, this.getTagBLength(i)));            
            DNASequence hit = new DNASequence(hitS);
            SequencePair<DNASequence, NucleotideCompound> psa;
            psa = Alignments.getPairwiseAlignment(query, hit, Alignments.PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
            int queryStart = psa.getIndexInQueryAt(1);
            int queryEnd = psa.getIndexInQueryAt(psa.getLength());
            int hitStart = psa.getIndexInTargetAt(1);
            int hitEnd = psa.getIndexInTargetAt(psa.getLength());
            int overlap = psa.getLength();
            int idenNum = psa.getNumIdenticals();
            
            if (i%10000 == 0) {
                System.out.println("Contigged " + i + " PETags. Total: " + this.getTagCount());
            }
            
            if (hitStart != 1) continue;
            if (queryEnd != this.getTagFLength(i)) continue;
            if (overlap < minOverlap) continue;
            if ((double)(idenNum/overlap) < minIden) continue;
            StringBuilder sb = new StringBuilder();
            sb.append(queryS.substring(0, queryStart-1));
            if (queryEnd < halfLength) {
                sb.append(queryS.substring(queryStart-1, queryEnd));
                int startPos = queryEnd - queryStart + hitStart;
                sb.append(hitS.substring(startPos, hitS.length()));
            }
            else {
                if (queryStart < halfLength) {
                    sb.append(queryS.substring(queryStart-1, halfLength));
                    int startPos = halfLength-queryStart+hitStart;
                    sb.append(hitS.substring(startPos, hitS.length()));
                }
                else {
                    sb.append(hitS.substring(hitStart-1, hitS.length()));
                }
            }
//            String contigS = sb.toString();
//            System.out.println(i+"\t"+queryStart+"\t"+queryEnd+"\t"+this.getTagFLength(i)+"\t"+hitStart+"\t"+hitEnd+"\t"+this.getTagBLength(i)+"\t"+psa.getLength());
//            System.out.println(psa.toString(300));
//            System.out.println(contigS);
//            System.out.println(queryS);
//            System.out.println(hitS);
//            System.out.println("\n\n");
            contigLength[i] = (short)sb.length();
            int leftover = contigLength[i] % BaseEncoder.chunkSize;
            if (leftover == 0) {
                contigLengthInLong[i] = (byte)(contigLength[i] / BaseEncoder.chunkSize);
            }
            else {
                contigLengthInLong[i] = (byte)(contigLength[i] / BaseEncoder.chunkSize + 1);
                for (int j = 0; j < BaseEncoder.chunkSize-leftover; j++) {
                    sb.append("A");
                }
            }
            String contigS = sb.toString();
            long[] temp = BaseEncoder.getLongArrayFromSeq(contigS);
            contig[i] = temp;
        }
    }
    
    /**
     * Compare forward tag and backward tag
     * @param index
     * @return result of comparison from forward tag and backward tag -1(<), 0(=), 1(>) 
     */
    private int compareTagFTagB (int index) {
        for (int i = 0; i < tagLengthInLong; i++) {
            if (tagsF[index][i] < tagsB[index][i]) return -1;
            if (tagsF[index][i] > tagsB[index][i]) return 1;
        }
        return 0;
    }
    
    /**
     * Compare PE tag with another PE tag
     * @param index
     * @param tagF
     * @param tagB
     * @return result of comparison from different PE tags -1(<), 0(=), 1(>) 
     */
    private int compareTags(int index, long[] tagF, long[] tagB) {
        for (int i = 0; i < tagLengthInLong; i++) {
            if (tagsF[index][i] < tagF[i]) {
                return -1;
            }
            if (tagsF[index][i] > tagF[i]) {
                return 1;
            }
        }
        for (int i = 0; i < tagLengthInLong; i++) {
            if (tagsB[index][i] < tagB[i]) {
                return -1;
            }
            if (tagsB[index][i] > tagB[i]) {
                return 1;
            }
        }
        return 0;
    }

    
    @Override
    public short getTagFLength(int index) {
        return tagFLength[index];
    }
    
    @Override
    public short getTagBLength(int index) {
        return tagBLength[index];
    }

    @Override
    public int getTagSizeInLong() {
        return tagLengthInLong;
    }

    /**
     * Initialize the matrix of data structure extending AbstractPETags
     * @param tagLengthInLong
     * @param tagNum 
     */
    protected void iniMatrix (int tagLengthInLong, int tagNum) {
        this.tagLengthInLong = tagLengthInLong;
        tagsF = new long[tagNum][tagLengthInLong];
        tagsB = new long[tagNum][tagLengthInLong];
        tagFLength = new short[tagNum];
        tagBLength = new short[tagNum];
        contigLengthInLong = new byte[tagNum];
        contig = new long[tagNum][];
        contigLength = new short[tagNum];
    }
    
    /**
     * Swap two PE tags
     * @param index1
     * @param index2
     */
    @Override
    public void swap(int index1, int index2) {
        long[] temp;
        temp = tagsF[index1];
        tagsF[index1] = tagsF[index2];
        tagsF[index2] = temp;
        temp = tagsB[index1];
        tagsB[index1] = tagsB[index2];
        tagsB[index2] = temp;
        short tl;
        tl = tagFLength[index1];
        tagFLength[index1] = tagFLength[index2];
        tagFLength[index2] = tl;
        tl = tagBLength[index1];
        tagBLength[index1] = tagBLength[index2];
        tagBLength[index2] = tl;
        byte temByte = contigLengthInLong[index1];
        contigLengthInLong[index1] = contigLengthInLong[index2];
        contigLengthInLong[index2] = temByte;
        temp = contig[index1];
        contig[index1] = contig[index2];
        contig[index2] = temp;
        tl = contigLength[index1];
        contigLength[index1] = contigLength[index2];
        contigLength[index2] = tl;
    }
    
    /**
     * Compare two PE tags
     * @param index1
     *        index of the first PE tag
     * @param index2
     *        index of the second PE tag
     * @return result of comparison from two PE tags -1(<), 0(=), 1(>)
     */
    @Override
    public int compare(int index1, int index2) {
        for (int i = 0; i < tagLengthInLong; i++) {
            if (tagsF[index1][i] < tagsF[index2][i]) {
                return -1;
            }
            if (tagsF[index1][i] > tagsF[index2][i]) {
                return 1;
            }
        }
        for (int i = 0; i < tagLengthInLong; i++) {
            if (tagsB[index1][i] < tagsB[index2][i]) {
                return -1;
            }
            if (tagsB[index1][i] > tagsB[index2][i]) {
                return 1;
            }
        }
        return 0;
    }

    /**
     *  sort PE tags by Long primitive data type value
     */
    public void sort() {
        System.out.println("Position index sort begin.");
        GenericSorting.quickSort(0, this.getTagCount(), this, this);
        System.out.println("Position index sort end.");
    }
}
