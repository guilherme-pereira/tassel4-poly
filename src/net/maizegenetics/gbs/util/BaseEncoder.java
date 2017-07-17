package net.maizegenetics.gbs.util;

/**
 * Utility class for encoding tags into longs.
 * <p>
 * Sequencing reads are chunked into 32bp and recorded in a 64-bit long.  Only 
 * A (00), C (01), G (10), T (11) are encoded.  Any other character sets the entire long to -1.
 * Missing data at the end is padded with poly-A or (0).  This missing end, is tracked
 * by the tag length attribute.
 * <p>
 * Some of these methods should be transitioned to {@link net.maizegenetics.pal.alignment.NucleotideAlignmentConstants},
 * however, BaseEncoder only supports four states, while NucleotideAlignment includes gaps, insertions, and missing.
 * 
 * @author Ed Buckler
 */
public class BaseEncoder {

    /** defines the number of bases fitting with a long */
    public static final int chunkSize = 32;
    public static final int chunkSizeForInt = 16;
    /** defines the base order */
    public static final char[] bases = {'A', 'C', 'G', 'T'};

    private BaseEncoder() {
    }

    /**
     * Returns a long for a sequence in a String
     * @param seq
     * @return 2-bit encode sequence (-1 if an invalid sequence state is provided e.g. N)
     */
    public static long getLongFromSeq(String seq) {
        int seqLength = seq.length();
        long v = 0;
        for (int i = 0; i < seqLength; i++) {
            switch (seq.charAt(i)) {
                case 'A':
                case 'a':
                    v = v << 2;
                    break;
                case 'C':
                case 'c':
                    v = (v << 2) + (byte) 1;
                    break;
                case 'G':
                case 'g':
                    v = (v << 2) + (byte) 2;
                    break;
                case 'T':
                case 't':
                    v = (v << 2) + (byte) 3;
                    break;
                default:
                    return -1;
            }
        }
        if (seqLength == chunkSize) {
            return v;
        }
        if (seqLength > chunkSize) {
            return -1;
        }
        v = (v << (2 * (chunkSize - seqLength))); //if shorter fill with AAAA
        return v;
    }

    /**
     * @param seq A String containing a DNA sequence.
     * @return result A array of Long containing the binary representation of the sequence.
     * null if sequence length is not a multiple of BaseEncoder.chunksize.
     */
    public static long[] getLongArrayFromSeq(String seq) {
        if (seq.length() % chunkSize != 0) {
            return null;
        }
        long[] result = new long[seq.length() / chunkSize];
        for (int i = 0; i < result.length; i++) {
            result[i] = getLongFromSeq(seq.substring(i * chunkSize, (i + 1) * chunkSize));
        }
        return result;
    }
    
    /**
     * @param seq A String containing a DNA sequence.
     * @return result A array of Long containing the binary representation of the sequence.
     * if sequence length is shorter than padded Length adds A to the end.
     */
    public static long[] getLongArrayFromSeq(String seq, int paddedLength) {
        if(seq.length()<paddedLength) {
            seq=seq+"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".substring(0, paddedLength-seq.length());
        }
    //    System.out.println("PadLength:"+seq.length());
        return getLongArrayFromSeq(seq);
    }



    /**
     * Returns the reverse complement of a sequence already encoded in a 2-bit long.
     * <p>
     * Note: polyA is used represent unknown, but reverse complement will change it to polyT which does not mean the same
     * sometimes it is best to reverseComplement by text below
     * @param seq  2-bit encoded sequence
     * @param len  length of the sequence
     * @return  2-bit reverse complement
     */
    public static long getReverseComplement(long seq, byte len) {
        // if(seq==-1) return -1;
        long rev = 0;
        // byte b=0;
        long mask = 3;
        seq = ~seq;
        for (int i = 0; i < len; i++) {
            rev = (rev << 2) + (seq & mask);
            seq = seq >> 2;
            // System.out.println("v = " + v);
        }
        return rev;
    }

    /**
     * Returns the reverse complement of a sequence already encoded in a 2-bit long.
     * The entire long (32-bp) is reverse complemented.
     * <p>
     * Note: polyA is used represent unknown, but reverse complement will change it to polyT which does not mean the same
     * sometimes it is best to reverseComplement by text below
     * @param seq  2-bit encoded sequence
     * @return  2-bit reverse complement
     */
    public static long getReverseComplement(long seq) {
        return getReverseComplement(seq, (byte) chunkSize);
    }


    /**
     * Returns the reverse complement of a arrays of sequences already encoded in a 2-bit long.
     * <p>
     * Note: polyA is used represent unknown, but reverse complement will change it to polyT which does not mean the same
     * sometimes it is best to reverseComplement by text below
     * @param seq  array of 2-bit encoded sequences
     * @return  array of 2-bit reverse complements
     */
    public static long[] getReverseComplement(long[] seq) {
        long[] rev = new long[seq.length];
        for (int i = 0; i < rev.length; i++) {
            rev[i] = getReverseComplement(seq[seq.length - i - 1], (byte) chunkSize);
        }
        return rev;
    }

    /**
     * Returns a string based reverse complement.  Get around issues with the poly-A tailing in the 2-bit encoding approach.
     *
     * @param seq  DNA sequence
     * @return  reverse complement DNA sequence
     */
    public static String getReverseComplement(String seq) {
        StringBuilder sb = new StringBuilder(seq.length());
        for (int i = seq.length() - 1; i >= 0; i--) {
            sb.append(getComplementBase(seq.charAt(i)));
        }
        return sb.toString();
    }

    /**
     * Returns reverse complement for a sequence.
     * @param base
     * @return  reverse complement of base
     */
    public static char getComplementBase(char base) {
        switch (base) {
            case 'A':
                return 'T';
            case 'C':
                return 'G';
            case 'G':
                return 'C';
            case 'T':
                return 'A';
        }
        return 'N';
    }

    /**
     * Returns the byte {@link net.maizegenetics.pal.alignment.NucleotideAlignmentConstants} representation
     * used by TASSEL for the 2-bit encoded long.
     * <p>
     * e.g. A > 2-bit encode 00 > byte (0)
     * @param val 2-bit encoded DNA sequence
     * @return array of bytes for the DNA sequence
     */
    public static byte[] getByteSeqFromLong(long val) {
        byte[] b = new byte[chunkSize];
        long mask = 3;
        for (int i = 0; i < chunkSize; i++) {
            b[chunkSize - i - 1] = (byte) (val & mask);
            val = val >> 2;
        }
        return b;
    }

    /**
     * Returns the byte {@link net.maizegenetics.pal.alignment.NucleotideAlignmentConstants} representation
     * used by TASSEL for the 2-bit encoded long.
     * <p>
     * e.g. A > 2-bit encode 00 > byte (0)
     * @param valA array of 2-bit encoded DNA sequence
     * @return array of bytes for the DNA sequence
     */
    public static byte[] getByteSeqFromLong(long[] valA) {
        byte[] b = new byte[chunkSize * valA.length];
        long mask = 3;
        long val;
        for (int j = 0; j < valA.length; j++) {
            val = valA[j];
            for (int i = 0; i < chunkSize; i++) {
                b[(j * chunkSize) + chunkSize - i - 1] = (byte) (val & mask);
                val = val >> 2;
            }
        }
        return b;
    }

    /**
     * Returns the 2-bit encoded long represented by 32 bytes representing {@link net.maizegenetics.pal.alignment.NucleotideAlignmentConstants} 
     * representation
     * <p>
     * @param b array of bytes encoding NucleotideAlignmentConstants
     * @return 2-bit encoded long
     */
    public static long getLongSeqFromByteArray(byte[] b) {
        //the byte array must be in 0-3 coding for A, C, G, T
        long v = 0;
        if (b.length != chunkSize) {
            return -1;
        }
        for (int i = 0; i < b.length; i++) {
            v = (v << 2) + b[i];
        }
        return v;
    }
    
     /**
     * Return a string representation of the 2-bit encoded long.
     * @param val 2-bit encoded sequence
     * @param len length of the sequence
     * @return DNA sequence as a string
     */
    public static String getSequenceFromLong(long val, byte len) {
        StringBuilder seq = new StringBuilder(chunkSize + 4);
        long mask = 3;
        for (int i = 0; i < len; i++) {
            byte base = (byte) (val & mask);
            seq.insert(0, bases[base]);
            val = val >> 2;
        }
        return seq.toString();
    }

     /**
     * Return a string representation of an array of 2-bit encoded longs.
     * @param val array of 2-bit encoded sequences
     * @return DNA sequence as a string
     */
    public static String getSequenceFromLong(long[] val) {
        StringBuilder seq = new StringBuilder();
        for (long v : val) {
            seq.append(getSequenceFromLong(v));
        }
        return seq.toString();
    }

    /**
     * Split a 2-bit encoded long into 2 integers.
     * @param val 2-bit encoded long sequence
     * @return array of 2-bit encoded integers
     */
    public static int[] getIntFromLong(long val) {
        int[] ival = new int[2];
        ival[0] = (int) (val >> chunkSize);
        ival[1] = (int) (val);
        return ival;
    }

    /**
     * Return a string representation of the 2-bit encoded Integer (16bp).
     * @param val 2-bit encoded sequence
     * @return DNA sequence as a string
     */
    public static String getSequenceFromInt(int val) {
        StringBuilder seq = new StringBuilder(chunkSizeForInt + 1);
        long mask = 3;
        for (int i = 0; i < chunkSizeForInt; i++) {
            byte base = (byte) (val & mask);
            seq.insert(0, bases[base]);
            val = val >> 2;
        }
        return seq.toString();
    }

    /**
     * Returns the position of the first low quality positions based on a quality
     * fastq (?) string.
     * @param quality fastq quality string
     * @param minQual minimum quality threshold
     * @return position of first low quality position (quality length is returned is not low 
     * quality base is found.
     */
    public static int getFirstLowQualityPos(String quality, int minQual) {
        int qualInt = 0;
        for (int i = 0; i < quality.length(); i++) {
            qualInt = (int) quality.charAt(i) - 64;
            if (qualInt < minQual) {
                return i;
            }
        }
        return quality.length();
    }


    /**
     * Return a string representation of the 2-bit encoded long.
     * @param val 2-bit encoded sequence
     * @return DNA sequence as a string
     */
    public static String getSequenceFromLong(long val) {
        return getSequenceFromLong(val, (byte) chunkSize);
    }

    /**
     * Returns the number of bp differences between two 2-bit encoded longs.
     * Maximum divergence is used to save time when only interested in very similar 
     * sequences.
     * @param seq1 2-bit encoded sequence
     * @param seq2 2-bit encoded sequence
     * @param maxDivergence threshold for counting divergence upto
     * @return count of the divergence (above the maxDivergence, chunkSize is returned)
     */
    public static byte seqDifferences(long seq1, long seq2, int maxDivergence) {
        long mask = 3;
        byte cnt = 0;
        long diff = seq1 ^ seq2;
        for (int x = 0; x < chunkSize && cnt <= maxDivergence; x++) {
            if ((diff & mask) > 0) {
                cnt++;
            }
            diff = diff >> 2;
            // System.out.println("v = " + v);
        }
        if (cnt > maxDivergence) {
            cnt = (byte) chunkSize;
        }
        // if(x<(chunkSize-1)) cnt=(byte)chunkSize;  //if didn't get to the end of the sequence set to maximum
        return cnt;
    }

    
    /**
     * Returns the number of bp differences between two 2-bit encoded longs.
     * @param seq1 2-bit encoded sequence
     * @param seq2 2-bit encoded sequence
     * @return count of the divergence 
     */
    public static byte seqDifferences(long seq1, long seq2) {
        long mask = 3;
        byte cnt = 0;
        long diff = seq1 ^ seq2;
        for (int x = 0; x < chunkSize; x++) {
            if ((diff & mask) > 0) {
                cnt++;
            }
            diff = diff >> 2;
            // System.out.println("v = " + v);
        }
        return cnt;
    }
    /**
     * Returns the number of sequencing differences between two 2-bit encoded longs.
     * Maximum divergence is used to save time when only interested in very similar 
     * sequences.
     * @param seq1 2-bit encoded sequence
     * @param seq2 2-bit encoded sequence
     * @param lengthOfComp number of sites to compare
     * @param maxDivergence threshold for counting divergence upto
     * @return count of the divergence (above the maxDivergence, chunkSize is returned)
     */
    public static byte seqDifferencesForSubset(long seq1, long seq2, int lengthOfComp, int maxDivergence) {
        long mask = 3;
        byte cnt = 0;
        long diff = seq1 ^ seq2;
        diff = diff >> (2 * (chunkSize - lengthOfComp));  //shift to 5' end of sequence
        for (int x = 0; x < lengthOfComp && cnt < maxDivergence; x++) {
            if ((diff & mask) > 0) {
                cnt++;
            }
            diff = diff >> 2;
        }
        return cnt;
    }

    /**
     * Trim the poly-A off the sequence string
     * @param s input sequence
     * @return sequence with polyA removed
     */
    public static String removePolyAFromEnd(String s) {
        int index = s.length() - 1;
        while (s.charAt(index) == 'A') {
            index--;
            if (index < 1) {
                return null;
            }
        }
        return s.substring(0, index + 1);
    }
}
