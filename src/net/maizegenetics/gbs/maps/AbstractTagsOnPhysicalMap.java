/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.maps;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;
import net.maizegenetics.gbs.tagdist.AbstractTags;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;

/**
 * Abstract TagsOnPhysicalMap object.  Abstract alignment is implemented by two types of TOPM classes
 * {@link TagsOnPhysicalMap} and {@link TagsOnPhysMapHDF5}.   TagsOnPhysicalMap is fully loaded into memory when
 * read, and it can only support one mapping position (the files are either binary or text).  TagsOnPhysMapHDF5
 * can support a large number of mapping position, and only a portion of the data is loaded into memory.
 *
 * In most cases, TagsOnPhysMapHDF5 should be preferred going forward.
 *
 * @author edbuckler
 */
public abstract class AbstractTagsOnPhysicalMap extends AbstractTags implements TOPMInterface {
    protected int[] bestChr; // 4 bytes
    // 4 bytes
    //if these disagree with the location, then set the p to negative
    // 1+4+1+4+4+1+8+8+1+1 = 33 bytes per position + 16 bytes for a two long tag + 1 byte for tagLength in bases = 50 bytes
    // ~50 bytes per position.  If we have 10 million tags then this will be at 500M byte data structure.
    protected int[] indicesOfSortByPosition;
    protected int myMaxVariants = 8;
    protected byte[] multimaps = null; // number of locations this tagSet maps to; unknown = Byte.MIN_VALUE; multiple, but unknown number = 99
    // number of locations this tagSet maps to; unknown = Byte.MIN_VALUE; multiple, but unknown number = 99
    protected int[] bestStartPos = null; // chromosomal position of the barcoded end of the tag  // 4 bytes
    // chromosomal position of the barcoded end of the tag  // 4 bytes
    protected byte[] bestStrand = null; // 1 = same sense as reference FASTA file.  -1 = opposite sense.  unknown = Byte.MIN_VALUE  // 1 byte
    // 1 = same sense as reference FASTA file.  -1 = opposite sense.  unknown = Byte.MIN_VALUE  // 1 byte
    protected int myNumTags = 0;
    //TODO
//    public int tagNum;  //remove this and set to above
//    public int maxVariants;
    
    protected byte[][] variantDefs; // allele state - A, C, G, T or some indel definition  // myMaxVariants bytes [tag][variant]
    // allele state - A, C, G, T or some indel definition  // myMaxVariants bytes [tag][variant]
    protected byte[][] variantOffsets; // offset from position minimum, maximum number of variants is defined above  // myMaxVariants bytes [tag][variant]
    // offset from position minimum, maximum number of variants is defined above  // myMaxVariants bytes [tag][variant]
    protected int[] myChromosomes = null; //sort ascending
    protected int[][] myUniquePositions = null; //dimensions [chromosome][positions] note ragged, and sort ascending

    public AbstractTagsOnPhysicalMap() {
    }

    
    @Override
    public int getMaxNumVariants() {
        return myMaxVariants;
    }
    
    @Override
    public byte getMultiMaps(int index) {
        return multimaps[index];
    }
    
    @Override
    public int getChromosome(int index) {
        return bestChr[index];
    }

    @Override
    public int getSize() {
        return myNumTags;
    }

    @Override
    public int getStartPosition(int index) {
        return bestStartPos[index];
    }

    @Override
    public byte getStrand(int tagIndex) {
        return bestStrand[tagIndex];
    }

    @Override
    public byte[][] getVariantDef() {
        byte[][] result = new byte[getTagCount()][myMaxVariants];
        for (int i = 0; i < getTagCount(); i++) {
            for (int j = 0; j < myMaxVariants; j++) {
                result[i][j] = getVariantDef(i, j);
            }
        }
        return result;
    }
    
    /**
     * Provides the reference
     * @return 
     */
    protected byte[][] getVariantDefByReference() {
        return variantDefs;
    }
    
    protected byte[][] getVariantOffByReference() {
        return variantOffsets;
    }
    

    @Override
    public byte getVariantDef(int tagIndex, int variantIndex) {
        if((variantDefs[tagIndex]==null)||(variantDefs[tagIndex].length<=variantIndex)) return TOPMInterface.BYTE_MISSING;
        return variantDefs[tagIndex][variantIndex];
    }

    /**
     * Returns an array containing all variant definitions for the tag at the
     * supplied index.
     */
    @Override
    public byte[] getVariantDefArray(int tagIndex) {
        return variantDefs[tagIndex];
//        if(variantDefs[tagIndex]==null) return null;
//        byte[] result = new byte[variantDefs[tagIndex].length];
//        for (int i = 0; i < variantDefs[tagIndex].length; i++) {
//            result[i] = getVariantDef(tagIndex, i);
//        }
//        return result;
    }

    @Override
    public byte[][] getVariantOff() {
        byte[][] result = new byte[getTagCount()][myMaxVariants];
        for (int i = 0; i < getTagCount(); i++) {
            System.arraycopy(variantOffsets[i], 0, result[i], 0, myMaxVariants);
        }
        return result;
    }

    @Override
    public byte getVariantPosOff(int tagIndex, int variantIndex) {
        if((variantOffsets[tagIndex]==null)||(variantOffsets[tagIndex].length<=variantIndex)) return TOPMInterface.BYTE_MISSING;
        return variantOffsets[tagIndex][variantIndex];
    }

    /**
     * Returns an array containing all variant position offsets for the tag at
     * the supplied index.
     */
    @Override
    public byte[] getVariantPosOffArray(int tagIndex) {
        return variantOffsets[tagIndex];
    }
    
    public String printRow(int row) {
        StringBuilder sb = new StringBuilder();
        sb.append(sb);
        //long
        sb.append(BaseEncoder.getSequenceFromLong(this.getTag(row)) + "\t");
        sb.append(printWithMissing(tagLength[row])); sb.append("\t");
        sb.append(printWithMissing(multimaps[row]) + "\t");
        sb.append(printWithMissing(bestChr[row]) + "\t");
        sb.append(printWithMissing(bestStrand[row]) + "\t");
        sb.append(printWithMissing(bestStartPos[row]) + "\t");
        sb.append(printWithMissing(getEndPosition(row)) + "\t");
        sb.append(printWithMissing(getDivergence(row)) + "\t");
        for (int j = 0; j < myMaxVariants; j++) {
            sb.append(printWithMissing(getVariantPosOff(row, j)) + "\t");
            byte vd=getVariantDef(row, j);
            if(vd==TOPMInterface.BYTE_MISSING) {sb.append(printWithMissing(vd) + "\t");}
            else {
                byte genotype=AlignmentUtils.getDiploidValue(vd, vd);
                sb.append(NucleotideAlignmentConstants.getNucleotideIUPAC(genotype) + "\t");
            }
        }
        sb.append(printWithMissing(getDcoP(row)) + "\t");
        sb.append(printWithMissing(getMapP(row)) + "\t");
 //       System.out.println("Line:"+row+":"+sb.toString());
        return sb.toString();
    }
    
    public static String printWithMissing(byte b) {
        if (b == Byte.MIN_VALUE) {
            return "*";
        }
        return Byte.toString(b);
    }

    public static String printWithMissing(int i) {
        if (i == Integer.MIN_VALUE) {
            return "*";
        }
        return Integer.toString(i);
    }

    @Override
    public void writeTextFile(File outfile) {
        try {
            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfile), 65536));
            fw.writeBytes(myNumTags + "\t" + tagLengthInLong + "\t" + myMaxVariants + "\n");
            for (int row = 0; row < myNumTags; row++) {
                fw.writeBytes(printRow(row) + "\n");
            }
            fw.flush();
            fw.close();
        } catch (Exception e) {
            System.out.println("Catch in writeTextFile file e=" + e);
            e.printStackTrace();
        }
        System.out.println("Number of tags in file:" + myNumTags);
    }

    @Override
    public int[] getChromosomes() {
        if(myChromosomes==null) populateChrAndVarPositions();
        return myChromosomes;
    }
        
    @Override
    public int[] getUniquePositions(int chromosome) {
        if(myUniquePositions==null) populateChrAndVarPositions();
        return myUniquePositions[chromosome];
    }

    /**
     * Creates the arrays for the all the positions for each chromosome as defined
     * by the variants.  Call this method after remapping or after loading the file, and
     * it will be used for SNP calling.
     */
    protected void populateChrAndVarPositions() {
        long chrSum = 0;
        System.out.println("chrSum" + chrSum);
        TreeMap<Integer, TreeSet<Integer>> theChrs = new TreeMap<Integer, TreeSet<Integer>>();
        for (int i = 0; i < myNumTags; i++) {
            int chr = getChromosome(i);
            if (chr != TOPMInterface.INT_MISSING) {
                if (!theChrs.containsKey(chr)) {
                    theChrs.put(chr, new TreeSet<Integer>());
                }
                TreeSet<Integer> thePos = theChrs.get(chr);
                int startPos = getStartPosition(i);
                byte[] varOffs = getVariantPosOffArray(i);
                if(varOffs==null) continue;
                for (byte b : varOffs) {
                    thePos.add((int) (startPos + b));
                }
            }
        }
        myChromosomes = new int[theChrs.size()];
        myUniquePositions = new int[theChrs.size()][];
        int cnt = 0;
        for (Entry<Integer, TreeSet<Integer>> aChr : theChrs.entrySet()) {
            myUniquePositions[cnt] = new int[aChr.getValue().size()];
            int p = 0;
            for (int ls : aChr.getValue()) {
                myUniquePositions[cnt][p++] = ls;
            }
            myChromosomes[cnt++] = aChr.getKey();
//            System.out.printf("Chr:%d TagStart:%d %n", myChromosomes[cnt - 1], myUniquePositions[cnt - 1].length);
        }
    }
    
    @Override
    public Locus[] getLoci() {
        int[] chrs = getChromosomes();
        Locus[] result = new Locus[chrs.length];

        for (int i = 0; i < result.length; i++) {
            result[i] = new Locus(chrs[i] + "", chrs[i] + "", -1, -1, null, null);
        }
        return result;
    }

    @Override
    public Locus getLocus(int tagIndex) {
        if (bestChr[tagIndex] == TOPMInterface.INT_MISSING) {
            return null;
        } //Return null for unmapped tags
        return new Locus(bestChr[tagIndex] + "", bestChr[tagIndex] + "", -1, -1, null, null);
    }
    
   void initPhysicalSort() {
        System.out.println("initPhysicalSort");
        indicesOfSortByPosition = new int[myNumTags];
        for (int i = 0; i < indicesOfSortByPosition.length; i++) {
            indicesOfSortByPosition[i] = i;
        }
        Swapper swapperPos = new Swapper() {
            public void swap(int a, int b) {
                int t1;
                t1 = indicesOfSortByPosition[a];
                indicesOfSortByPosition[a] = indicesOfSortByPosition[b];
                indicesOfSortByPosition[b] = t1;
            }
        };
        IntComparator compPos = new IntComparator() {
            public int compare(int a, int b) {
                int index1 = indicesOfSortByPosition[a];
                int index2 = indicesOfSortByPosition[b];
                if (bestChr[index1] < bestChr[index2]) {
                    return -1;
                }
                if (bestChr[index1] > bestChr[index2]) {
                    return 1;
                }
                if (bestStartPos[index1] < bestStartPos[index2]) {
                    return -1;
                }
                if (bestStartPos[index1] > bestStartPos[index2]) {
                    return 1;
                }
                if (bestStrand[index1] < bestStrand[index2]) {
                    return -1;
                }
                if (bestStrand[index1] > bestStrand[index2]) {
                    return 1;
                }
                for (int i = 0; i < tagLengthInLong; i++) {
                    if (tags[i][index1] < tags[i][index2]) {
                        return -1;
                    }
                    if (tags[i][index1] > tags[i][index2]) {
                        return 1;
                    }
                }
                return 0;
            }
        };
        System.out.println("Position index sort begin.");
        GenericSorting.quickSort(0, indicesOfSortByPosition.length, compPos, swapperPos);
        System.out.println("Position index sort end.");
    }

    public void writeBinaryFile(File outFile) {
        writeBinaryFile(outFile, Integer.MAX_VALUE, false, false, Float.NaN, true);
    }

    protected void writeBinaryFile(File outFile, boolean binary) {
        writeBinaryFile(outFile, Integer.MAX_VALUE, false, false, Float.NaN, binary);
    }

    /**
     * TODO need to add LocusList to test, whether to include in output, move this to abstract
     * @param outFile
     * @param minResolution
     * @param requirePhysPosition
     * @param requireDCOMap
     * @param minDCOP
     * @param binary
     */
    public void writeBinaryFile(File outFile, int minResolution, boolean requirePhysPosition, boolean requireDCOMap, float minDCOP, boolean binary) {
        int hapsOutput = 0;
        try {
            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile), 4000000));
            if (requirePhysPosition) {
                fw.writeInt(mappedTags()[0]);
            } // the index 0 provides the number of tags with unique positions
            else {
                fw.writeInt(myNumTags);
            }
            fw.writeInt(tagLengthInLong);
            fw.writeInt(myMaxVariants);
            for (int row = 0; row < myNumTags; row++) {
                if ((requirePhysPosition == true) && (bestChr[row] == Integer.MIN_VALUE)) {
                    continue;
                }
                for (int j = 0; j < tagLengthInLong; j++) {
                    fw.writeLong(tags[j][row]);
                }
                fw.writeByte(tagLength[row]);
                fw.writeByte(multimaps[row]);
                fw.writeInt(bestChr[row]);
                fw.writeByte(bestStrand[row]);
                fw.writeInt(bestStartPos[row]);
                fw.writeInt(getEndPosition(row));
                fw.writeByte(getDivergence(row));
                for (int j = 0; j < myMaxVariants; j++) {
                    fw.writeByte(getVariantPosOff(row,j));
                    fw.writeByte(getVariantDef(row,j));
                }
                fw.writeByte(getDcoP(row));
                fw.writeByte(getMapP(row));
                hapsOutput++;
            }
            fw.flush();
            fw.close();
            System.out.println("Tag positions written to:" + outFile.toString());
            System.out.println("Number of tags in file:" + hapsOutput);
        } catch (Exception e) {
            System.err.println("Catch in writing output file e=" + e);
        }
    }

    /**
     * @return An int[] result where : result[0] = The number of tags with a
     * unique physical positions in this file (i.e. , tags for which the
     * bestChr number is known). result[1] = The number of tags which align
     * to multiple positions (i.e., where multimaps[tagIndex] > 0)
     *
     */
    public int[] mappedTags() {
        int[] result = {0, 0};
        int unique = 0;
        int multi = 1; // the indices of result
        for (int row = 0; row < myNumTags; row++) {
            if (bestChr[row] == Integer.MIN_VALUE) {
                if (multimaps[row] > 0) {
                    result[multi]++;
                }
            } else {
                result[unique]++;
            }
        }
        return result;
    }

    public void printRows(int numRows) {
        for (int i = 0; i < numRows; i++) {
            System.out.println(printRow(i));
        }
    }
    
    public String printRow(int row, boolean byPosition) {
        if (byPosition) {
            return printRow(indicesOfSortByPosition[row]);
        }
        return printRow(row);
    }

    public void printRows(int numRows, boolean requirePhysPosition, boolean byPosition) {
        int outCount = 0;
        for (int i = 0; outCount < numRows; i++) {
            int r = (byPosition) ? indicesOfSortByPosition[i] : i;
            if ((requirePhysPosition == true) && (bestChr[r] < 1)) {
                continue;
            }
            System.out.println(printRow(r));
            outCount++;
        }
    }

    public void printRows(int numRows, boolean requirePhysPosition, int printChr) {
        int outCount = 0;
        boolean byPosition = true;
        for (int i = 0; outCount < numRows; i++) {
            int r = (byPosition) ? indicesOfSortByPosition[i] : i;
            if ((requirePhysPosition == true) && (bestChr[r] != printChr)) {
                continue;
            }
            System.out.println(printRow(r));
            outCount++;
        }
    }
    
    protected long[][] getTagsArray() {
        return tags;
    }
    
    protected byte[] getTagLengthArray() {
        return tagLength;
    }
    
}
