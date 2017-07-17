/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.homology;

import net.maizegenetics.gbs.util.BaseEncoder;

/**
 *
 * @author Fei Lu
 */
public class PEParseBarcodeRead extends ParseBarcodeRead {
    /**Likely read end of backward read, including the second cut site and the barcoded adaptor*/
    String[] bLikelyReadEnd;
    
    /**
     * Construct empty object from key file, enzyme and Illumina data from a lane. Set up likely end and barcode for parsing data
     * @param keyFile
     *        the key file
     * @param enzyme
     *        enzyme used
     * @param flowcell
     *        flowcell name
     * @param lane 
     *        lane number
     */
    public PEParseBarcodeRead(String keyFile, String enzyme, String flowcell, String lane) {
        super(keyFile, enzyme, flowcell, lane);
        bLikelyReadEnd = new String[likelyReadEnd.length/2];
        for (int i = 0; i < bLikelyReadEnd.length; i++) {
            bLikelyReadEnd[i] = likelyReadEnd[i];
        }
    }
    
    /**
     * Parse sequence pair to a PE tag and taxon.
     * @param seqSF
     *        Forward sequence
     * @param qualSF
     *        Quality string of forward sequence
     * @param seqSB
     *        Backward sequence
     * @param qualSB
     *        Quality string of backward sequence
     * @param fastq
     *        Boolean value of if the data is fastq format
     * @param minQual
     *        Minimum of quality value of sequence, default is 0 (quality is not controled)
     * @param tagLengthInLong
     *        Tag length in Long primitive data type, defaut is 8 (8 * 32 = 256 bp)
     * @return Parsed PE read (forward tag, backward tag and taxon)
     */
    public PEReadBarcodeResult parseReadIntoTagAndTaxa(String seqSF, String qualSF, String seqSB, String qualSB, boolean fastq, int minQual, int tagLengthInLong) {
        boolean ifMatch = false;
        for (int i = 0; i < initialCutSiteRemnant.length; i++) {
            if (seqSB.startsWith(initialCutSiteRemnant[i])) {
                ifMatch = true;
                break;
            }
        }
        if (!ifMatch) return null;
        if ((minQual > 0) && (qualSF != null)) {
            int firstBadBase = BaseEncoder.getFirstLowQualityPos(qualSF, minQual);
            if (firstBadBase < (maxBarcodeLength + tagLengthInLong * BaseEncoder.chunkSize)) {
                return null;
            }
        }
        if ((minQual > 0) && (qualSB != null)) {
            int firstBadBase = BaseEncoder.getFirstLowQualityPos(qualSB, minQual);
            if (firstBadBase < (maxBarcodeLength + tagLengthInLong * BaseEncoder.chunkSize)) {
                return null;
            }
        }
        int miss1 = -1, miss2 = -1;
        if (fastq) {
            miss1 = seqSF.indexOf('N');
            miss2 = seqSB.indexOf('N');
        } else {
            miss1 = seqSF.indexOf('.');
            miss2 = seqSB.indexOf('.');
        }
        if ((miss1 != -1) && (miss1 < (maxBarcodeLength + tagLengthInLong * BaseEncoder.chunkSize))) {
            return null;  //bad sequence so skip
        }
        if ((miss2 != -1) && (miss2 < (tagLengthInLong * BaseEncoder.chunkSize))) {
            return null;  //bad sequence so skip
        }
        Barcode bestBarcode = findBestBarcode(seqSF, maximumMismatchInBarcodeAndOverhang);
        if (bestBarcode == null) {
            return null;  //overhang missing so skip
        }
        String genomicSeqF = seqSF.substring(bestBarcode.barLength, seqSF.length());
        ShortReadBarcodeResult tagProcessingResultsF = removeSeqAfterSecondCutSite(genomicSeqF, (short)(tagLengthInLong * BaseEncoder.chunkSize), tagLengthInLong, bestBarcode.getTaxaName(), true);
        String genomicSeqB = this.removeSeqAfterBarcode(bestBarcode.barcodeS, genomicSeqF, seqSB);
        ShortReadBarcodeResult tagProcessingResultsB = removeSeqAfterSecondCutSite(genomicSeqB, (short)(tagLengthInLong * BaseEncoder.chunkSize), tagLengthInLong, bestBarcode.getTaxaName(), false);
        return new PEReadBarcodeResult(tagProcessingResultsF, tagProcessingResultsB);
    }
    
    /**
     * Remove the sequence after barcode adaptor for the backward sequence
     * @param bestBarcodeS
     * @param genomicSeqF
     * @param seqSB
     * @return 
     */
    private String removeSeqAfterBarcode (String bestBarcodeS, String genomicSeqF, String seqSB) {
        String barcodeGenomicS = bestBarcodeS+genomicSeqF.substring(0, readEndCutSiteRemnantLength);
        String rBarcodeGenomicS = BaseEncoder.getReverseComplement(barcodeGenomicS);
        int index = seqSB.indexOf(rBarcodeGenomicS);
        if (index == -1) return seqSB;
        return seqSB.substring(0, index+readEndCutSiteRemnantLength);
    }
    
    /**
     * Remove the sequence of the second cutsite
     * @param seq
     *        sequence need to be processed
     * @param maxLength
     *        maximum length of padded sequence, e.g. 32 * 8 = 256 bp. 256 is the maxLength
     * @param tagLengthInLong
     *        PE Tag in Long primitive data type, e.g. 8
     * @param taxaName
     *        taxon name
     * @param ifForward
     *        Forward or backward sequence
     * @return 
     */
    public ShortReadBarcodeResult removeSeqAfterSecondCutSite(String seq, short maxLength, int tagLengthInLong, String taxaName, boolean ifForward) {
        //this looks for a second restriction site or the common adapter start, and then turns the remaining sequence to AAAA
        int cutSitePosition = 9999;
        ShortReadBarcodeResult returnValue = new ShortReadBarcodeResult(seq);
        returnValue.taxonName = taxaName;
        //Look for cut sites, starting at a point past the length of the initial cut site remnant that all reads begin with
        String match = null;
        for (String potentialCutSite : likelyReadEnd) {
            int p = seq.indexOf(potentialCutSite, 1);
            if ((p > 1) && (p < cutSitePosition)) {
                cutSitePosition = p;
                match = potentialCutSite;
            }
        }
        if (theEnzyme.equalsIgnoreCase("ApeKI") && cutSitePosition == 2
                && (match.equalsIgnoreCase("GCAGC") || match.equalsIgnoreCase("GCTGC"))) {  // overlapping ApeKI cut site: GCWGCWGC
            seq = seq.substring(3, seq.length());  // trim off the initial GCW from GCWGCWGC
            cutSitePosition = 9999;
            returnValue.unprocessedSequence = seq;
            if (ifForward) {
                for (String potentialCutSite : likelyReadEnd) {
                    int p = seq.indexOf(potentialCutSite, 1);
                    if ((p > 1) && (p < cutSitePosition)) {
                        cutSitePosition = p;
                    }
                }
            }
            else {
                for (String potentialCutSite : bLikelyReadEnd) {
                    int p = seq.indexOf(potentialCutSite, 1);
                    if ((p > 1) && (p < cutSitePosition)) {
                        cutSitePosition = p;
                    }
                }
            }
            
        }

        if (cutSitePosition < maxLength) {  // Cut site found
            //Trim tag to sequence up to & including the cut site
            returnValue.length = (short) (cutSitePosition + readEndCutSiteRemnantLength);
            returnValue.processedSequence = seq.substring(0, cutSitePosition + readEndCutSiteRemnantLength);
        } else {
            if (seq.length() <= 0) {
                //If cut site is missing because there is no sequence
                returnValue.processedSequence = "";
                returnValue.length = 0;
            } else {
                //If cut site is missing because it is beyond the end of the sequence (or not present at all)
                returnValue.length = (short) Math.min(seq.length(), maxLength);
                returnValue.processedSequence = (seq.substring(0, returnValue.length));
            }
        }

        //Pad sequences shorter than max. length with A
        if (returnValue.length < maxLength) {
            returnValue.paddedSequence = returnValue.processedSequence + nullS;
            returnValue.paddedSequence = returnValue.paddedSequence.substring(0, maxLength);
        } else {
            //Truncate sequences longer than max. length
            returnValue.paddedSequence = returnValue.processedSequence.substring(0, maxLength);
            returnValue.length = maxLength;
        }
        returnValue.read = BaseEncoder.getLongArrayFromSeq(returnValue.paddedSequence);
        return returnValue;
    }
}
