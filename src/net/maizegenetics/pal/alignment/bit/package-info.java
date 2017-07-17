/**
 * Background on Bit Encoding in TASSEL:<p></p>
 TASSEL has been optimized for several things.  One of TASSEL&rsquo;s optimization is for very fast comparison
 of the sequence of two taxa for millions of sites (e.g. genetic distance or imputation), or the comparison of
 two sites for thousands of taxa (e.g. LD).
 <p></p>
 To do this, the presence of an allele is coded with a single bit  (presence= 1,  absence=0).
 With this encoding, then bit level operations can be used calculate genetic distance or LD.
 <p></p>
 So how is this done?
 The main assumptions:  (1) we generally only care about the two most alleles at any given locus.
 This may not be true for GWAS at a specific site, but across a genome it is a close enough approximation.
 (2) phase is generally not tracked in the bit array, however, if everything is phased, then the haplotypes can be used.
 (3) If allele 2 is unknown, the genotype is assumed to be homozygous for allele 1.
 e.g. A/N assumed to be A/A.  Unknown heterozygotes are not support.  Completely missing genotypes are supported.
 <p></p>
 So the approach is slightly lossy &ndash; a bit of information is lots &ndash; just like some compression algorithms,
 in exchange for some very nice features.
 <p></p>
 <img src="bit\BitEncodingTable.png"/>
 <p></p>
  While the bit level conserves memory space, it is not the main purpose.  The above sequence if store as text would
 take 16 bytes (8 sites x 2 alleles in a diploid).  Using TASSEL half byte approach, it fits into 8 bytes.  The
 bit version would take one byte for the major allele and one for the minor allele for a total of 2 bytes.
 <p></p>
 In practice, we store these bits of allele presence in sets of 64 (a long).  This works very efficiently as most
 computer architectures are also 64-bit, so when we process the sequence everything is being done at 64 sites at a time.
 THESE RESULTS IN THE ALGORITHMS THAT USE BIT LEVEL ENCODING TO BE 25 TO 50 TIMES FASTER.

 */
package net.maizegenetics.pal.alignment.bit;