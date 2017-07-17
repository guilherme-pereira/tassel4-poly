/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

/**
 *
 * @author terry
 */
public class TerryPipelines {

    public static void runDiscoverySNPCallerPlugin() {

        String baseDirTestFuzzyPoz282 = "/Users/terry/users/jeff/";
        String[] testFuzzyPoz282Args = new String[]{
            "-i", baseDirTestFuzzyPoz282 + "C08L7ACXX_6.tbt.byte",
            "-y", // use TagsByTaxaByte
            "-o", baseDirTestFuzzyPoz282 + "hapmap",
            "-m", baseDirTestFuzzyPoz282 + "MGP1_low_vol_min3_wPosit.topm.bin",
            //            "-mUpd", baseDir+"",
            //"-ref", "maize_agp_v2.fasta",
            //"-LocusBorder", "150",
            "-mnF", "0.8",
            "-mnMAF", "0.01",
            "-mnMAC", "99999", // this will never be satified: this way -mnMAF overrides it
            "-mnLCov", "0.1", // Minimum locus coverage (proportion of Taxa)
            //            "-inclGaps",  // Include sites where major or minor allele is a GAP
            //            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-s", "10", // Start chromosome
            "-e", "10" // End chromosome
        };

        DiscoverySNPCallerPlugin plugin = new DiscoverySNPCallerPlugin();
        plugin.setParameters(testFuzzyPoz282Args);
        plugin.performFunction(null);

    }

    public static void main(String[] args) {
        runDiscoverySNPCallerPlugin();
    }
}
