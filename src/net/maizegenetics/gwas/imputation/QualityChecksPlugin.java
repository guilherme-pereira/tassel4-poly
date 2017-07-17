package net.maizegenetics.gwas.imputation;

import java.awt.Frame;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import javax.swing.ImageIcon;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.DefaultXYDataset;

import net.maizegenetics.baseplugins.ExportPlugin;
import net.maizegenetics.baseplugins.GenotypeSummaryPlugin;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.Utils;

public class QualityChecksPlugin extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(QualityChecksPlugin.class);
	private String pedigreeFile;
	private int windowSizeForR2 = 15;
	private double minNonMissingProportionForTaxon = 0.1;
	private double minNonMissingProportionForSNP = 0.1;
	private String avgr2Filename = null;
	private String avgr2Plotname = null;
	private String propNonconsensusFilename = null;
	private String summaryFilename = null;
	private String logfileName = null;
	private boolean hasFileAppender = false;
	
//	public enum checkType {AVERAGE_R2, NONCONSENSUS_PROPORTION, NONCONSENSUS_SITES};
//	private ArrayList<checkType> analysisList = new ArrayList<checkType>();
	
	public QualityChecksPlugin(Frame parentFrame) {
		super(parentFrame, false);
		BasicConfigurator.configure();
	}
	
	@Override
	public DataSet performFunction(DataSet input) {

		if (!hasFileAppender) {
			try{
				myLogger.addAppender(new FileAppender(new PatternLayout("%-5p [%t]: %m%n"), logfileName));
			} catch(Exception e) {
				e.printStackTrace();
			}
			hasFileAppender = true;
		}
		
		List<Datum> datumList = input.getDataOfType(Alignment.class);
		ArrayList<PopulationData> familyList;
		if (pedigreeFile != null) {
			familyList = PopulationData.readPedigreeFile(pedigreeFile);
		} else {
			familyList = new ArrayList<PopulationData>();
		}
		
		for (Datum datum : datumList) {
			Alignment anAlignment = (Alignment) datum.getData();
			if (pedigreeFile == null) {
				processFamily(anAlignment, null);
			} else {
				for (PopulationData family : familyList) {
					String[] names = new String[family.members.size()];
					family.members.toArray(names);
					Alignment align = FilterAlignment.getInstance(anAlignment, new SimpleIdGroup(names), false);
					processFamily(align, family.name);
				}
			}
			
		}
		return null;
	}
	
	private void processFamily(Alignment align, String familyname) {
		myLogger.info("\nResults for chromosome " + align.getLocusName(0) + ", family " + familyname);
		align = preFilterAlignment(align);

		if (avgr2Filename != null || avgr2Plotname != null) {
			calculateAverageR2ForSnps(align, familyname);
		}
		
		if (propNonconsensusFilename != null) {
			double[] proportion = calculateProportionNonConsensusPerTaxon(align);
			saveProportionNonConsensusToFile(proportion, align, addFamilyToFilename(propNonconsensusFilename, familyname, align.getLocusName(0), ".txt"));
		}
		
		if (summaryFilename != null) {
			runAndExportGenotypeSummaryForTaxa(align, addFamilyToFilename(summaryFilename, familyname, align.getLocusName(0), ".txt"));
		}
	}
	
	private String addFamilyToFilename(String filename, String family, String chr, String extension) {
		if (!extension.startsWith(".")) extension = "." + extension;
		
		StringBuilder sb = new StringBuilder();
		if (filename.endsWith(extension)) {
			sb.append(filename.substring(0, filename.length() - extension.length()));
		} else {
			sb.append(filename);
		}
		
		if (family == null) {
			if (chr != null) sb.append(".chr").append(chr);
			sb.append(extension);
		} else {
			family = family.replace('/', '_');
			sb.append(".").append(family);
			if (chr != null) sb.append(".chr").append(chr);
			sb.append(extension);
		}
		
		return sb.toString();
	}
	
	public Alignment preFilterAlignment(Alignment align) {
		int ntaxa = align.getSequenceCount();
		int nsites = align.getSiteCount();
		int nTaxaGametes = 2 * nsites;
		int nSiteGametes = 2 * ntaxa;
		int minTaxaGametes = (int) Math.ceil(nTaxaGametes * minNonMissingProportionForTaxon);
		int minSiteGametes = (int) Math.ceil(nSiteGametes * minNonMissingProportionForSNP);
		
		myLogger.info("preFilter Alignment: ");
		myLogger.info("Before filter alignment has " + ntaxa + " taxa and " + nsites + " sites.");
		
		//create list of taxa with too much missing data
		LinkedList<Identifier> taxaDiscardList = new LinkedList<Identifier>();
		for (int t = 0; t < ntaxa; t++) {
			if (align.getTotalGametesNotMissingForTaxon(t) < minTaxaGametes) taxaDiscardList.add(align.getIdGroup().getIdentifier(t));
		}
		if (taxaDiscardList.size() > 0) {
			myLogger.info("\nThe following taxa will not be included in the analysis because the proportion of nonMissing data is below " + minNonMissingProportionForTaxon + ":\n");
			for (Identifier id:taxaDiscardList) myLogger.info(id.getFullName());
			
			Identifier[] ids = new Identifier[taxaDiscardList.size()];
			taxaDiscardList.toArray(ids);
			
			align = FilterAlignment.getInstanceRemoveIDs(align, new SimpleIdGroup(ids));
		}
		
		myLogger.info("After filtering for taxa, there are " + align.getSequenceCount() + " taxa.");
		
		//number of non-missing values per site
		int[] sitesToKeep = new int[nsites];
		int nsitesKept = 0;
		for (int s = 0; s < nsites; s++) {
			if (align.getTotalGametesNotMissing(s) >= minSiteGametes) sitesToKeep[nsitesKept++] = s;
		}
		
		if (nsitesKept < nsites) {
			myLogger.info(nsitesKept + " sites had more than " + minSiteGametes + " and were retained.");
			sitesToKeep = Arrays.copyOf(sitesToKeep, nsitesKept);
			align = FilterAlignment.getInstance(align, sitesToKeep);
		}
		
		return align;
	}
	
	private void calculateAverageR2ForSnps(Alignment align, String familyname) {
		
		//first filter out monomorphic sites
		int nsites = align.getSiteCount();
		int[] polysites = new int[nsites];
		int sitecount = 0;
		for (int s = 0; s < nsites; s++) {
			if (align.getMinorAlleleFrequency(s) > 0.15) polysites[sitecount++] = s;  
		}
		polysites = Arrays.copyOf(polysites, sitecount);
		align = FilterAlignment.getInstance(align, polysites);
		align = BitAlignment.getInstance(align, true);
		
		myLogger.info("Chromosome " + align.getLocusName(0) + ", family " + familyname + " has " + sitecount + " polymorphic snps.");
		
		nsites = align.getSiteCount();
		double[] avgRsq = new double[nsites];
		
		for (int s = 0; s < nsites; s++) {
			int start = Math.max(s - windowSizeForR2, 0);
			int end = Math.min(nsites - 1, s + windowSizeForR2);
			double sum = 0;
			double count = 0;
            BitSet sMj = align.getAllelePresenceForAllTaxa(s, 0);
            BitSet sMn = align.getAllelePresenceForAllTaxa(s, 1);

			for (int i = start; i <= end; i++) {
				if (i != s) {
					int[][] contig = new int[2][2];
		            BitSet iMj = align.getAllelePresenceForAllTaxa(i, 0);
		            BitSet iMn = align.getAllelePresenceForAllTaxa(i, 1);
		            contig[0][0] = (int) OpenBitSet.intersectionCount(sMj, iMj);
		            contig[1][0] = (int) OpenBitSet.intersectionCount(sMn, iMj);
		            contig[0][1] = (int) OpenBitSet.intersectionCount(sMj, iMn);
		            contig[1][1] = (int) OpenBitSet.intersectionCount(sMn, iMn);
		            double rsq = calculateRSqr(contig[0][0], contig[1][0], contig[0][1], contig[1][1], 4);
		            if (!Double.isNaN(rsq)) {
		            	sum += Math.sqrt(rsq);
		            	count++;
		            }
				}
			}
			
			if (count > 0) {
				avgRsq[s] = sum / count;
			} else {
				avgRsq[s] = Double.NaN;
			}
			
		}
		
		String chrname = align.getLocusName(0);
		if (avgr2Filename != null) saveToFileAverageR2(avgRsq, align, addFamilyToFilename(avgr2Filename, familyname, chrname, ".txt"));
		if (avgr2Plotname != null) plotAverageR2(avgRsq, align, addFamilyToFilename(avgr2Plotname, familyname, chrname, ".png"));

	}
	
    static double calculateRSqr(int countAB, int countAb, int countaB, int countab, int minTaxaForEstimate) {
        //this is the Hill & Robertson measure as used in Awadella Science 1999 286:2524
        double freqA, freqB, rsqr, nonmissingSampleSize;
        nonmissingSampleSize = countAB + countAb + countaB + countab;
        if (nonmissingSampleSize < minTaxaForEstimate) {
            return Double.NaN;
        }
        freqA = (double) (countAB + countAb) / nonmissingSampleSize;
        freqB = (double) (countAB + countaB) / nonmissingSampleSize;

        //Through missing data & incomplete datasets some alleles can be fixed this returns missing value
        if ((freqA == 0) || (freqB == 0) || (freqA == 1) || (freqB == 1)) {
            return Double.NaN;
        }

        rsqr = ((double) countAB / nonmissingSampleSize) * ((double) countab / nonmissingSampleSize);
        rsqr -= ((double) countaB / nonmissingSampleSize) * ((double) countAb / nonmissingSampleSize);
        rsqr *= rsqr;
        rsqr /= freqA * (1 - freqA) * freqB * (1 - freqB);
        return rsqr;
    }

    private void saveToFileAverageR2(double[] avgr2, Alignment align, String saveFilename) {
    		BufferedWriter bw = Utils.getBufferedWriter(saveFilename);
    		int nsites = align.getSiteCount();
    		try {
    			bw.write("Site\tchr\tpos\tr2");
    			bw.newLine();
    			
    			for (int s = 0; s < nsites; s++) {
    				bw.write(align.getSNPID(s));
    				bw.write("\t");
    				bw.write(align.getLocusName(s));
    				bw.write("\t");
    				bw.write(Integer.toString(align.getPositionInLocus(s)));
    				bw.write("\t");
    				bw.write(Double.toString(avgr2[s]));
    				bw.newLine();
    			}
    			bw.close();
    		} catch(IOException e) {
    			myLogger.error("error opening file for avgr2data\n" + e.getMessage() + e.getStackTrace());
    		}
    }
    
    private void plotAverageR2(double[] avgr2, Alignment align, String saveFilename) {
    		int nsites = align.getSiteCount();
    		String title = "Average R2 in " + windowSizeForR2 + " bp window, chromosome " + align.getLocusName(0);
    		String xLabel = "position(Mbp)";
    		String yLabel ="Average R-squared";
    		DefaultXYDataset xydata = new DefaultXYDataset();
    		double[][] dataset = new double[2][nsites];
    		for (int s = 0; s < nsites; s++) {
    			dataset[0][s] = ((double) align.getPositionInLocus(s)) / 1000000.0 ;
    		}
    		dataset[1] = avgr2;
    		xydata.addSeries("avgr2", dataset);
    		JFreeChart chart = ChartFactory.createScatterPlot(title, xLabel, yLabel, xydata, PlotOrientation.VERTICAL, false, false, false);
    		try {
    			ChartUtilities.saveChartAsPNG(new File(saveFilename), chart, 800, 300);
    		} catch (IOException e) {
    			myLogger.error("error saving png in plotAverageR2\n" + e.getMessage() + e.getStackTrace());
    		}
    }
    
	private double[] calculateProportionNonConsensusPerTaxon(Alignment align) {
		double maxMaf = 0.05;
		int ntaxa = align.getSequenceCount();
		int nsites = align.getSiteCount();
		
		OpenBitSet lowmaf = new OpenBitSet(nsites);
		for (int s = 0; s < nsites; s++) {
			if (align.getMinorAlleleFrequency(s) < maxMaf) lowmaf.set(s);
		}
		
//		myLogger.info("taxa = " + align.getSequenceCount() + ", sites = " + align.getSiteCount());
//		myLogger.info("lowmaf count = " + lowmaf.cardinality());
		
		try {
			align.optimizeForTaxa(null);
		} catch(Exception e) {
			align = BitAlignment.getInstance(align, false);
		}
		
		
		double[] proportionMinor = new double[ntaxa];
//		double nmono = lowmaf.cardinality();
		for (int t = 0; t < ntaxa; t++) {
			BitSet major = align.getAllelePresenceForAllSites(t, 0);
			BitSet minor = align.getAllelePresenceForAllSites(t, 1);
			long minorCount = OpenBitSet.intersectionCount(lowmaf, minor);
			 
			OpenBitSet notMissing = new OpenBitSet(major.getBits(), major.getNumWords());
			notMissing.union(minor);
			long notMissingCount = OpenBitSet.unionCount(lowmaf, notMissing);
			proportionMinor[t] = ((double) minorCount) / ((double) notMissingCount);
		}
		
		return proportionMinor;
	}
	
	private void saveProportionNonConsensusToFile(double[] propNonconsensus, Alignment align, String saveFilename) {
    		BufferedWriter bw = Utils.getBufferedWriter(saveFilename);
    		int ntaxa = align.getSequenceCount();
    		try {
    			bw.write("Taxon\tchr\tpropNC");
    			bw.newLine();
    			String chr = align.getLocusName(0);
    			for (int t = 0; t < ntaxa; t++) {
    				bw.write(align.getFullTaxaName(t));
    				bw.write("\t");
    				bw.write(chr);
    				bw.write("\t");
    				bw.write(Double.toString(propNonconsensus[t]));
    				bw.newLine();
    			}
    			bw.close();
    		} catch(IOException e) {
    			myLogger.error("error opening file for proportion nonconsensus\n" + e.getMessage() + e.getStackTrace());
    		}
	}
	
	private void runAndExportGenotypeSummaryForTaxa(Alignment align, String outfile) {
		
		GenotypeSummaryPlugin gsp = new GenotypeSummaryPlugin(null,false);
		gsp.setCaculateOverview(false);
		gsp.setCalculateSiteSummary(false);
		gsp.setCalculateTaxaSummary(true);
		DataSet result = gsp.performFunction(new DataSet(new Datum("alignment", align, "alignment"), this));
		ExportPlugin exporter = new ExportPlugin(null, false);
		exporter.setSaveFile(outfile);
		exporter.performFunction(result);
	}
	
//	private void saveNonConsensusSites(Alignment align) {
//		
//	}
	
	@Override
	public void setParameters(String[] args) {
		if (args == null || args.length == 0) {
			myLogger.error(getUsage());
			return;
		}
		
		int narg = args.length;
		for (int i = 0; i < narg - 1; i++) {
			if (args[i].equals("-p") || args[i].equalsIgnoreCase("-pedigrees")) {
				pedigreeFile = args[++i];
			}
			else if (args[i].equals("-w") || args[i].equalsIgnoreCase("-window")) {
				windowSizeForR2 = Integer.parseInt(args[++i]);
			}
			else if (args[i].equals("-s") || args[i].equalsIgnoreCase("-nmsnp")) {
				minNonMissingProportionForSNP = Double.parseDouble(args[++i]);
			}
			else if (args[i].equals("-t") || args[i].equalsIgnoreCase("-nmtaxon")) {
				minNonMissingProportionForTaxon = Double.parseDouble(args[++i]);
			}
			else if (args[i].equals("-r") || args[i].equalsIgnoreCase("-r2file")) {
				avgr2Filename = args[++i];
			}
			else if (args[i].equals("-x") || args[i].equalsIgnoreCase("-r2xyplot")) {
				avgr2Plotname = args[++i];
				if (!avgr2Plotname.endsWith(".png")) avgr2Plotname += ".png";
			}
			else if (args[i].equals("-c") || args[i].equalsIgnoreCase("-confile")) {
				propNonconsensusFilename = args[++i];
			}
			else if (args[i].equals("-s") || args[i].equalsIgnoreCase("-summaryfile")) {
				propNonconsensusFilename = args[++i];
			}
			else if (args[i].equals("-l") || args[i].equalsIgnoreCase("-logfile")) {
				logfileName = args[++i];
			}
			else if (args[i].equals("?")) myLogger.error(getUsage());
		}
	}

	public String getUsage() {
		StringBuilder usage = new StringBuilder("The QualityChecksPlugin can take the following parameters:\n");
		usage.append("-p or -pedigrees : a file containing pedigrees of the individuals to be evaluated\n");
		usage.append("-w or -window : use a window of +/- this size to evaluate average R-square (default = 25).\n");
		usage.append("-s or -nmsnp : the minimum proportion of non-missing values allowed for a snp (default = 0.1)\n");
		usage.append("-t or -nmtaxon : the minimum proportion of non-missing values allowed for a taxon (default = 0.1)\n");
		usage.append("-r or -r2file : the name of the file to save the average R2 value for each SNP\n");
		usage.append("-x or -r2xyplot : name of the png file of the average R2 of each SNP, .png will be appended\n");
		usage.append("-c or -confile : name of the file to save the proportion of nonConsensus SNPs for each taxon\n");
		usage.append("-s or -summaryfile : name of the file for the taxa summary\n");
		usage.append("-l or -logfile: name of the file to which the log will be appended\n");
		usage.append("? : print the parameter list.\n");

		return usage.toString();
	}
	
	@Override
	public ImageIcon getIcon() {
		return null;
	}

	@Override
	public String getButtonName() {
		return null;
	}

	@Override
	public String getToolTipText() {
		return null;
	}

	public void setWindowSizeForR2(int windowSizeForR2) {
		this.windowSizeForR2 = windowSizeForR2;
	}

	public void setMinNonMissingProportionForTaxon(
			double minNonMissingProportionForTaxon) {
		this.minNonMissingProportionForTaxon = minNonMissingProportionForTaxon;
	}

	public void setMinNonMissingProportionForSNP(
			double minNonMissingProportionForSNP) {
		this.minNonMissingProportionForSNP = minNonMissingProportionForSNP;
	}

	public void setAvgr2Filename(String avgr2Filename) {
		this.avgr2Filename = avgr2Filename;
	}

	public void setAvgr2Plotname(String avgr2Plotname) {
		this.avgr2Plotname = avgr2Plotname;
	}

	public void setPropNonconsensusFilename(String propNonconsensusFilename) {
		this.propNonconsensusFilename = propNonconsensusFilename;
	}

	public void setPedigreeFile(String pedigreeFile) {
		this.pedigreeFile = pedigreeFile;
	}

	public void setLogfileName(String logfileName) {
		this.logfileName = logfileName;
	}

	public void setSummaryFilename(String summaryFilename) {
		this.summaryFilename = summaryFilename;
	}

	
}
