package net.maizegenetics.gwas.NAM;

import java.io.File;

public class FileNames {
	public File snps = null;
	public File chrmodel = null;
	public File chrsteps = null;
	public File model = null;
	public File steps = null;
	public File residuals = null;
	public File agpmap = null;
	public File namMarkersByChr = null;
	public File namMarkers = null;
	public File phenotypes = null;
	public File projectedFile = null;
	
	public double enterlimit = 1e-6;
	public double exitlimit = 2e-6;
	public double[] limits = null;
	public int traitnumber = -1;
	public int iterations = 100;
	public String replacement;
	public String analysis = "";
	public int maxsnps = 1000;
	public boolean threaded = false;
	public boolean permute = false;
	public boolean bootstrapPermutation = false;
	public boolean randomizeSnpOrder = false;
	public boolean fullModel = false;
	public int startIteration = 1;
	public int chromosome = 0;
	public boolean useB73asReference = true;
	
}
