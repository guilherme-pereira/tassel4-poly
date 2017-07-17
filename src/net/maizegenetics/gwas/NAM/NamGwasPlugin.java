package net.maizegenetics.gwas.NAM;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;

public class NamGwasPlugin extends AbstractPlugin {
	private boolean resample = true;
	private FileNames parameters = new FileNames();
	
	private static final Logger myLogger = Logger.getLogger(NamGwasPlugin.class);

    public NamGwasPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

	@Override
	public DataSet performFunction(DataSet input) {
		if (parameters.agpmap == null || parameters.residuals == null || parameters.namMarkersByChr == null || parameters.snps == null || parameters.chrmodel == null || parameters.chrsteps == null) {
			myLogger.info(getUsage());
			return null;
		}
		
		try {
			//read the chromosome from the snp file
			BufferedReader br = new BufferedReader(new FileReader(parameters.snps));
			br.readLine();
			String[] data = br.readLine().split("\t");
			parameters.chromosome = Integer.parseInt(data[2]);
			br.close();
		} catch(IOException e) {
			myLogger.error(e);
		}
		
		ModelFitterBootstrapForward fitter = new ModelFitterBootstrapForward(parameters);
		fitter.run();
		return null;
	}

    public void setMapFilename(String mapFilename) {
		parameters.agpmap = new File(mapFilename);
	}

	public void setResidualFilename(String residualFilename) {
		parameters.residuals = new File(residualFilename);
	}

	public void setRilmarkerFilename(String rilmarkerFilename) {
		parameters.namMarkersByChr = new File(rilmarkerFilename);
	}

	public void setFounderFilename(String founderFilename) {
		parameters.snps = new File(founderFilename);
	}

	public void setRandomizeSnpOrder(boolean randomizeSnpOrder) {
		parameters.randomizeSnpOrder = randomizeSnpOrder;
	}

	public void setEnterLimit(double enterLimit) {
		parameters.enterlimit = enterLimit;
	}

	public void setResample(boolean resample) {
		this.resample = resample;
	}

	public void setIterations(int iterations) {
		parameters.iterations = iterations;
	}

	public void setStart(int start) {
		parameters.startIteration = start;
	}

	public void setThreaded(boolean threaded) {
		parameters.threaded = threaded;
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

	@Override
	public void setParameters(String[] args) {
		
		int narg = args.length;
		for (int i = 0; i < narg; i++) {
			if (args[i].equals("-c") || args[i].equalsIgnoreCase("-map")) {
				parameters.agpmap = new File(args[++i]);
			}
			else if (args[i].equals("-t") || args[i].equalsIgnoreCase("-trait")) {
				parameters.residuals = new File(args[++i]);
			}
			else if (args[i].equals("-r") || args[i].equalsIgnoreCase("-rils")) {
				parameters.namMarkersByChr = new File(args[++i]);
			}
			else if (args[i].equals("-f") || args[i].equalsIgnoreCase("-founders")) {
				parameters.snps = new File(args[++i]);
			}
			else if (args[i].equals("-m") || args[i].equalsIgnoreCase("-model")) {
				parameters.chrmodel = new File(args[++i]);
			}
			else if (args[i].equals("-s") || args[i].equalsIgnoreCase("-steps")) {
				parameters.chrsteps = new File(args[++i]);
			}
			else if (args[i].equals("-a") || args[i].equalsIgnoreCase("-randomize")) {
				if (args[++i].toUpperCase().startsWith("T")) parameters.randomizeSnpOrder = true;
				else parameters.randomizeSnpOrder = false;
			}
			else if (args[i].equals("-e") || args[i].equalsIgnoreCase("-enterlimit")) {
				parameters.enterlimit = Double.parseDouble(args[++i]);
			}
			else if (args[i].equals("-i") || args[i].equalsIgnoreCase("-iterations")) {
				parameters.iterations = Integer.parseInt(args[++i]);
			}
			else if (args[i].equals("-d") || args[i].equalsIgnoreCase("-start")) {
				parameters.startIteration = Integer.parseInt(args[++i]);
			}
			else if (args[i].equals("-g") || args[i].equalsIgnoreCase("-maxsnps")) {
				parameters.maxsnps = Integer.parseInt(args[++i]);
			}
			else if (args[i].equals("-noresample")) {
				resample = false;
			}
			else if (args[i].equals("-enablethreads")) {
				parameters.threaded = true;
			}
			else if (args[i].equals("-fullmodel")) {
				parameters.fullModel = true;
			}
			else if (args[i].equals("-noref")) { //does not use B73 as reference, assumes alleles coded as 0 are from B73
				parameters.useB73asReference = false;
			}
			else if (args[i].equals("?")) myLogger.info(getUsage());
			else {
				myLogger.info(getUsage());
				return;
			}
		}
	}
	
	private String getUsage() {
		StringBuilder usage = new StringBuilder("The NamGwasPlugin takes the following parameters:\n");
		usage.append("-c or -map : file name of RIL markers with map coordinates (required)\n");
		usage.append("-t or -trait : file name of chromosome residuals for a trait (required)\n");
		usage.append("-r or -rils : file name of the ril markers for this chromosome (required)\n");
		usage.append("-f or -founders : file name of markers for this chromosome genotyped on founders to be projected on RILs (required)\n");
		usage.append("-m or -model : file name of output file for the final fitted model (required)\n");
		usage.append("-s or -steps : file name of output file for the model fitting steps as markers are added (required)\n");
		usage.append("-a or -randomize : true if snps should be tested in random order (default = false)\n");
		usage.append("-e or -enterlimit : the largest p-value for which a new term will be added to the model (default = 1e-6)\n");
		usage.append("-i or -iterations : the number of resampling iterations (default = 100)\n");
		usage.append("-d or -start : the number of the first iteration in this sequence (default = 1)\n");
		usage.append("-g or -maxsnps : the maximum number of snps that will be fit (default = 1000)\n");
//		usage.append("-noresample : do not resample (default = resample)\n");
//		usage.append("-enablethreads : have application use multiple cores if available. (default is single threaded.)\n");
		usage.append("-fullmodel : test snps for entry using the full model (default = use residuals from the previous model)\n");
		usage.append("? : print the parameter list.\n");
		
		return usage.toString();
	}
	
	 
}
