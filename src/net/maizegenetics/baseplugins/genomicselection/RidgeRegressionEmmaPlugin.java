package net.maizegenetics.baseplugins.genomicselection;

import java.awt.Frame;
import java.net.URL;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import javax.swing.ImageIcon;
import javax.swing.JOptionPane;

import org.apache.log4j.Logger;

import net.maizegenetics.jGLiM.dm.FactorModelEffect;
import net.maizegenetics.jGLiM.dm.ModelEffectUtils;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.pal.alignment.MarkerPhenotype;
import net.maizegenetics.pal.alignment.MarkerPhenotypeAdapter;
import net.maizegenetics.pal.alignment.MarkerPhenotypeAdapterUtils;
import net.maizegenetics.pal.alignment.Phenotype;
import net.maizegenetics.pal.alignment.SimplePhenotype;
import net.maizegenetics.pal.alignment.Trait;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.report.SimpleTableReport;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

public class RidgeRegressionEmmaPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(RidgeRegressionEmmaPlugin.class);

    public RidgeRegressionEmmaPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet performFunction(DataSet input) {
        try {
            List<Datum> datasets = input.getDataOfType(Phenotype.class);
            if (datasets.size() < 1) {
                String msg = "No datasets of an appropriate type were selected for the GS analysis.";
                myLogger.error(msg);
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), msg, "GS Error", JOptionPane.ERROR_MESSAGE);
                }
                return null;
            }

            LinkedList<Datum> results = new LinkedList<Datum>();
            for (Datum dataset : datasets) {
                try {
                    LinkedList<Datum> aResult = null;
                    aResult = processData(dataset);
                    if (aResult != null) {
                        results.addAll(aResult);
                        fireDataSetReturned(new DataSet(aResult, this));
                    }
                } catch (Exception e) {
                    StringBuilder msg = new StringBuilder("Error in GS processing " + dataset.getName());
                    msg.append(". ").append(e.getMessage());
                    myLogger.error(msg.toString());
                    e.printStackTrace();
                    if (isInteractive()) {
                        JOptionPane.showMessageDialog(getParentFrame(), msg.toString(), "GS Error", JOptionPane.ERROR_MESSAGE);
                    }
                }
            }

            return new DataSet(results, this);
        } finally {
            fireProgress(100);
        }
    }

    public LinkedList<Datum> processData(Datum dataset) {
        DoubleMatrix phenotype;
        DoubleMatrix genotype;
        DoubleMatrix fixedEffects;
        LinkedList<Datum> theResults = new LinkedList<Datum>();

        MarkerPhenotypeAdapter theAdapter;
        if (dataset.getDataType().equals(MarkerPhenotype.class)) {
            String msg = "Ridge Regression has only been implemented for numeric genotypes. No analysis will be run on this data.";
            if (isInteractive()) {
                JOptionPane.showMessageDialog(getParentFrame(), msg, "Ridge Regression Error", JOptionPane.ERROR_MESSAGE);
            } else {
                myLogger.error(msg);
            }
            return null;
        } else {
            theAdapter = new MarkerPhenotypeAdapter((Phenotype) dataset.getData());
        }

        //numbers of different things
        int numberOfMarkers = theAdapter.getNumberOfMarkers();
        int numberOfPhenotypes = theAdapter.getNumberOfPhenotypes();

        //iterate through the phenotypes
        for (int ph = 0; ph < numberOfPhenotypes; ph++) {

            //get phenotype data
            double[] phenotypeData = theAdapter.getPhenotypeValues(ph);
            int nObs = phenotypeData.length;
            phenotype = DoubleMatrixFactory.DEFAULT.make(nObs, 1, phenotypeData);

            //get the taxa
            Identifier[] taxaIDs = theAdapter.getTaxa(ph);

            //keep track of missing rows
            boolean[] missing = theAdapter.getMissingPhenotypes(ph);

            //get factors
            ArrayList<String[]> factorList = MarkerPhenotypeAdapterUtils.getFactorList(theAdapter, ph, missing);

            //get covariates
            ArrayList<double[]> covariateList = MarkerPhenotypeAdapterUtils.getCovariateList(theAdapter, ph, missing);

            //make the fixed effect matrix
            int numberOfFactors;
            if (factorList == null) {
                numberOfFactors = 0;
            } else {
                numberOfFactors = factorList.size();
            }
            int numberOfCovariates;
            if (covariateList == null) {
                numberOfCovariates = 0;
            } else {
                numberOfCovariates = covariateList.size();
            }
            int numberOfEffects = numberOfFactors + numberOfCovariates + 1;
            if (numberOfEffects > 1) {
                DoubleMatrix[][] effects = new DoubleMatrix[1][numberOfEffects];
                effects[0][0] = DoubleMatrixFactory.DEFAULT.make(nObs, 1, 1);
                for (int i = 0; i < numberOfFactors; i++) {
                    FactorModelEffect fme = new FactorModelEffect(ModelEffectUtils.getIntegerLevels(factorList.get(i)), true);
                    effects[0][i + 1] = fme.getX();
                }
                for (int i = 0; i < numberOfCovariates; i++) {
                    effects[0][i + numberOfFactors + 1] = DoubleMatrixFactory.DEFAULT.make(nObs, 1, covariateList.get(i));
                }
                fixedEffects = DoubleMatrixFactory.DEFAULT.compose(effects);
            } else {
                fixedEffects = DoubleMatrixFactory.DEFAULT.make(nObs, 1, 1);
            }

            //get the genotypes
            genotype = DoubleMatrixFactory.DEFAULT.make(nObs, numberOfMarkers);
            String[] markerNames = new String[numberOfMarkers];
            for (int m = 0; m < numberOfMarkers; m++) {
                Object[] markerValue = theAdapter.getMarkerValue(ph, m);
                for (int i = 0; i < nObs; i++) {
                    genotype.set(i, m, ((Double) markerValue[i]).doubleValue());
                }
                markerNames[m] = theAdapter.getMarkerName(m);
            }

            RegRidgeEmmaDoubleMatrix ridgeRegression = new RegRidgeEmmaDoubleMatrix(phenotype, fixedEffects, genotype);
            ridgeRegression.solve();

            //output the gebv's for the taxa
            double[] gebv = ridgeRegression.getBlups();
            double[][] traitTable = new double[nObs][1];
            for (int i = 0; i < nObs; i++) {
                traitTable[i][0] = gebv[i];
            }
            LinkedList<Trait> traitList = new LinkedList<Trait>();
            String phenoName = theAdapter.getPhenotypeName(ph);
            traitList.add(new Trait(phenoName + "_GEBV", false, Trait.TYPE_DATA));
            Phenotype outGebv = new SimplePhenotype(new SimpleIdGroup(taxaIDs), traitList, traitTable);
            String datumName = dataset.getName() + "_GEBVs_" + phenoName;
            StringBuilder comment = new StringBuilder("Ridge Regression from ");
            comment.append(dataset.getName()).append(":\n");
            comment.append("Genomic Estimated Breeding Values (GEBVs)\n");
            comment.append("trait = ").append(phenoName).append("\n");
            comment.append(nObs).append(" lines");
            theResults.add(new Datum(datumName, outGebv, comment.toString()));

            //output the marker blups for the markers as a report
            double[] markerBlups = ridgeRegression.getMrkBlups();
            Object[][] blupTable = new Object[numberOfMarkers][2];
            for (int i = 0; i < numberOfMarkers; i++) {
                blupTable[i][0] = markerNames[i];
                blupTable[i][1] = new Double(markerBlups[i]);
            }
            SimpleTableReport str = new SimpleTableReport("Marker BLUPs for " + dataset.getName(), new String[]{"Marker", phenoName + "_BLUP"}, blupTable);
            datumName = dataset.getName() + "_marker BLUPs_" + phenoName;
            comment = new StringBuilder("Ridge Regression from ");
            comment.append(dataset.getName()).append(":\n");
            comment.append("Marker BLUPs\n");
            comment.append("trait = ").append(phenoName).append("\n");
            comment.append(numberOfMarkers).append(" markers");
            theResults.add(new Datum(datumName, str, comment.toString()));

        }


        return theResults;
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = RidgeRegressionEmmaPlugin.class.getResource("LinearAssociation.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Genomic Selection";
    }

    @Override
    public String getToolTipText() {
        return "Predict Phenotypes using Ridge Regression for Genomic Selection";
    }
}
