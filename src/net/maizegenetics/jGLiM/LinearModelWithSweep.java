package net.maizegenetics.jGLiM;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.DoubleFactory2D;

import java.util.ArrayList;
import java.util.Iterator;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: PeterLocal
 * Date: Nov 24, 2007
 * Time: 10:15:49 AM
 * To change this template use File | Settings | File Templates.
 */
public class LinearModelWithSweep {

    public DoubleMatrix2D[][] xtxmatrices;
    public DoubleMatrix2D[][] xtymatrices;
    double[] data;
    public Sweep sweep;
    public ArrayList<Double> effectSS = new ArrayList<Double>();  //incremental or type I SS for model effects
    public ArrayList<Double> effectdf = new ArrayList<Double>();  //incremental of type I df for model effects
    public double totalSS;

    public static void main(String[] args) {
        LinearModelWithSweep lm = new LinearModelWithSweep();
        lm.test();
    }

    public void test() {
        ArrayList<String> A = new ArrayList<String>();
        ArrayList<String> B = new ArrayList<String>();
        ArrayList<String> C = new ArrayList<String>();
        ArrayList<String> D = new ArrayList<String>();
        ArrayList<Double> values = new ArrayList<Double>();
        try {
            BufferedReader br = new BufferedReader(new FileReader("c:/temp/testdata.txt"));
            //br.readLine();
            String line;
            while ((line = br.readLine()) != null) {
                String[] splitline = line.split("\t");
                A.add(splitline[0]);
                B.add(splitline[1]);
                C.add(splitline[2]);
                D.add(splitline[3]);
                values.add(Double.parseDouble(splitline[4]));
            }
            br.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        data = new double[values.size()];
        Iterator<Double> dit = values.iterator();
        int count = 0;
        while (dit.hasNext()) {
            data[count++] = dit.next().doubleValue();
        }

        double[] dblC = new double[C.size()];
        Iterator<String> sit = C.iterator();
        count = 0;
        while (sit.hasNext()) {
            dblC[count++] = Double.parseDouble(sit.next());
        }

        double[] dblD = new double[D.size()];
        sit = D.iterator();
        count = 0;
        while (sit.hasNext()) {
            dblD[count++] = Double.parseDouble(sit.next());
        }

        int[] meanLevel = new int[data.length];
        ModelEffect memean = new ModelEffect(meanLevel);
        ModelEffect meA = new ModelEffect(ModelEffect.getIntegerLevels(A));
        ModelEffect meB = new ModelEffect(ModelEffect.getIntegerLevels(B));
        ModelEffect meC = new CovariateModelEffect(dblC);
        ModelEffect meD = new CovariateModelEffect(dblD);
        ModelEffect[] theEffects = new ModelEffect[]{memean, meA, meB, meC, meD};
        xtxmatrices = modelEffectsToxtx(theEffects);
        xtymatrices = modelEffectsToxty(theEffects, data);

        initialSweep();

        double cfm = effectSS.get(0);
        for (int i = 0; i < effectSS.size(); i++) {
            System.out.println(i + ", " + effectSS.get(i) + ", " + effectdf.get(i));
        }
        System.out.println("model ss = " + (getModelSS() - cfm) + ", model df = " + (getModeldf() - 1));
        System.out.println("error ss = " + getErrorSS() + ", error df = " + getErrordf());
        double[] result = marginalEffectSSdf(1);
        System.out.println("marginal A ss = " + result[0] + ", df = " + result[1]);
        result = marginalEffectSSdf(2);
        System.out.println("marginal B ss = " + result[0] + ", df = " + result[1]);
        result = marginalEffectSSdf(3);
        System.out.println("marginal C ss = " + result[0] + ", df = " + result[1]);
        result = marginalEffectSSdf(4);
        System.out.println("marginal D ss = " + result[0] + ", df = " + result[1]);

        System.out.println("finished");
    }

    public LinearModelWithSweep(DoubleMatrix2D[][] xtxmatrices, DoubleMatrix2D[][] xtymatrices, double[] data) {
        this.xtxmatrices = xtxmatrices;
        this.xtymatrices = xtymatrices;
        this.data = data;
        sweep = new Sweep();
        initialSweep();
    }

    public LinearModelWithSweep(ArrayList<ModelEffect> effects, double[] data) {
        int numberOfEffects = effects.size();
        this.data = data;
        xtxmatrices = new DoubleMatrix2D[numberOfEffects][numberOfEffects];
        xtymatrices = new DoubleMatrix2D[numberOfEffects][1];
        for (int i = 0; i < numberOfEffects; i++) {
            ModelEffect me = effects.get(i);
            xtxmatrices[i][i] = me.getXTX();
            xtymatrices[i][0] = DoubleFactory2D.dense.make(me.getXTy(data).toArray(), me.getNumberOfLevels());
            for (int j = 0; j < i; j++) {
                xtxmatrices[i][j] = ModelEffect.getX1TX2(me, effects.get(j));
                xtxmatrices[j][i] = xtxmatrices[i][j].viewDice();
            }
        }
        sweep = new Sweep();
        initialSweep();
    }

    public LinearModelWithSweep() {
        sweep = new Sweep();
    }

    private void initialSweep() {
        sweep.XTX = DoubleFactory2D.dense.compose(xtxmatrices);
        sweep.XTy = DoubleFactory2D.dense.compose(xtymatrices);

        totalSS = 0;
        for (int i = 0; i < data.length; i++) {
            totalSS += data[i] * data[i];
        }
        sweep.yTy = DoubleFactory2D.dense.make(1, 1, totalSS);
        sweep.makeA();
        int lastrow = sweep.A.rows() - 1;
        int col = 0;
        double prevSS = totalSS;

        for (int i = 0; i < xtxmatrices.length; i++) {
            double df = 0;
            for (int j = 0; j < xtxmatrices[i][i].rows(); j++) {
                if (sweep.revg2sweep(col++)) {
                    df++;
                }
            }
            effectdf.add(new Double(df));
            double newSS = sweep.A.getQuick(lastrow, lastrow);
            double difSS = prevSS - newSS;
            effectSS.add(new Double(difSS));
            prevSS = newSS;
        }
    }

    public double getTotalSS() {
        return totalSS - effectSS.get(0);
    }

    public double getModelSS() {
        return getTotalSS() - getErrorSS();
    }

    public double getErrorSS() {
        int last = sweep.A.rows() - 1;
        return sweep.A.getQuick(last, last);
    }

    public double getTotaldf() {
        return data.length - 1;
    }

    public double getErrordf() {
        return getTotaldf() - getModeldf();
    }

    public double getModeldf() {
        Iterator<Double> dit = effectdf.iterator();
        double df = 0;
        while (dit.hasNext()) {
            df += dit.next();
        }
        return df - 1;
    }

    public DoubleMatrix1D getBeta() {
        return sweep.A.viewColumn(sweep.A.columns() - 1).viewPart(0, sweep.A.rows() - 1).copy();
    }

    public ArrayList<Double> incrementalEffectSS() {
        return effectSS;
    }

    public ArrayList<Double> incrementalEffectdf() {
        return effectdf;
    }

    public double[] marginalEffectSSdf(int index) {
        int error = sweep.A.rows() - 1;
        double errorss = sweep.A.getQuick(error, error);

        int first = 0;
        for (int i = 0; i < index; i++) {
            first += xtxmatrices[i][i].rows();
        }
        int last = first + xtxmatrices[index][index].rows();
        int df = 0;
        for (int i = first; i < last; i++) {
            if (sweep.revg2sweep(i)) {
                df++;
            }
        }

        int nonsingular = sweep.sweepSingularColumns();
        df -= nonsingular;

        double reducedError = sweep.A.getQuick(error, error);
        double[] result = new double[]{(reducedError - errorss), df};
        for (int i = first; i < last; i++) {
            sweep.revg2sweep(i);
        }

        return result;
    }

    public DoubleMatrix2D getInverse() {
        int last = sweep.A.rows() - 1;
        return sweep.A.viewPart(0, 0, last, last).copy();
    }

    public static DoubleMatrix2D[][] modelEffectsToxtx(ModelEffect[] mes) {
        int dim = mes.length;
        DoubleMatrix2D[][] xtx = new DoubleMatrix2D[dim][dim];
        for (int i = 0; i < dim; i++) {
            xtx[i][i] = mes[i].getXTX();
            for (int j = i + 1; j < dim; j++) {
                xtx[i][j] = ModelEffect.getX1TX2(mes[i], mes[j]);
                xtx[j][i] = xtx[i][j].viewDice();
            }
        }
        return xtx;
    }

    public static DoubleMatrix2D[][] modelEffectsToxtx(ArrayList<ModelEffect> mes) {
        ModelEffect[] mearray = new ModelEffect[mes.size()];
        mes.toArray(mearray);
        return modelEffectsToxtx(mearray);
    }

    public static DoubleMatrix2D[][] modelEffectsToxty(ModelEffect[] mes, double[] data) {
        int dim = mes.length;
        DoubleMatrix2D[][] xty = new DoubleMatrix2D[dim][1];
        for (int i = 0; i < dim; i++) {
            DoubleMatrix1D dm = mes[i].getXTy(data);
            xty[i][0] = DoubleFactory2D.dense.make(dm.toArray(), dm.size());
        }
        return xty;
    }

    public static DoubleMatrix2D[][] modelEffectsToxty(ArrayList<ModelEffect> mes, double[] data) {
        ModelEffect[] mearray = new ModelEffect[mes.size()];
        mes.toArray(mearray);
        return modelEffectsToxty(mearray, data);
    }

    public double[] getFp() {
        double F = getModelSS() / getModeldf() / getErrorSS() * getErrordf();
        double p;
        try {
            p = AbstractLinearModel.Ftest(F, getModeldf(), getErrordf());
        } catch (Exception e) {
            p = 1.0;
        }
        return new double[]{F, p};
    }
}
