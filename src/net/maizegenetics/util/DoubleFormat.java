package net.maizegenetics.util;

import java.text.DecimalFormat;
import java.text.NumberFormat;

public class DoubleFormat {

    private static NumberFormat nfe = null;
    private static NumberFormat nf = null;
    private static String NAN_STRING = "NaN";

    public static String format(double value) {
        if (Double.isNaN(value)) {
            return NAN_STRING;
        }
        if (Double.isInfinite(value)) {
            return "Infinity";
        }
        if (nfe == null) {
            nfe = NumberFormat.getInstance();
            if (nfe instanceof DecimalFormat) {
                ((DecimalFormat) nfe).applyPattern("0.0###E0");
            }
        }
        if (nf == null) {
            nf = NumberFormat.getInstance();
            if (nf instanceof DecimalFormat) {
                ((DecimalFormat) nf).applyPattern("0.#####");
            }
        }
        if (value == 0) {
            return nf.format(value);
        }
        if (value < .001 || value >= 1000000) {
            return nfe.format(value);
        }
        return nf.format(value);
    }

    public static String format(Double value) {
        if (value == null) {
            return NAN_STRING;
        }
        return format(value.doubleValue());
    }
}
