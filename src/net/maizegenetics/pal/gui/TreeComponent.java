// TreeComponent.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.gui;

import net.maizegenetics.pal.alignment.Phenotype;
import net.maizegenetics.pal.ids.LabelMapping;
import net.maizegenetics.pal.tree.TimeOrderCharacterData;
import net.maizegenetics.pal.tree.Tree;

import java.awt.*;

/**
 * An AWT Component for displaying a tree.
 *
 * @author Alexei Drummond
 * @version $Id: TreeComponent.java,v 1.2 2007/01/21 02:52:14 buckler4 Exp $
 */
public class TreeComponent extends Component {

    boolean circular_ = false;
    boolean isQuantitative = false;
    public final static int NORMAL_BW = 0;
    public final static int CIRCULAR_BW = 1;
    public final static int NORMAL_COLOR = 2;
    public final static int CIRCULAR_COLOR = 3;
    public final static int QUANTITATIVE = 4;
    int mode_;
    public final static String[] MODE_NAMES = new String[5];

    // unfortunate but necessary to avoid Java 1.1 language features


    static {
        MODE_NAMES[NORMAL_BW] = "Normal (bw)";
        MODE_NAMES[CIRCULAR_BW] = "Circular (bw)";
        MODE_NAMES[NORMAL_COLOR] = "Normal (color)";
        MODE_NAMES[CIRCULAR_COLOR] = "Circular (color)";
        MODE_NAMES[QUANTITATIVE] = "Quantitative";

    }
    TreePainterCircular circlePainter_;
    TreePainter painter_;
    TreePainterQuant quantitativePainter;
    boolean invertCiruclar_;

    // constructors
    public TreeComponent(Tree tree, boolean usingSymbols) {
        this(tree, (TimeOrderCharacterData) null, usingSymbols);
    }

    public TreeComponent(Tree tree, TimeOrderCharacterData tocd, boolean usingSymbols) {
        this(tree);

        if (tocd != null) {
            painter_.setTimeOrderCharacterData(tocd);
        }
        painter_.setUsingSymbols(usingSymbols);
    }

    public TreeComponent(Tree tree) {
        this(tree, "", false);
    }

    public TreeComponent(Tree tree, String title) {
        this(tree, title, true);
    }

    public TreeComponent(Tree tree, String title, boolean showTitle) {
        painter_ = new TreePainterNormal(tree, title, showTitle);
        circlePainter_ = new TreePainterCircular(tree, title, showTitle);
        quantitativePainter = new TreePainterQuant(tree, title, showTitle);
        setMode(NORMAL_COLOR);
        setSize(getPreferredSize());
    }

    public void setLabelFontSize(int size) {
        painter_.setLabelFontSize(size);
        circlePainter_.setLabelFontSize(size);
        quantitativePainter.setLabelFontSize(size);
    }

    public final void setLabelMapping(LabelMapping lm) {
        painter_.setLabelMapping(lm);
        circlePainter_.setLabelMapping(lm);
    }

    public void setColouriser(NameColouriser nc) {
        painter_.setColouriser(nc);
        circlePainter_.setColouriser(nc);
    }

    public void setAttributeName(String name) {
        painter_.setAttributeName(name);
        repaint();
    }

    public void setMaxHeight(double maxHeight) {
        painter_.setMaxHeight(maxHeight);
    }

    public void setTree(Tree tree) {
        painter_.setTree(tree);
        circlePainter_.setTree(tree);
        quantitativePainter.setTree(tree);
    }

    public void setTrait(Phenotype ca, int trait) {
        quantitativePainter.setTrait(ca, trait);
    }

    public void setTitle(String title) {
        painter_.setTitle(title);
        circlePainter_.setTitle(title);
    }

    public Dimension getPreferredSize() {
        return painter_.getPreferredSize();
    }

    public void setInvertCircular(boolean invert) {
        this.invertCiruclar_ = invert;
    }

    public Dimension getMinimumSize() {
        return getPreferredSize();
    }

    public String getTitle() {
        return painter_.getTitle();
    }

    /**
     * Sets the mode of display for this component.
     * @param mode to switch to. Valid arguments are NORMAL (for normal tree
     * display), and CIRCULAR (for a circular view of the trees)
     */
    public void setMode(int mode) {
        this.mode_ = mode;
        switch (mode) {
            case NORMAL_BW: {
                circular_ = false;
                painter_.setUsingColor(false);
                isQuantitative = false;
                break;
            }
            case NORMAL_COLOR: {
                circular_ = false;
                painter_.setUsingColor(true);
                isQuantitative = false;
                break;
            }
            case CIRCULAR_BW: {
                circular_ = true;
                circlePainter_.setUsingColor(false);
                isQuantitative = false;
                break;
            }
            case CIRCULAR_COLOR: {
                circular_ = true;
                circlePainter_.setUsingColor(true);
                isQuantitative = false;
                break;
            }
            case QUANTITATIVE: {
                isQuantitative = true;
                quantitativePainter.setUsingColor(false);
                break;
            }
        }
    }

    public void paint(Graphics g) {
        if (isQuantitative) {
            quantitativePainter.paint(g, getSize().width, getSize().height);
        } else if (!circular_) {
            painter_.paint(g, getSize().width, getSize().height);
        } else {
            circlePainter_.paint(g, getSize().width, getSize().height, invertCiruclar_);
        }

    }
}
