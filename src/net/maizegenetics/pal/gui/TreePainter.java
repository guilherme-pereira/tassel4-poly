// TreePainter.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.pal.gui;

import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.LabelMapping;
import net.maizegenetics.pal.tree.*;

import java.awt.*;

/**
 * A class that can paint a tree into a Graphics object.
 *
 * @version $Id: TreePainter.java,v 1.2 2007/01/21 02:52:14 buckler4 Exp $
 *
 * @author Alexei Drummond
 * @note
 *     <ul>
 *       <li> 14 August 2003 - Changed to reflect NameColouriser changes
 *     </ul>
 */
abstract public class TreePainter implements Painter {

	public static final Color BACKGROUND = Color.white;
	public static final Color FOREGROUND = Color.black;
	public static final Color NORMAL_LABEL_COLOR = Color.black;
	public static final LabelDisplayer NORMAL_LABEL_DISPLAY = LabelDisplayer.Utils.buildDisplay(Color.black);
	public static final LabelDisplayer HILITED_LABEL_DISPLAY = LabelDisplayer.Utils.buildDisplay(Color.red.darker(), Font.BOLD);
	public static final Color BOOTSTRAP_SUPPORT_COLOUR = Color.black;

	public static final String BOOTSTRAP_ATTRIBUTE_NAME = "bootstrap";

	public static final int DEFAULT_FONT_SIZE = 15;
	public static final int DEFAULT_FONT_STYLE = Font.PLAIN;
	public static final String DEFAULT_FONT_NAME = "times";
	public static final Font DEFAULT_FONT = new Font(DEFAULT_FONT_NAME, DEFAULT_FONT_STYLE, DEFAULT_FONT_SIZE);


	public PositionedNode treeNode;

	public String title_;

	protected String attName = null;

	boolean showTitle_;

	/** The tree being painted */
	private Tree tree;

	/** The time order character data used for determining symbols. */
	private TimeOrderCharacterData tocd = null;

	/** the number of leaves in the tree. */
	double width;

	/** the height of the root */
	double height;

	/**
	 * the maximum height of the display area in the same units as tree.
	 * if this is not -1.0 then it overrides the natural height.
	 */
	double maxHeight = -1.0;

	double maxLeafTime = 0.0;
	double sizeOfScale = 0.0;

	/** Width of pen used to paint lines */
	private int penWidth = 2;

	/** determines whether colors are used to distinguish branch depth */
	private boolean usingColor = true;

	/** determines whether node heights are displayed on the tree */
	private boolean showingNodeHeights = false;

	/** determines whether internal nodes are labelled */
	protected boolean showingInternalLabels = true;

	/** determines whether symbols are used instead of names */
	private boolean usingSymbols = false;

	private NameColouriser colouriser_ = null;
	private LabelMapping labelMapping_ = null;

	Font labelFont_ = DEFAULT_FONT;

	/**
	 * @param toDisplay the tree being painted.
	 * @param title the title of the tree.
	 * @param showTitle true if a title is being displayed.
	 * @param nc A name colouriser, may be null.
	 */
	public TreePainter(Tree toDisplay, String title, boolean showTitle) {

		this.title_ = title;
		this.showTitle_ = showTitle;
		this.tree = toDisplay;
//		if (toDisplay instanceof DatedTipsClockTree) {
//			tocd = ((DatedTipsClockTree)toDisplay).getTimeOrderCharacterData();
//		}
		standardTreePrep();
	}
	/**
	 * Returns -1 if no
	 * trap value available
	 */
	protected int getBootstrapValue(PositionedNode node){
		Object o = tree.getAttribute(node.getPeer(),BOOTSTRAP_ATTRIBUTE_NAME);
		if(o==null) {
			return -1;
		}
		return ((Integer)o).intValue();
	}

	/**
	 * Rotates the tree by leaf count, creates a positioned node version of the
	 * trees root and calculates postions and width and height information.
	 */
	protected void standardTreePrep() {
		TreeUtils.rotateByLeafCount(tree);
		treeNode = new PositionedNode(tree.getRoot());
		treeNode.calculatePositions();

		width = NodeUtils.getLeafCount(treeNode);
		height = treeNode.getNodeHeight();

		maxLeafTime = 0.0;
		maxLeafTime = getMaxLeafTime(treeNode);
		maxLeafTime *= 1.5;

		sizeOfScale = getSizeOfScale( height / 5.0);
	}

	/**
	 * sets the maximum height of plot.
	 * if this height is smaller than root height then
	 * only a proportion of tree is drawn.
	 */
	public final void setMaxHeight(double maxHeight) {
		this.maxHeight = maxHeight;
		sizeOfScale = getSizeOfScale( maxHeight / 5.0);
	}

	public final void setAttributeName(String name) {
		attName = name;
	}

	public final void setPenWidth(int p) {
		penWidth = p;
	}

	public final int getPenWidth() {
		return penWidth;
	}

	public final void setTree(Tree tree) {
		this.tree = tree;
		standardTreePrep();
		setTreeImpl(tree);
	}

	/**
	 * may be implemented by sub classes
	 */
	public void setTreeImpl(Tree tree) { }
	public final void setUsingColor(boolean use) {
		usingColor = use;
	}

	public final boolean isUsingColor() {
		return usingColor;
	}

	public final void setShowingNodeHeights(boolean s) {
		showingNodeHeights = s;
	}

	public final boolean isShowingNodeHeights() {
		return showingNodeHeights;
	}

	public final boolean isShowingInternalLabels() {
		return showingInternalLabels;
	}

	public final TimeOrderCharacterData getTimeOrderCharacterData() {
		return this.tocd;
	}

	public final void setTimeOrderCharacterData(TimeOrderCharacterData tocd) {
		this.tocd = tocd;
		usingSymbols = true;
	}

	public final boolean isUsingSymbols() {
		return usingSymbols;
	}

	/**
	 * Sets whether the tree is painted with symbols. This can
	 * only be set to true of a TimeOrderCharacterData has been set.
	 */
	public final void setUsingSymbols(boolean use) {
		usingSymbols = use;
		if (tocd == null) usingSymbols = false;
	}
	protected final Tree getTree() {
		return tree;
	}
	protected final double getSizeOfScale(double target) {

		double sos = 0.1;
		boolean accept = false;
		boolean divideByTwo = true;

		while (!accept) {
			if ((sos / target) >= 5.0) {
				sos /= (divideByTwo ? 2.0 : 5.0);
				divideByTwo = !divideByTwo;
			} else if ((sos / target) < 0.2) {
				sos *= (divideByTwo ? 5.0 : 2.0);
				divideByTwo = !divideByTwo;
			} else accept = true;
		}

		return sos;
	}


	protected static final double getMaxLeafTime(Node node) {

		if (!node.isLeaf()) {
		double max = getMaxLeafTime(node.getChild(0));
			double posmax = 0.0;
			for (int i = 1; i < node.getChildCount(); i++) {
				posmax = getMaxLeafTime(node.getChild(i));
				if (posmax > max) max = posmax;
			}
			return max;
		} else {
			return node.getNodeHeight();
		}
	}

	public final static void drawSymbol(Graphics g, int x, int y, int width, int index) {

		int halfWidth = width / 2;

		switch (index % 6) {
			case 0: g.fillRect(x, y, width, width); break;
			case 1: g.drawRect(x, y, width, width); break;
			case 2: g.fillOval(x, y, width, width); break;
			case 3: g.drawOval(x, y, width, width); break;
			case 4: // draw triangle
				g.drawLine(x, y + width, x + halfWidth, y);
				g.drawLine(x + halfWidth, y, x + width, y + width);
				g.drawLine(x, y + width, x + width, y + width);
				break;
			case 5: // draw X
				g.drawLine(x, y, x + width, y + width);
				g.drawLine(x, y + width, x + width, y);
				break;
		}
	}

	public final boolean isShowTitle() {
		return showTitle_;
	}
	public final void setColouriser(NameColouriser nc) {
		this.colouriser_ = nc;
	}
	public final void setLabelMapping(LabelMapping lp) {
		this.labelMapping_ = lp;
	}
	public final void setTitle(String title) {
		this.title_ = title;
		showTitle_ = true;
	}
	protected final String getNodeName(Node node) {
		if(labelMapping_!=null) {
			return labelMapping_.getLabel(node.getIdentifier());
		}
		return node.getIdentifier().getName();
	}
	public final String getTitle() {
		return title_;
	}

	public final void doTitle(Graphics g, int x, int y) {
		if(showTitle_) {
			g.drawString(title_, x,y);
		}
	}

	protected final void doScale(Graphics g, double xScale, int x, int y) {
		g.setColor(FOREGROUND);
		g.drawLine(x,
				 y,
				 x+ (int)Math.round(sizeOfScale * xScale),
				 y);

		g.drawString(sizeOfScale + " "+Units.UNIT_NAMES[tree.getUnits()],
				 x,y-12);
	}

	protected final LabelDisplayer getNodeDisplay(Node n, LabelDisplayer defaultDisplay) {
		if(colouriser_!=null) {
			return colouriser_.getDisplay(n.getIdentifier(),defaultDisplay);
		}
		return defaultDisplay;
	}
	private int getLongestIdentifierPixelWidth(FontMetrics fm, Node n) {
		Identifier id = n.getIdentifier();
		int myWidth = 0;
		if(id!=null) {
			String name = id.getName();
			if(name!=null) {
				if(labelMapping_!=null) {
					name = labelMapping_.getLabel(name,name);
				}
				myWidth = fm.stringWidth(name);
			}
		}
		int numberChild = n.getChildCount();
		for(int i = 0 ; i < numberChild ; i++) {
			myWidth = Math.max(myWidth, getLongestIdentifierPixelWidth(fm, n.getChild(i)));
		}
		return myWidth;
	}
	protected final int getLongestIdentifierPixelWidth(FontMetrics fm) {
		return getLongestIdentifierPixelWidth(fm,treeNode);
	}
// ============================================================================
// == Font Stuff
/**
 * Set the font used to display labels
 */
	public final void setLabelFont(Font f) { this.labelFont_ = f; }
	/**
	 * Set the font used to display labels
	 */
	public final void setLabelFontSize(int size)
        { this.labelFont_ = new Font(labelFont_.getFontName(),labelFont_.getStyle(), size); }
	protected final Font getLabelFont() {	return labelFont_;	}

}
