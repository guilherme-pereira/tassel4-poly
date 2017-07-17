// TreePainterNormal.java
//


// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package net.maizegenetics.pal.gui;

import net.maizegenetics.pal.io.FormattedOutput;
import net.maizegenetics.pal.tree.Tree;

import java.awt.*;

/**
 * A class that can paint a tree into a Graphics object.
 *
 * @version $Id: TreePainterNormal.java,v 1.2 2007/01/21 02:52:14 buckler4 Exp $
 *
 * @author Alexei Drummond
 * @note
 *     <ul>
 *       <li> 14 August 2003 - Changed to reflect NameColouriser changes
 *     </ul>
 */
public class TreePainterNormal extends TreePainter {

	public static final int RIGHTBORDER = 10;
	public static final int LEFTBORDER = 10;
	public static final int TOPBORDER = 20;
	public static final int BOTTOMBORDER = 30;

	public static final int FONT_SIZE = 15;
	public static final int YSPACER = 20;
	public static final int XSPACER = 4;

	private double xScale = 1.0;
	private double yScale = 1.0;

//	private Font labelFont_ = new Font("Times", Font.PLAIN, FONT_SIZE);
	int maxLabelWidth_ = -1;

	public TreePainterNormal(Tree toDisplay, String title, boolean showTitle) {
		super(toDisplay,title,showTitle);
	}

	/**
	 * Returns the preferred size for drawing
	 * (that is the size that will show everything nicely)
	 */
	public Dimension getPreferredSize() {
		return new Dimension(100 + LEFTBORDER + RIGHTBORDER,
			(int)Math.round(width * FONT_SIZE) + TOPBORDER + BOTTOMBORDER);
	}

	protected void paint(PositionedNode node, Graphics g,
		int displayWidth, int displayHeight, LayoutTracker lt, boolean isRoot) {

		Point p = getPoint(node,displayWidth, displayHeight);
		g.setColor(FOREGROUND);
		if(isRoot) {
			g.fillRect(p.x-4, p.y-1, 4, 3); //Cheap hack!
		}
		if (node.hasChildren()) {
			for (int i = 0; i < node.getChildCount(); i++) {
				paintLeafBranch(p, getPoint((PositionedNode)node.getChild(i),
					displayWidth,displayHeight), node, g,lt);
			}

			for (int i = 0; i < node.getChildCount(); i++) {
				paint((PositionedNode)node.getChild(i), g,displayWidth, displayHeight,lt,false);
			}
			int bootStrapValue = getBootstrapValue(node);
			if(bootStrapValue>=50) {
				g.setColor(BOOTSTRAP_SUPPORT_COLOUR);
				g.drawString(bootStrapValue+"", p.x + XSPACER,
							p.y + (FONT_SIZE / 2));
			}
		} else {

			if ((maxLeafTime > 0.0) && isUsingColor()) {
				g.setColor(Color.getHSBColor((float)(maxLeafTime - node.getNodeHeight())/(float)maxLeafTime, 1.0f, 1.0f));
			} else {
				g.setColor(NORMAL_LABEL_COLOR);
			}

			if (isUsingColor()) {
				int halfWidth = getPenWidth() / 2;
				g.fillRect(p.x - halfWidth, p.y - halfWidth, getPenWidth(), getPenWidth());
			}
			if (isUsingSymbols()&&getTimeOrderCharacterData()!=null) {

				drawSymbol(g, p.x + XSPACER, p.y - (FONT_SIZE / 2), FONT_SIZE,
					getTimeOrderCharacterData().getTimeOrdinal(getTimeOrderCharacterData().whichIdNumber(node.getIdentifier().getName())));
			} else {
				String name = getNodeName(node);
				int width = g.getFontMetrics().stringWidth(name);
				if(isUsingColor()) {
					g.drawString(name, p.x + XSPACER,
							p.y + (FONT_SIZE / 2));
					if(node.isHighlighted()) {
						g.setColor(Color.red);
						g.drawOval(p.x - 4+XSPACER, p.y-FONT_SIZE/2-5, width +10, FONT_SIZE+8 );
					}
				} else {
					LabelDisplayer defaultDisplay =(node.isHighlighted() ? HILITED_LABEL_DISPLAY : NORMAL_LABEL_DISPLAY );
					getNodeDisplay(node,defaultDisplay).display(g,name, p.x + XSPACER, p.y + (FONT_SIZE / 2));
				}
				//Inform layout tracker of new String
				if(lt!=null) {
					lt.addMapping(name,new Rectangle(p.x+XSPACER, p.y - (FONT_SIZE / 2), width,FONT_SIZE));
				}
			}
		}

	}

	public Point getPoint(PositionedNode node, int displayWidth, int displayHeight) {

		return new Point(displayWidth -
			(int)Math.round(node.getNodeHeight() * xScale) - RIGHTBORDER,
			(int)Math.round(node.x * yScale) + TOPBORDER);
	}

	private void paintLeafBranch(Point p, Point lp, PositionedNode node, Graphics g, LayoutTracker lt) {

		int halfWidth = getPenWidth() / 2;

		// paint join to parent
		g.fillRect(p.x - halfWidth, Math.min(p.y, lp.y) - halfWidth,
				getPenWidth(), Math.abs(lp.y - p.y) + getPenWidth());

		// paint branch
		g.fillRect(Math.min(p.x, lp.x) - halfWidth, lp.y - halfWidth,
				 Math.abs(lp.x - p.x) + getPenWidth(), getPenWidth());

		if (isShowingNodeHeights()) {

			String label = FormattedOutput.getInstance().getDecimalString(node.getNodeHeight(), 4);
			int width = g.getFontMetrics().stringWidth(label);

			int x = Math.min(p.x, lp.x) - (width / 2);

			g.drawString(label, x, p.y - halfWidth - 1);

		}

		if (isShowingInternalLabels()) {
			String label = getNodeName(node);
			int width = g.getFontMetrics().stringWidth(label);

			int x = Math.min(p.x, lp.x) - (width / 2);

			g.drawString(label, x, p.y - halfWidth - 1);
		}

		Object att = null;
		if (attName != null) {
			if (attName.equals("node height")) {
				att = new Double(node.getNodeHeight());
			} else if (attName.equals("branch length")) {
				att = new Double(node.getBranchLength());
			} else {
				att = node.getAttribute(attName);
			}
			if (att != null) {
				String label = null;
				if (att instanceof Double) {
					label = FormattedOutput.getInstance().getDecimalString(((Double)att).doubleValue(), 3);
				} else label = att.toString();

				int width = g.getFontMetrics().stringWidth(label);
				int height = g.getFontMetrics().getAscent();
				int x = Math.min(p.x, lp.x) + halfWidth + 1;
				g.drawString(label, x, p.y + (height / 2));
			}
		}
	}

	public void paint(Graphics g, int displayWidth, int displayHeight) {
		paint(g,displayWidth,displayHeight,null);
	}

	public void paint(Graphics g, int displayWidth, int displayHeight, LayoutTracker lt) {
		g.setFont(labelFont_);

		if(maxLabelWidth_<0) {
			maxLabelWidth_ = getLongestIdentifierPixelWidth(g.getFontMetrics());
		}

		double h = height;
		if (maxHeight != -1.0) { h = maxHeight; }
		xScale = (double)(displayWidth - LEFTBORDER - RIGHTBORDER  - maxLabelWidth_) / h;
		yScale = (double)(displayHeight - TOPBORDER - BOTTOMBORDER) / width;

		g.setColor(BACKGROUND);
		g.fillRect(0, 0, displayWidth, displayHeight);
		paint(treeNode, g, displayWidth-maxLabelWidth_, displayHeight,lt, true);

		doTitle(g,LEFTBORDER, TOPBORDER - 8);
		doScale(g,xScale,LEFTBORDER,displayHeight - BOTTOMBORDER + 12);
	}
}
