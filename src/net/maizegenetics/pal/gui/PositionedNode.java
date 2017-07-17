// PositionedNode.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.pal.gui;

import net.maizegenetics.pal.tree.Node;
import net.maizegenetics.pal.tree.SimpleNode;


/**
 * A tree node that has a scalar position for the purposes of drawing the tree.
 *
 * @author Alexei Drummond
 * @version $Id: PositionedNode.java,v 1.1 2007/01/12 03:26:14 tcasstevens Exp $
 */
public class PositionedNode extends SimpleNode {

	protected double x; //Please excuse this - it will be returned to its non public state eventually - MG

	boolean highlight_;
	Node peer_;

	/** Builds a tree based on node, but highlights highlightNode */
	public PositionedNode(Node node, Node highlightNode) {

		init(node);
		this.peer_ = node;
		if (!node.isLeaf()) {
			for (int i = 0; i < node.getChildCount(); i++) {
				addChild(new PositionedNode(node.getChild(i),highlightNode));
			}
		}

		highlight_ = (node==highlightNode);
	}

	public PositionedNode(Node node) {

		init(node);
		this.peer_ = node;

		if (!node.isLeaf()) {
			for (int i = 0; i < node.getChildCount(); i++) {
				addChild(new PositionedNode(node.getChild(i)));
			}
		}
	}

	public void calculatePositions() {

		double[] currentXPos = {0.0};
		calculateXPositions(currentXPos);
	}
	public Node getPeer() {
		return peer_;
	}
	private double calculateXPositions(double[] currentXPos) {

		if (!isLeaf()) {
			// find average x position
			x = ((PositionedNode)getChild(0)).calculateXPositions(currentXPos);
			for (int i = 1; i < getChildCount(); i++) {
				x += ((PositionedNode)getChild(i)).calculateXPositions(currentXPos);
			}
			x /= getChildCount();
		} else {
			x = currentXPos[0];
			currentXPos[0] += 1.0;
		}

		return x;
	}

	public boolean isHighlighted() {
			return highlight_;
	}

	public double getX() {
			return x;
	}
}
