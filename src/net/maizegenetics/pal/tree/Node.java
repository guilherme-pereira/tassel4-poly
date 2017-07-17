// Node.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package net.maizegenetics.pal.tree;

import net.maizegenetics.pal.ids.Identifier;

import java.io.Serializable;


/**
 * interface for a node (includes branch) in a binary/non-binary
 * rooted/unrooted tree
 *
 * @version $Id: Node.java,v 1.1 2007/01/12 03:26:17 tcasstevens Exp $
 *
 * @author Alexei Drummond
 * @author Korbinian Strimmer
 *
 */

public interface Node extends Serializable {

	/** Returns the parent node of this node. */
	Node getParent();

	/** Set the parent node of this node. */
	void setParent(Node node);

	/** Returns the sequence at this node, in the form an array of bytes. */
	byte[] getSequence();

	/** Sets the sequence using an array of bytes. */
	void setSequence(byte[] array);

	/** return the index of this node */
	int getNumber();

	/** set the index of this node */
	void setNumber(int number);

	/** Get the length of the branch attaching this node to its parent. */
	double getBranchLength();

	/**
	 * Set the length of the branch attaching this node to its parent.
	 */
	void setBranchLength(double value);

	/** Get the length SE of the branch attaching this node to its parent. */
	double getBranchLengthSE();

	/** Set the length SE of the branch attaching this node to its parent. */
	void setBranchLengthSE(double value);

	/** Get the height of this node relative to the most recent node. */
	double getNodeHeight();

	/**
	 * Set the height of this node relative to the most recent node.
	 */
	void setNodeHeight(double value);

	/**
	 * Set the height of this node relative to the most recent node.
	 * @param adjustChildBranchLengths if true
	 */
	void setNodeHeight(double value,boolean adjustChildBranchLengths);

	///** Set the height SE of this node relative to the most recent node. */
	//void setNodeHeightSE(double value);

	///** Get the height SE of this node relative to the most recent node. */
	//double getNodeHeightSE();


	/** Returns the identifier for this node. */
	Identifier getIdentifier();

	/** Set identifier for this node. */
	void setIdentifier(Identifier id);

	/**
	 * Returns the number of children this node has.
	 */
	int getChildCount();

	/**
	 * check whether this node is an external node
	 *
	 * @return result (true or false)
	 */
	boolean isLeaf();

	/**
	 * check whether this node is a root node
	 *
	 * @return result (true or false)
	 */
	boolean isRoot();

	/**
	 * get child node
	 *
	 * @param n number of child
	 *
	 * @return child node
	 */
	Node getChild(int n);

	/**
	 * set child node
	 *
	 * @param n number
	 * @node node new child node
	 */
	void setChild(int n, Node node);

	/**
	 * add new child node
	 *
	 * @param c new child node
	 */
	void addChild(Node c);

	/**
	 * add new child node (insertion at a specific position)
	 *
	 * @param c new child node
	 + @param pos position
	 */
	void insertChild(Node c, int pos);


	/**
	 * remove child
	 *
	 * @param n number of child to be removed
	 */
	Node removeChild(int n);
}
