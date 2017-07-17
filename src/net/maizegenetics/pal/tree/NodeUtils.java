// NodeUtils.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package net.maizegenetics.pal.tree;

import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.io.FormattedOutput;
import net.maizegenetics.pal.util.BranchLimits;
import net.maizegenetics.pal.util.Log;

import java.io.PrintWriter;
import java.util.Vector;

/**
 * Helper routines for dealing with nodes.
 *
 * @version $Id: NodeUtils.java,v 1.1 2007/01/12 03:26:17 tcasstevens Exp $
 *
 * @author Alexei Drummond
 * @author Korbinian Strimmer
 * @author Matthew Goode
 */
public class NodeUtils {
	/**
	 * Appends all external nodes from tree defined by root to Vector store
	 * @param root The root node defining tree
	 * @param store Where leaf nodes are stored (original contents is not touched)
	 */
	public static void getExternalNodes(Node root, Vector store) {
		if(root.isLeaf()) {
			store.addElement(root);
		} else {
			for(int i = 0 ; i < root.getChildCount() ; i++) {
				getExternalNodes(root.getChild(i),store);
			}
		}
	}
	/**
	 * Obtains all external nodes from tree defined by root and returns as an array
	 * @param root The root node defining tree
	 * @return an array of nodes where each node is a leaf node, and is a member
	 * of the tree defined by root
	 */
	public static Node[] getExternalNodes(Node root) {
		Vector v = new Vector();
		getExternalNodes(root,v);
		Node[] result = new Node[v.size()];
		v.copyInto(result);
		return result;
	}
	/**
	 * Appends all internal nodes from tree defined by root to Vector store
	 * @param root The root node defining tree
	 * @param store Where internal nodes are stored (original contents is not touched)
	 * @note Root will be the first node added
	 */
	public static void getInternalNodes(Node root, Vector store) {
		if(!root.isLeaf()) {
			store.addElement(root);
			for(int i = 0 ; i < root.getChildCount() ; i++) {
				getInternalNodes(root.getChild(i),store);
			}
		}
	}
	/**
	 * Obtains all internal nodes from tree defined by root and returns as an array
	 * @param root The root node defining tree
	 * @return an array of nodes where each node is a internal node, and is a member
	 * of the tree defined by root
	 * @note Root will be the first node added (if included)
	 */
	public static Node[] getInternalNodes(Node root, boolean includeRoot) {
		Vector v = new Vector();
		getInternalNodes(root,v);
		Node[] result = new Node[v.size()];
		v.copyInto(result);
		if(includeRoot) {
			return result;
		}
		Node[] adjustedResult = new Node[result.length-1];
		System.arraycopy(result, 1,adjustedResult,0,adjustedResult.length);
		return adjustedResult;
	}

	/**
	 * @return the maxium distance in nodes from root to a leaf
	 * @note this is a CompSci depth measure and has nothing to do with branch lengths/node heights :)
	 */
	public static int getMaxNodeDepth(Node root) {
		int max = 0;
		for(int i = 0 ; i < root.getChildCount() ; i++) {
			int depth = getMaxNodeDepth(root.getChild(i));
			if(depth>max) {
				max = depth;
			}
		}
		return max+1;
	}

	/**
	 * Calculates max/min lengths of paths from root to leaf, taking into account branch lengths
	 *
	 * @return a double array where
	 *  <ul>
	 *    <li>The first element is the minimum length path from root to leaf.</li>
	 *    <li>The second element is the second most minimum length path from root to leaf.</li>
	 *    <li>The third element is the second most maximum length path from root to leaf.</li>
	 *    <li>The forth element is the maximum length path from root to leaf.</li>
	 *  </ul>
	 * @see getMaxNodeDepth()
	 */
	public static final double[] getPathLengthInfo(Node root) {
		double[] lengthInfo = new double[] { Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY };
		getLengthInfo(root,0,lengthInfo);
		return lengthInfo;
	}
	/**
	 * @return the maximum length of a path from root to a leaf
	 */
	public static final double getMaximumPathLengthLengthToLeaf(Node root) {
		if(root.isLeaf()) { return 0; }
		double maxLength = Double.NEGATIVE_INFINITY;
		for(int i = 0 ; i < root.getChildCount() ; i++) {
			Node c = root.getChild(i);
			double length = c.getBranchLength()+getMaximumPathLengthLengthToLeaf(c);
			maxLength=Math.max(length,maxLength);
		}
		return maxLength;
	}
	/**
	 * @return the minimum length of a path from root to a leaf
	 */
	public static final double getMinimumPathLengthLengthToLeaf(Node root) {
		if(root.isLeaf()) { return 0; }
		double minLength = Double.POSITIVE_INFINITY;
		for(int i = 0 ; i < root.getChildCount() ; i++) {
			Node c = root.getChild(i);
			double length = c.getBranchLength()+getMinimumPathLengthLengthToLeaf(c);
			minLength=Math.min(length,minLength);
		}
		return minLength;
	}

	/**
	 * Traverses tree (recursively) updating lengthInfo each time a leaf node is met.
	 * @note: first value is assumed to be minimum length, second value is assumed to be maximum length
	 */
	private static final void getLengthInfo(Node root, double lengthFromRoot, double[] lengthInfo) {
		if(root.isLeaf()) {
			if(lengthFromRoot<lengthInfo[1]) {
				if(lengthFromRoot<lengthInfo[0]) {
					lengthInfo[1] = lengthInfo[0];
					lengthInfo[0] = lengthFromRoot;
				} else {
					lengthInfo[1] = lengthFromRoot;
				}
			}
			if(lengthFromRoot>lengthInfo[2]) {
				if(lengthFromRoot>lengthInfo[3]) {
					lengthInfo[2] = lengthInfo[3];
					lengthInfo[3] = lengthFromRoot;
				} else {
					lengthInfo[2] = lengthFromRoot;
				}
			}
		} else {
			for(int i = 0 ; i < root.getChildCount() ; i++) {
				Node c = root.getChild(i);
				getLengthInfo(c, lengthFromRoot+c.getBranchLength(),lengthInfo);
			}
		}
	}

	/**
	 * Converts lengths to heights, *without* assuming contemporaneous
	 * tips.
	 */
	public static void lengths2Heights(Node root) {

		lengths2Heights(root, getMaximumPathLengthLengthToLeaf(root));
	}

	/**
	 * @return the number of internal nodes from a root node
	 */
	public static int getInternalNodeCount(Node root) {
		if(root.isLeaf()) {
			return 0;
		}
		int count = 0;
		for(int i = 0 ; i < root.getChildCount() ; i++) {
			count+=getInternalNodeCount(root.getChild(i));
		}
		return count+1;
	}

	/**
	 * Converts lengths to heights, but maintains tip heights.
	 */
	public static void lengths2HeightsKeepTips(Node node, boolean useMax) {

		if (!node.isLeaf()) {
			for (int i = 0; i < node.getChildCount(); i++) {
				lengths2HeightsKeepTips(node.getChild(i), useMax);
			}

			double totalHL = 0.0;
			double maxHL = 0.0;
			double hl = 0.0;
			double maxH = 0.0;
			double h = 0.0;
			for (int i = 0; i < node.getChildCount(); i++) {
				h = node.getChild(i).getNodeHeight();
				hl = node.getChild(i).getBranchLength() + h;
				if (hl > maxHL) maxHL = hl;
				if (h > maxH) maxH = h;
				totalHL += hl;
			}
			if (useMax) {
				hl = maxHL; // set parent height to maximum parent height implied by children
			} else {
				hl = totalHL /  node.getChildCount(); // get mean parent height
				if (hl < maxH) hl = maxHL; // if mean parent height is not greater than all children height, fall back on max parent height.
			}
			node.setNodeHeight(hl); // set new parent height

			// change lengths in children to reflect changes.
			for (int i = 0; i < node.getChildCount(); i++) {
				h = node.getChild(i).getNodeHeight();
				node.getChild(i).setBranchLength(hl - h);
			}
		}
	}


	/**
	 * sets this nodes height value to newHeight and all children's
	 * height values based on length of branches.
	 */
	private static void lengths2Heights(Node node, double newHeight) {

		if (!node.isRoot()) {
			newHeight -= node.getBranchLength();
			node.setNodeHeight(newHeight);
		} else {
			node.setNodeHeight(newHeight);
		}

		for (int i = 0; i < node.getChildCount(); i++) {
			lengths2Heights(node.getChild(i), newHeight);
		}
	}

	/**
	 * Exchange field info between two nodes. Specifically
	 * identifiers, branch lengths, node heights and branch length
	 * SEs.
	 */
	public static void exchangeInfo(Node node1, Node node2) {

		Identifier swaps;
		double swapd;

		swaps = node1.getIdentifier();
		node1.setIdentifier(node2.getIdentifier());
		node2.setIdentifier(swaps);

		swapd = node1.getBranchLength();
		node1.setBranchLength(node2.getBranchLength());
		node2.setBranchLength(swapd);

		swapd = node1.getNodeHeight();
		node1.setNodeHeight(node2.getNodeHeight());
		node2.setNodeHeight(swapd);

		swapd = node1.getBranchLengthSE();
		node1.setBranchLengthSE(node2.getBranchLengthSE());
		node2.setBranchLengthSE(swapd);
	}

	/**
	 * determines branch lengths of this and all descendent nodes
	 * from heights
	 */
	public static void heights2Lengths(Node node) {
		heights2Lengths(node, true); //respect minimum
	}

	/**
	 * determines branch lengths of this and all descendent nodes
	 * from heights
	 */
	public static void heights2Lengths(Node node, boolean respectMinimum) {
		for (int i = 0; i < node.getChildCount(); i++) {
			heights2Lengths(node.getChild(i));
		}

		if (node.isRoot()) {
			node.setBranchLength(0.0);
		}
		else {
			node.setBranchLength(node.getParent().getNodeHeight() - node.getNodeHeight());
			if (respectMinimum && (node.getBranchLength() < BranchLimits.MINARC))
			{
				node.setBranchLength(BranchLimits.MINARC);
			}
		}
	}

	/**
	 * determines branch lengths of this node and its immediate descendent nodes
	 * from heights.
	 */
	public static void localHeights2Lengths(Node node, boolean respectMinimum) {

		for (int i = 0; i < node.getChildCount(); i++) {
			Node child = node.getChild(i);

			child.setBranchLength(node.getNodeHeight() - child.getNodeHeight());
		}

		if (node.isRoot()) {
			node.setBranchLength(0.0);
		}
		else {
			node.setBranchLength(node.getParent().getNodeHeight() - node.getNodeHeight());
			if (respectMinimum && (node.getBranchLength() < BranchLimits.MINARC))
			{
				node.setBranchLength(BranchLimits.MINARC);
			}
		}
	}


	/**
	 * Finds the largest child (in terms of node height).
	 */
	public static double findLargestChild(Node node) {
		// find child with largest height
		double max = node.getChild(0).getNodeHeight();
		for (int j = 1; j < node.getChildCount(); j++){
			double h = node.getChild(j).getNodeHeight();
			if (h > max)
			{
				max = h;
			}
		}
		return max;
	}

	/**
	 * remove child
	 *
	 * @param node child node to be removed
	 */
	public static void removeChild(Node parent, Node child)
	{
		int rm = -1;
		for (int i = 0; i < parent.getChildCount(); i++)
		{
			if (child == parent.getChild(i))
			{
				rm = i;
				break;
			}
		}

		parent.removeChild(rm);
	}

	/**
	 * remove internal branch (collapse node with its parent)
	 *
	 * @param node node associated with internal branch
	 */
	public static void removeBranch(Node node)
	{
		if (node.isRoot() || node.isLeaf())
		{
			throw new IllegalArgumentException("INTERNAL NODE REQUIRED (NOT ROOT)");
		}

		Node parent = node.getParent();

		// add childs of node to parent
		// (node still contains the link to childs
		// to allow later restoration)
		int numChilds = node.getChildCount();
		for (int i = 0; i < numChilds; i++)
		{
			parent.addChild(node.getChild(i));
		}

		// remove node from parent
		// (link to parent is restored and the
		// position is stored)
		int rm = -1;
		for (int i = 0; i < parent.getChildCount(); i++)
		{
			if (node == parent.getChild(i))
			{
				rm = i;
				break;
			}
		}
		parent.removeChild(rm);
		node.setParent(parent);
		node.setNumber(rm);
	}

	/**
	 * restore internal branch
	 *
	 * @param node node associated with internal branch
	 */
	public static void restoreBranch(Node node)
	{
		if (node.isRoot() || node.isLeaf())
		{
			throw new IllegalArgumentException("INTERNAL NODE REQUIRED (NOT ROOT)");
		}

		Node parent = node.getParent();

		// remove childs of node from parent and make node their parent
		int numChilds = node.getChildCount();
		for (int i = 0; i < numChilds; i++)
		{
			Node c = node.getChild(i);
			removeChild(parent, c);
			c.setParent(node);
		}

		// insert node into parent
		parent.insertChild(node, node.getNumber());
	}



	/**
	 * join two childs, introducing a new node/branch in the tree
	 * that replaces the first child
	 *
	 * @param n1 number of first child
	 * @param n2 number of second child
	 */
	public static void joinChilds(Node node, int n1, int n2) {

		if (n1 == n2) {
			throw new IllegalArgumentException("CHILDREN MUST BE DIFFERENT");
		}

		int c1, c2;
		if (n2 < n1)
		{
			c1 = n2;
			c2 = n1;
		}
		else
		{
			c1 = n1;
			c2 = n2;
		}

		Node newNode = NodeFactory.createNode();

		Node child1 = node.getChild(c1);
		Node child2 = node.getChild(c2);

		node.setChild(c1, newNode);
		newNode.setParent(node);
		node.removeChild(c2); // now parent of child2 = null

		newNode.addChild(child1);
		newNode.addChild(child2);
	}

	/**
	 * determine preorder successor of this node
	 *
	 * @return next node
	 */
	public static Node preorderSuccessor(Node node) {

		Node next = null;

		if (node.isLeaf()) {
			Node cn = node, ln = null; // Current and last node

			// Go up
			do
			{
				if (cn.isRoot())
				{
					next = cn;
					break;
				}
				ln = cn;
				cn = cn.getParent();
			}
			while (cn.getChild(cn.getChildCount()-1) == ln);

			// Determine next node
			if (next == null)
			{
				// Go down one node
				for (int i = 0; i < cn.getChildCount()-1; i++)
				{
					if (cn.getChild(i) == ln)
					{
						next = cn.getChild(i+1);
						break;
					}
				}
			}
		}
		else
		{
			next = node.getChild(0);
		}

		return next;
	}

	/**
	 * determine postorder successor of a node
	 *
	 * @return next node
	 */
	public static Node postorderSuccessor(Node node) {

		Node cn = null;
		Node parent = node.getParent();

		if (node.isRoot()){
			cn = node;
		}	else{

			// Go up one node
			if (parent.getChild(parent.getChildCount()-1) == node) {
				return parent;
			}
			// Go down one node
			for (int i = 0; i < parent.getChildCount()-1; i++)	{
				if (parent.getChild(i) == node)	{
					cn = parent.getChild(i+1);
					break;
				}
			}
		}
		// Go down until leaf
		while (cn.getChildCount() > 0)
		{

			cn = cn.getChild(0);
		}


		return cn;
	}

	/**
	 * prints node in New Hamshire format.
	 */
	public static void printNH(PrintWriter out, Node node,
		boolean printLengths, boolean printInternalLabels) {

		printNH(out, node, printLengths, printInternalLabels, 0, true);
	}


	public static int printNH(PrintWriter out, Node node,
		boolean printLengths, boolean printInternalLabels, int column, boolean breakLines) {

		if (breakLines) column = breakLine(out, column);

		if (!node.isLeaf())
		{
			out.print("(");
			column++;

			for (int i = 0; i < node.getChildCount(); i++)
			{
				if (i != 0)
				{
					out.print(",");
					column++;
				}

				column = printNH(out, node.getChild(i), printLengths, printInternalLabels, column, breakLines);
			}

			out.print(")");
			column++;
		}

		if (!node.isRoot())
		{
			if (node.isLeaf() || printInternalLabels)
			{
				if (breakLines) column = breakLine(out, column);

				String id = node.getIdentifier().toString();
				out.print(id);
				column += id.length();
			}

			if (printLengths)
			{
				out.print(":");
				column++;

				if (breakLines) column = breakLine(out, column);

				column += FormattedOutput.getInstance().displayDecimal(out, node.getBranchLength(), 7);
			}
		}

		return column;
	}

	private static int breakLine(PrintWriter out, int column)
	{
		if (column > 70)
		{
			out.println();
			column = 0;
		}

		return column;
	}
	/**
	 * Returns the first nodes in this tree that has the
	 * required identifiers.
	 * @return null if none of the identifiers names match nodes in tree, else return array which may have
	 *  null "blanks" for corresponding identifiers that do not match any node in the tree
	 */
	public static final Node[] findByIdentifier(Node node, String[] identifierNames) {
		Node[] nodes = new Node[identifierNames.length];
		boolean foundSomething = false;
		for(int i = 0 ; i < nodes.length ; i++) {
			nodes[i] = findByIdentifier(node,identifierNames[i]);
			foundSomething = foundSomething||(nodes[i]!=null);
		}
		if(!foundSomething) {
			return null;
		}
		return nodes;
	}
	/**
	 * Returns the first nodes in this tree that has the
	 * required identifiers.
	 */
	public static final Node[] findByIdentifier(Node node, Identifier[] identifiers) {
		Node[] nodes = new Node[identifiers.length];
		for(int i = 0 ; i < nodes.length ; i++) {
			nodes[i] = findByIdentifier(node,identifiers[i]);
		}
		return nodes;
	}
	/**
	 * Returns the first node in this tree that has the
	 * required identifier.
	 */
	public static final Node findByIdentifier(Node node, Identifier identifier) {
		return findByIdentifier(node,identifier.getName());
	}
	/**
	 * Returns the first node in this tree that has the
	 * required identifier.
	 */
	public static final Node findByIdentifier(Node node, String identifierName) {

		Log.getDefaultLogger().debug("node identifier = " + node.getIdentifier());
		Log.getDefaultLogger().debug("target identifier name = " + identifierName);

		if (node.getIdentifier().getName().equals(identifierName)) {
			return node;
		} else {
			Node pos = null;
			for (int i = 0; i < node.getChildCount(); i++) {
				pos = findByIdentifier(node.getChild(i), identifierName);
				if (pos != null) return pos;
			}
			//if (pos == null && !node.isRoot()) {
			//	pos = findByIdentifier(node.getParent(), identifier);
			//}
			if (pos != null) return pos;
			return null;
		}
	}



	/**
	 * determine distance to root
	 *
	 * @return distance to root
	 */
	public static double getDistanceToRoot(Node node)
	{
		if (node.isRoot())
		{
			return 0.0;
		}
		else
		{
			return node.getBranchLength() + getDistanceToRoot(node.getParent());
		}
	}

	/**
	 * Return the number of terminal leaves below this node or 1 if this is
	 * a terminal leaf.
	 */
	public static int getLeafCount(Node node) {

		int count = 0;
		if (!node.isLeaf()) {
			for (int i = 0; i < node.getChildCount(); i++) {
				count += getLeafCount(node.getChild(i));
			}
		} else {
			count = 1;
		}
		return count;
	}
	/**
	 * For two nodes in the tree true if the first node is the ancestor of the second
	 *
	 * @param possibleAncestor the node that may be the ancestor of the other node
	 * @param node the node that may have the other node as it's ancestor
	 */
	public static boolean isAncestor(Node possibleAncestor, Node node) {
		if(node==possibleAncestor) {
			return true;
		}
		while(!node.isRoot()){
			node = node.getParent();
			if(node==possibleAncestor) {
				return true;
			}
		}
		return false;
	}

	/**
	 * For a set of nodes in the tree returns the common ancestor closest to all nodes (most recent common ancestor)
	 *
	 * @param nodes the nodes to check, is okay if array elements are null!
	 * @returns null if a at least one node is disjoint from the others nodes disjoint
	 */
	public static Node getFirstCommonAncestor(Node[] nodes) {
		Node currentCA = nodes[0];
		for(int i = 1; i < nodes.length ;i++) {
			if(currentCA!=null&&nodes[i]!=null) {
				currentCA = getFirstCommonAncestor(currentCA,nodes[i]);
				if(currentCA==null) {
					return null;
				}
			}
		}
		return currentCA;
	}
	/**
	 * For two nodes in the tree returns the common ancestor closest to both nodes (most recent common ancestor)
	 *
	 * @param nodeOne
	 * @param nodeTwo
	 * @returns null if two nodes disjoint (from different trees). May also return either nodeOne or nodeTwo if one node is an ancestor of the other
	 */
	public static Node getFirstCommonAncestor(Node nodeOne, Node nodeTwo) {
		if(isAncestor(nodeTwo, nodeOne)) {
			return nodeTwo;
		}
		if(isAncestor(nodeOne, nodeTwo)) {
			return nodeOne;
		}
		while(!nodeTwo.isRoot()) {
			nodeTwo = nodeTwo.getParent();
			if(isAncestor(nodeTwo, nodeOne)) {
				return nodeTwo;
			}
		}
		return null;
	}

		/** returns number of branches centered around an internal node in an unrooted tree */
	public static final int getUnrootedBranchCount(Node center) {
		if (center.isRoot()) 	{
			return center.getChildCount();
		}
		else {
			return center.getChildCount()+1;
		}
	}


	/**
	 * If the given node or the sub tree defined by that node have negative branch lengths, they'll have
	 * zeron branch lengths after a call to this function!
	 */
	public static final void convertNegativeBranchLengthsToZeroLength(Node node) {
		convertNegativeBranchLengthsToZeroLengthImpl(node);
		lengths2Heights(node);
	}
	private static final void convertNegativeBranchLengthsToZeroLengthImpl(Node node) {
		if(node.getBranchLength()<0) {
			node.setBranchLength(0);
		}
		for(int i = 0 ; i < node.getChildCount() ; i++) {
			convertNegativeBranchLengthsToZeroLengthImpl(node.getChild(i));
		}
	}
}
