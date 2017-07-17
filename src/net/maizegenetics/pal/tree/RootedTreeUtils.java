// RootedTreeUtils.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.tree;

import java.util.Hashtable;
import java.util.Vector;

/**
 * This class contains utility methods. These include: <BR>
 * 1. gathering information about subtrees from a set of trees <BR>
 * 2. comparing subtrees and clades. <BR>
 * All these methods assume rooted trees!
 *
 * @author Alexei Drummond
 * @version  $Id: RootedTreeUtils.java,v 1.1 2007/01/12 03:26:17 tcasstevens Exp $
 */
public class RootedTreeUtils {

	/**
	 * @return true if the first node contains a subtree identical to the second node
	 * or is identical to the second node.
	 * @param root the root of the tree in which search for a subtree
	 * @param node the subtree to search for.
	 */
	public static boolean containsSubtree(Node root, Node node) {
		return (getSubtree(root, node) != null);
	}

	/**
	 * @return true if the given tree contains a clade holding all the
	 * taxa in the given subtree.
	 *
	 * @param root the root of the tree in which search for a subtree
	 * @param taxa the hashtable of taxa.
	 */
	public static boolean containsClade(Node root, Node clade) {
		return (getClade(root, clade) != null);
	}

	/**
	 * @return a subtree within the first node with the same
	 * labelled topology as the second node or null if it doesn't exist.
	 */
	public static Node getSubtree(Node root, Node node) {
		if (equal(root, node)) return root;
		for (int i =0; i < root.getChildCount(); i++) {
			Node match = getSubtree(root.getChild(i), node);
			if (match != null) return match;
		}
		return null;
	}

	/**
	 * @return a subtree within the first node with the same
	 * labels as the second node or null if it doesn't exist.
	 */
	public static Node getClade(Node root, Node clade) {
		if (sameTaxa(root, clade)) return root;
		for (int i =0; i < root.getChildCount(); i++) {
			Node match = getClade(root.getChild(i), clade);
			if (match != null) return match;
		}
		return null;
	}



	/**
	 * @return true if the trees have the same tip-labelled
	 * structure. Child order is not important.
	 */
	public static boolean equal(Node node1, Node node2) {
		int nodeCount1 = node1.getChildCount();
		int nodeCount2 = node2.getChildCount();

		// if different childCount not the same
		if (nodeCount1 != nodeCount2) return false;

		if (nodeCount1 == 0) {
			return (node1.getIdentifier().getName().equals(node2.getIdentifier().getName()));
		} else {
			// ASSUMES BIFURCATING TREES
			// CHILD ORDER DIFFERENCES ARE ALLOWED!
			if (equal(node1.getChild(0), node2.getChild(0))) {
				return (equal(node1.getChild(1), node2.getChild(1)));
			} else if (equal(node1.getChild(0), node2.getChild(1))) {
				return (equal(node1.getChild(1), node2.getChild(0)));
			} else return false;
		}
	}

	/**
	 * @return true if the trees have the same tip labels. topology unimportant.
	 */
	public static boolean sameTaxa(Node node1, Node node2) {
		int leafCount1 = NodeUtils.getLeafCount(node1);
		int leafCount2 = NodeUtils.getLeafCount(node2);

		if (leafCount1 != leafCount2) return false;

		Hashtable table = new Hashtable(leafCount1+1);
		collectTaxa(node1, table);
		return !containsNovelTaxa(node2, table);
	}

	/**
	 * Collects all of the names of the taxa in the tree into a hashtable.
	 * @return the number of new taxa added to the hashtable from this tree.
	 * @param root the root node of the tree.
	 * @param taxa a hashtable to hold the taxa names, may already hold some taxa names.
	 */
	public static int collectTaxa(Node root, Hashtable table) {
		int nc = root.getChildCount();
		if (nc == 0) {
			String name = root.getIdentifier().getName();
			if (table.containsKey(name)) {
				return 0;
			} else {
				table.put(name, name);
				return 1;
			}
		} else {
			int newTaxaCount = 0;
			for (int i = 0; i < nc; i++) {
				newTaxaCount += collectTaxa(root.getChild(i), table);
			}
			return newTaxaCount;
		}
	}

	/**
	 * @return true if the given tree contains taxa not already in the given hashtable.
	 * @param root the root node of the tree.
	 * @param taxa a hashtable holding taxa names.
	 */
	public static boolean containsNovelTaxa(Node root, Hashtable taxa) {
		int nc = root.getChildCount();
		if (nc == 0) {
			return !taxa.containsKey(root.getIdentifier().getName());
		} else {
			for (int i = 0; i < nc; i++) {
				if (containsNovelTaxa(root.getChild(i), taxa)) return true;
			}
			return false;
		}
	}



	/**
	 * @return the number of taxa in the given tree that are NOT in the given hashtable.
	 * @param root the root node of the tree.
	 * @param taxa a hashtable holding taxa names.
	 */
	private static int newTaxaCount(Node root, Hashtable table) {
		int nc = root.getChildCount();
		if (nc == 0) {
			return (table.containsKey(root.getIdentifier().getName()) ? 0 : 1);
		} else {
			int newTaxaCount = 0;
			for (int i = 0; i < nc; i++) {
				newTaxaCount += newTaxaCount(root.getChild(i), table);
			}
			return newTaxaCount;
		}
	}

	/**
	 * @return the number of times the subtree was found in the
	 * given list of trees. If a subtree occurs more than once in a tree
	 * (for some bizarre reason) it is counted only once.
	 * @param subtree the subtree being searched for.
	 * @param trees a vector of trees to search for the subtree in.
	 */
	public static int subtreeCount(Node subtree, Vector trees) {
		int count = 0;
		Node root;
		for (int i = 0; i < trees.size(); i++) {
			root = ((Tree)trees.elementAt(i)).getRoot();
			if (containsSubtree(root, subtree)) {
				count += 1;
			}
		}
		return count;
	}

	/**
	 * @return the mean height of the given subtree in the
	 * given list of trees. If a subtree occurs more than once in a tree
	 * (for some bizarre reason) results are undefined.
	 * @param subtree the subtree being searched for.
	 * @param trees a vector of trees to search for the subtree in.
	 */
	public static double getMeanSubtreeHeight(Node subtree, Vector trees) {
		int count = 0;
		double totalHeight = 0.0;
		Node root;
		for (int i = 0; i < trees.size(); i++) {
			root = ((Tree)trees.elementAt(i)).getRoot();
			Node match = getSubtree(root, subtree);
			if (match != null) {
				count += 1;
				totalHeight += match.getNodeHeight();
			}
		}
		return totalHeight / (double)count;
	}

	/**
	 * @return the mean height of the given clade in the
	 * given list of trees. If a clade occurs more than once in a tree
	 * (for some bizarre reason) results are undefined.
	 * @param clade a node containing the clade being searched for.
	 * @param trees a vector of trees to search for the clade in.
	 */
	public static double getMeanCladeHeight(Node clade, Vector trees) {
		int count = 0;
		double totalHeight = 0.0;
		Node root;
		for (int i = 0; i < trees.size(); i++) {
			root = ((Tree)trees.elementAt(i)).getRoot();
			Node match = getClade(root, clade);
			if (match != null) {
				count += 1;
				totalHeight += match.getNodeHeight();
			}
		}
		return totalHeight / (double)count;
	}

	/**
	 * @return the number of times the clade was found in the
	 * given list of trees. If a clade occurs more than once in a tree
	 * (for some bizarre reason) it is counted only once.
	 * @param subtree a subtree containing the taxaset being searched for.
	 * @param trees a vector of trees to search for the clade in.
	 */
	public static int cladeCount(Node subtree, Vector trees) {
		int count = 0;
		Node root;
		for (int i = 0; i < trees.size(); i++) {
			root = ((Tree)trees.elementAt(i)).getRoot();
			if (containsClade(root, subtree)) {
				count += 1;
			}
		}
		return count;
	}

	public static void collectProportions(Tree tree, Vector trees) {

		for (int i =0; i < tree.getInternalNodeCount(); i++) {
			Node node = tree.getInternalNode(i);
			if (!node.isRoot()) {
				int cladeCount = cladeCount(node, trees);
				StringBuffer buffer = new StringBuffer();
				collectLeafNames(node, buffer);
				double pr = (double)cladeCount / (double)trees.size();
				tree.setAttribute(node, AttributeNode.CLADE_PROBABILITY, new Double(pr));

				double meanCladeHeight = getMeanCladeHeight(node, trees);
				tree.setAttribute(node, AttributeNode.MEAN_CLADE_HEIGHT, new Double(meanCladeHeight));

			}
			int subtreeCount = subtreeCount(node, trees);
			double pr = (double)subtreeCount / (double)trees.size();
			tree.setAttribute(node, AttributeNode.SUBTREE_PROBABILITY, new Double(pr));
		}
	}

	/**
	 * Fills given string buffer with preorder traversal space-delimited leaf names.
	 */
	private static void collectLeafNames(Node node, StringBuffer buffer) {

		if (node.isLeaf()) {
			buffer.append(node.getIdentifier().getName());
			buffer.append(' ');
		} else {
			for (int i = 0; i < node.getChildCount(); i++) {
				collectLeafNames(node.getChild(i), buffer);
			}
		}
	}

}

