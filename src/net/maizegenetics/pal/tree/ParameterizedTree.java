// ParameterizedTree.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package net.maizegenetics.pal.tree;

import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.math.OrthogonalHints;
import net.maizegenetics.pal.util.Parameterized;

import java.io.PrintWriter;
import java.io.StringWriter;


/**
 * abstract base class for a tree with an Parameterized interface
 *
 * @version $Id: ParameterizedTree.java,v 1.1 2007/01/12 03:26:17 tcasstevens Exp $
 *
 * @author Alexei Drummond
 * @author Korbinian Strimmer
 * @author Matthew Goode
 */
public interface ParameterizedTree extends Parameterized, Tree {
	OrthogonalHints getOrthogonalHints();
	String getParameterizationInfo();

	/**
	 * Factory interface
	 */
	static interface Factory {
		/**
		 * Generate a new parameterized tree wrapped around base
		 */
		ParameterizedTree generateNewTree(Tree base);
	}

	/**
	 * For parameterisations that work by adjusting a base tree (that is, they aren't really
	 * tree's themselves...)
	 * @note it should implment ParameterizedTree but, it that causes funny problems with my compiler ... I don't know why MG)
	 */
	static abstract class ParameterizedTreeBase implements  Parameterized, Tree{
		/**
		 * The non-parameterized tree that this parameterized tree is
		 * based on.
		 */
		private Tree tree;
		//
		// Public stuff
		//
		/**
		 * Cloning constructor
		 */
		protected ParameterizedTreeBase(ParameterizedTreeBase toCopy) {
			this.tree = toCopy.tree.getCopy();
		}

		public ParameterizedTreeBase() {}
		public ParameterizedTreeBase(Tree baseTree) {
			setBaseTree(baseTree);
		}

		protected void setBaseTree(Tree baseTree) {
			this.tree = baseTree;
			// make consistent
			createNodeList();
		}

		protected Tree getBaseTree() {		return tree;	}

		// interface tree

		/**
		 * Returns the root node of this tree.
		 */
		public final Node getRoot() { return tree.getRoot();	}

		public final void setRoot(Node root) {	tree.setRoot(root);	}

		//IdGroup stuff ==========

		public final Identifier getIdentifier(int i) {   return tree.getIdentifier(i);	  }
		public final void setIdentifier(int i,Identifier id) {   tree.setIdentifier(i,id);	  }
		public final int getIdCount() {		return tree.getIdCount();		}
		public final int whichIdNumber(String s) {		return tree.whichIdNumber(s);		}

	// ===========================================================================================

		/**
		 * returns a count of the number of external nodes (tips) in this
		 * tree.
		 */
		public final int getExternalNodeCount() {	return tree.getExternalNodeCount();		}

		/**
		 * returns a count of the number of internal nodes (and hence clades)
		 * in this tree.
		 */
		public final int getInternalNodeCount() {	return tree.getInternalNodeCount();}
		public final int getNodeCount() {	return tree.getInternalNodeCount()+tree.getExternalNodeCount();	}

		/**
		 * returns the ith external node in the tree.
		 */
		public final Node getExternalNode(int i) {	return tree.getExternalNode(i);	}

		/**
		 * returns the ith internal node in the tree.
		 */
		public final Node getInternalNode(int i) { return tree.getInternalNode(i);	}

		/**
		 * This method is called to ensure that the calls to other methods
		 * in this interface are valid.
		 */
		public final void createNodeList() { tree.createNodeList();	}
		public final int getUnits() {	return tree.getUnits();	}

		public final void setAttribute(Node node, String name, Object value) {
			tree.setAttribute(node, name, value);
		}
		public final Object getAttribute(Node node, String name) {
			return tree.getAttribute(node, name);
		}

		public String toString() {
			StringWriter sw = new StringWriter();
			NodeUtils.printNH(new PrintWriter(sw), getRoot(), true, false, 0, false);
			sw.write(";");
			return sw.toString();
		}

		/**
		 * The cheapy copy that just creates a SimpleTree
		 */
		public Tree getCopy() {
			return new SimpleTree((Tree)this);
		}
		// interface parameterized (remains abstract)

		/**
		 * @return null by default (implying not hint information)
		 */
		public OrthogonalHints getOrthogonalHints() {
			return null;
		}
	} //End of class Abstract
}
