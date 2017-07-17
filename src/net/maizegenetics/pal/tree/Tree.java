// Tree.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package net.maizegenetics.pal.tree;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.report.Report;
import net.maizegenetics.pal.ids.IdGroupUtils;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.io.StringWriter;
import java.util.Hashtable;



/**
 * Interface for a phylogenetic or genealogical tree.
 *
 * @version $Id: Tree.java,v 1.2 2009/06/24 06:20:05 tcasstevens Exp $
 *
 * @author Alexei Drummond
 */
public interface Tree extends IdGroup, Units, Serializable, UnitsProvider {

	/**
	 * @return the root node of this tree.
	 */
	Node getRoot();

	/**
	 * This method constructs a tree from the given root node.
	 * @param root the root node of the tree to construct.
	 */
	void setRoot(Node root);

	/**
	 * @return a count of the number of external nodes (tips) in this
	 * tree.
	 */
	int getExternalNodeCount();

	/**
	 * @return a count of the number of internal nodes (and hence clades)
	 * in this tree.
	 */
	int getInternalNodeCount();

	/**
	 * @return the ith external node in the tree.
	 */
	Node getExternalNode(int i);

	/**
	 * @return the ith internal node in the tree.
	 */
	Node getInternalNode(int i);

	/**
	 * This method is called to ensure that the calls to other methods
	 * in this interface are valid.
	 */
	void createNodeList();

	/**
	 * Gets the units that this tree's branch lengths and node
	 * heights are expressed in.
	 */
	int getUnits();

	/**
	 * Sets the units that this tree's branch lengths and node
	 * heights are expressed in.
	 */
	//void setUnits(int units);

	/**
	 * Sets an named attribute for a given node.
	 * @param node the node whose attribute is being set.
	 * @param name the name of the attribute.
	 * @param value the new value of the attribute.
	 */
	void setAttribute(Node node, String name, Object value);

	/**
	 * @return an object representing the named attributed for the numbered node.
	 * @param node the node being interrogated.
	 * @param name the name of the attribute of interest.
	 */
	Object getAttribute(Node node, String name);


	/**
	 * @return a clone of this tree
	 */
	 public Tree getCopy();

//====== Abstract Implementation ============================
	public static abstract class TreeBase implements Tree, Report, Units, Serializable{
		//
		// This class has explicit serialization code so if you alter any fields please alter
		// the serialization code too (make sure you use a new version number - see readObject/writeObject
		// Thanks, Matthew

		//
		// Public stuff
		//

		//
		// Private stuff
		/** root node */
		private Node root;

		/** list of internal nodes (including root) */
		private Node[] internalNode = null;

		/** number of internal nodes (including root) */
		private int numInternalNodes;

		/** list of external nodes */
		private Node[] externalNode = null;

		/** number of external nodes */
		private int numExternalNodes;

		/** attributes attached to this tree. */
		private Hashtable[] attributes = null;

		/** holds the units of the trees branches. */
		private int units = EXPECTED_SUBSTITUTIONS;

		private boolean setupLengthsAndHeights_ = true;
		//
		// Serialization Stuff
		//

	//static final long serialVersionUID=-7330318631600898531L;

	//serialver -classpath ./classes net.maizegenetics.pal.tree.Tree.AbstractTree

		/** I like doing things my self! */
		private void writeObject(java.io.ObjectOutputStream out) throws java.io.IOException {
			out.writeByte(2); //Version number
			out.writeObject(root);
			out.writeObject(attributes);
			out.writeInt(units);
			out.writeBoolean(setupLengthsAndHeights_);
		}

		private void readObject(java.io.ObjectInputStream in) throws IOException, ClassNotFoundException{
			byte version = in.readByte();
			switch(version) {
				case 1 : {
					root = (Node)in.readObject();
					createNodeList();
					attributes = (Hashtable[])in.readObject();
					units = in.readInt();
				}
				default : {
					root = (Node)in.readObject();
					createNodeList();
					attributes = (Hashtable[])in.readObject();
					units = in.readInt();
					setupLengthsAndHeights_ = in.readBoolean();
				}
			}
		}
		/**
		 * constructor
		 * @param theRoot the node to base tree on
		 */
		protected TreeBase(Node theRoot) {
			setRoot(theRoot);
		}

		/**
		 * constructor
		 * @note Need to call setRoot() at some point!
		 */
		protected TreeBase() {	}

		/**
		 * clone constructor
		 */
		protected TreeBase(TreeBase tree, Node newRoot) {
			setRoot(newRoot);
			setUnits(tree.getUnits());
		}
		/**
		 * Called by subclasses that can't give root at time of construction (of this class)
		 */
		public void setRoot(Node theRoot) {
			this.root = theRoot;
			createNodeList();
		}
		/**
			* Return the units that this tree is expressed in.
			*/
		public final int getUnits() { return units;	}

		/**
		 * Sets the units that this tree is expressed in.
		 */
		public final void setUnits(int units) { this.units = units;	}

		/**
		 * Returns the number of external nodes.
		 */
		public final int getExternalNodeCount() {
			if(externalNode==null) {	createNodeList();		}
			return numExternalNodes;
		}

		/**
		 * Returns the ith external node.
		 */
		public final Node getExternalNode(int i) {
			if(externalNode==null) {	createNodeList();		}
			return externalNode[i];
		}

		/**
		 * Returns the number of internal nodes.
		 */
		public final int getInternalNodeCount() {
			if(internalNode==null) {	createNodeList();		}
			return numInternalNodes;
		}

		/**
		 * Returns the ith internal node.
		 */
		public final Node getInternalNode(int i) {
			if(internalNode==null) {	createNodeList();	}
			return internalNode[i];
		}

		/**
		 * Returns the root node of this tree.
		 */
		public final Node getRoot() {	return root; }

		/**
		 * return node with number num (as displayed in ASCII tree)
		 *
		 * @param num number of node
		 *
		 * @return node
		 */
		public Node findNode(int num)	{
			createNodeList();
			if (num <= numExternalNodes) {
				return externalNode[num-1];
			} else {
				return internalNode[num-1-numExternalNodes];
			}
		}

		/**
		 * count and list external and internal nodes and
		 * compute heights of each node
		 */
		public void createNodeList()	{
			numInternalNodes = 0;
			numExternalNodes = 0;
			Node node = root;
			do
			{
				node = NodeUtils.postorderSuccessor(node);
				if (node.isLeaf())
				{
					node.setNumber(numExternalNodes);
					numExternalNodes++;
				}
				else
				{
					node.setNumber(numInternalNodes);
					numInternalNodes++;
				}
			}
			while(node != root);

			internalNode = new Node[numInternalNodes];
			externalNode = new Node[numExternalNodes];
			node = root;
			do
			{
				node = NodeUtils.postorderSuccessor(node);
				if (node.isLeaf())
				{
					externalNode[node.getNumber()] = node;
				}
				else
				{
					internalNode[node.getNumber()] = node;
				}
			}
			while(node != root);

			// compute heights if it seems necessary
			if (setupLengthsAndHeights_&&root.getNodeHeight() == 0.0) {
				NodeUtils.lengths2Heights(root);
			}
		}
		private int getIndex(Node node) {
			if (node.isLeaf()) return node.getNumber();
			return getExternalNodeCount() + node.getNumber();
		}

		/**
		 * Set whether or not this class should call methods to adjust
		 * either heights or lengths automatically
		 * @note by default this is true
		 */
		protected void setSetupLengthsAndHeights(boolean value) {
			this.setupLengthsAndHeights_ = value;
		}
		public String toString() {
			StringWriter sw = new StringWriter();
			NodeUtils.printNH(new PrintWriter(sw), getRoot(), true, false, 0, false);
			sw.write(";");
			return sw.toString();
		}

		/**
		 * Sets an named attribute for a given node.
		 * @param node the node whose attribute is being set.
		 * @param name the name of the attribute.
		 * @param value the new value of the attribute.
		 */
		public void setAttribute(Node node, String name, Object value) {
			if (node instanceof AttributeNode) {
				((AttributeNode)node).setAttribute(name, value);
			} else {
				int index = getIndex(node);
				if (attributes == null) {
					attributes = new Hashtable[getExternalNodeCount() + getInternalNodeCount()];
				}
				if (attributes[index] == null) {
					attributes[index] = new Hashtable();
				}
				attributes[index].put(name, value);
			}
		}

		/**
		 * @return an object representing the named attributed for the numbered node.
		 * @param node the node being interrogated.
		 * @param name the name of the attribute of interest.
		 */
		public Object getAttribute(Node node, String name) {
			if (node instanceof AttributeNode) {
				return ((AttributeNode)node).getAttribute(name);
			} else {
				int index = getIndex(node);
				if (attributes == null || attributes[index] == null) {
					return null;
				}
				return attributes[index].get(name);
			}
		}

		/**
		 * make node with number num to root node
		 *
		 * @param num number of node
		 */
		public void reroot(int num)	{
			TreeUtils.reroot(this, findNode(num));
		}

		/**
		 * make provided node the root node
		 *
		 * @param node the node to make the root.
		 */
		public void reroot(Node node) {
			TreeUtils.reroot(this, node);
		}

// ========= IdGroup stuff ===============================
	public int getIdCount() {
		return getExternalNodeCount();
	}
	public Identifier getIdentifier(int i) {
		return getExternalNode(i).getIdentifier();
	}
	public void setIdentifier(int i, Identifier id) {
		getExternalNode(i).setIdentifier(id);
	}
	public int whichIdNumber(String s) {
		return IdGroupUtils.whichIdNumber(this,s);
	}

//========================================================

		// interface Report

		public void report(PrintWriter out)	{
			TreeUtils.report(this, out);
		}
	} //End of class Tree.AbstractTree
}
