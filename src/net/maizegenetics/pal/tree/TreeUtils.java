// TreeUtils.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.tree;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.LabelMapping;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.io.FormattedOutput;
import net.maizegenetics.pal.math.MersenneTwisterFast;

import java.io.PrintWriter;
import java.util.Hashtable;

/**
 * various utility functions on trees.
 *
 * @author Alexei Drummond
 * @author Korbinian Strimmer
 * @author Matthew Goode
 */
public class TreeUtils {

    /**
     * computes Robinson-Foulds (1981) distance between two trees
     *
     * @param t1 tree 1
     * @param t2 tree 2
     *
     * Definition: Assuming that t1 is the reference tree, let fn be the
     * false negatives, i.e. the number of edges in t1 missing in t2,
     * and fp the number of false positives, i.e. the number of edges
     * in t2 missing in t1.  The RF distance is then (fn + fp)/2
     */
    public static double getRobinsonFouldsDistance(Tree t1, Tree t2) {
        SplitSystem s1 = SplitUtils.getSplits(t1);

        return getRobinsonFouldsDistance(s1, t2);
    }

    /**
     * computes Robinson-Foulds (1981) distance between two trees
     *
     * @param s1 tree 1 (as represented by a SplitSystem)
     * @param t2 tree 2
     */
    public static double getRobinsonFouldsDistance(SplitSystem s1, Tree t2) {
        IdGroup idGroup = s1.getIdGroup();
        SplitSystem s2 = SplitUtils.getSplits(idGroup, t2);

        if (s1.getLabelCount() != s2.getLabelCount()) {
            throw new IllegalArgumentException("Number of labels must be the same!");
        }

        int ns1 = s1.getSplitCount();
        int ns2 = s2.getSplitCount();

        // number of splits in t1 missing in t2
        int fn = 0;
        for (int i = 0; i < ns1; i++) {
            if (!s2.hasSplit(s1.getSplit(i))) {
                fn++;
            }
        }

        // number of splits in t2 missing in t1
        int fp = 0;
        for (int i = 0; i < ns2; i++) {
            if (!s1.hasSplit(s2.getSplit(i))) {
                fp++;
            }
        }


        return 0.5 * ((double) fp + (double) fn);
    }

    /**
     * computes Robinson-Foulds (1981) distance between two trees
     * rescaled to a number between 0 and 1
     *
     * @param t1 tree 1
     * @param t2 tree 2
     */
    public static double getRobinsonFouldsRescaledDistance(Tree t1, Tree t2) {
        SplitSystem s1 = SplitUtils.getSplits(t1);

        return getRobinsonFouldsRescaledDistance(s1, t2);
    }

    /**
     * computes Robinson-Foulds (1981) distance between two trees
     * rescaled to a number between 0 and 1
     *
     * @param s1 tree 1 (as represented by a SplitSystem)
     * @param t2 tree 2
     */
    public static double getRobinsonFouldsRescaledDistance(SplitSystem s1, Tree t2) {
        return getRobinsonFouldsDistance(s1, t2) / (double) s1.getSplitCount();
    }
    private static MersenneTwisterFast random = new MersenneTwisterFast();

    /**
     * Returns a uniformly distributed random node from the tree, including
     * both internal and external nodes.
     */
    public static Node getRandomNode(Tree tree) {
        int index = random.nextInt(tree.getExternalNodeCount() + tree.getInternalNodeCount());
        if (index >= tree.getExternalNodeCount()) {
            return tree.getInternalNode(index - tree.getExternalNodeCount());
        } else {
            return tree.getExternalNode(index);
        }
    }

    /**
     * @return the first found node that has a certain name (as determined by the nodes Identifier)
     * in the tree defined by a root node.
     * @param tree The Tree supposidly containing such a named node
     * @param name The name of the node to find.
     * @return The node with the name, or null if no such node exists
     * @see Identifier, Node
     */
    public static final Node getNodeByName(Tree tree, String name) {
        return getNodeByName(tree.getRoot(), name);
    }

    /**
     * @return the first found node that has a certain name (as determined by the nodes Identifier)
     * in the tree defined by a root node.
     * @param root The root node of a tree
     * @param name The name of the node to find.
     * @return The node with the name, or null if no such node exists
     * @see Identifier, Node
     */
    public static final Node getNodeByName(Node root, String name) {
        if (root.getIdentifier().getName().equals(name)) {
            return root;
        }
        for (int i = 0; i < root.getChildCount(); i++) {
            Node result = getNodeByName(root.getChild(i), name);
            if (result != null) {
                return result;
            }
        }
        return null;
    }

    /**
     * Takes a tree (in mutation units) and returns a scaled version of it (in generation units).
     * @param mutationRateModel the mutation rate model used for scaling
     * and the desired units are expected substitutions then this scale
     * factor should be equal to the mutation rate.
     * @param newUnits the new units of the tree.
     */
    public static Tree mutationsToGenerations(Tree mutationTree, MutationRateModel muModel) {
        return scale(mutationTree, muModel, Units.GENERATIONS);
        /*if (muModel.getMutationRate(0.0) <= 0.0) {
        throw new IllegalArgumentException("Non-positive mutation rate is not permitted!");
        }

        SimpleTree tree = new SimpleTree(mutationTree);

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
        double oldHeight = tree.getExternalNode(i).getNodeHeight();
        tree.getExternalNode(i).setNodeHeight(muModel.getTime(oldHeight));
        }
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
        double oldHeight = tree.getInternalNode(i).getNodeHeight();
        tree.getInternalNode(i).setNodeHeight(muModel.getTime(oldHeight));
        }
        NodeUtils.heights2Lengths(tree.getRoot());
        tree.setUnits(Units.GENERATIONS);

        return tree;*/
    }

    /**
     * Takes a tree (in generation units) and returns a scaled version of it (in mutation units).
     * @param mutationRateModel the mutation rate model.
     * The mutation rate must be in units of mutations per site per generation.
     */
    public static Tree generationsToMutations(Tree generationTree, MutationRateModel muModel) {

        if (muModel.getUnits() != Units.GENERATIONS) {
            throw new IllegalArgumentException("Mutation rate must be per generation!");
        }

        return generationsToMutations(generationTree, muModel, 1.0);
    }

    /**
     * Takes a tree (in generation units) and returns a scaled version of it (in mutation units).
     * @param mutationRateModel the mutation rate model in calendar units.
     * @param generationTime the length of a generation in calendar units.
     * If the mutation rate is in mutations per site per year, then the
     * generation time will be in generations per year.
     */
    public static Tree generationsToMutations(Tree generationTree, MutationRateModel muModel, double generationTime) {

        if (generationTree.getUnits() != Units.GENERATIONS) {
            throw new IllegalArgumentException("Tree must be in units of generations!");
        }

        if (muModel.getMutationRate(0.0) <= 0.0) {
            throw new IllegalArgumentException("Non-positive mutation rate is not permitted!");
        }

        SimpleTree tree = new SimpleTree(generationTree);

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            double oldHeight = tree.getExternalNode(i).getNodeHeight();
            tree.getExternalNode(i).setNodeHeight(muModel.getExpectedSubstitutions(oldHeight) * generationTime);
        }
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            double oldHeight = tree.getInternalNode(i).getNodeHeight();
            tree.getInternalNode(i).setNodeHeight(muModel.getExpectedSubstitutions(oldHeight) * generationTime);
        }
        //Don't respect minimum branch lengths
        NodeUtils.heights2Lengths(tree.getRoot(), false);
        tree.setUnits(Units.EXPECTED_SUBSTITUTIONS);

        return tree;
    }

    /**
     * @deprecated use getScaled()
     */
    public static Tree scale(Tree oldTree, double rate, int newUnits) {
        return getScaled(oldTree, rate, newUnits);
    }

    /**
     * Takes a tree and returns a scaled version of it.
     * Scales a tree keeping old units
     * @param rate scale factor.
     *
     */
    public static final Tree getScaled(Tree oldTree, double rate) {
        return getScaled(oldTree, rate, oldTree.getUnits());
    }

    /**
     * Takes a tree and returns a scaled version of it.
     * @param rate scale factor. If the original tree is in generations
     * and the desired units are expected substitutions then this scale
     * factor should be equal to the mutation rate.
     * @param newUnits the new units of the tree.
     */
    public static final Tree getScaled(Tree oldTree, double rate, int newUnits) {
        SimpleTree tree = new SimpleTree(oldTree);
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            Node n = tree.getExternalNode(i);
            n.setNodeHeight(rate * n.getNodeHeight());
        }
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            Node n = tree.getInternalNode(i);
            n.setNodeHeight(rate * n.getNodeHeight());
        }
        NodeUtils.heights2Lengths(tree.getRoot());
        tree.setUnits(newUnits);
        return tree;
    }

    /**
     * @deprecated use getScaled()
     */
    public static Tree scale(Tree mutationRateTree, MutationRateModel muModel) {
        return getScaled(mutationRateTree, muModel);
    }

    /**
     * Takes a tree and returns a scaled version of it.
     * @param rate scale factor. If the original tree is in generations
     * and the desired units are expected substitutions then this scale
     * factor should be equal to the mutation rate.
     * @note resulting units is defined by muModel's units
     */
    public static Tree getScaled(Tree mutationRateTree, MutationRateModel muModel) {
        return getScaled(mutationRateTree, muModel, muModel.getUnits());
    }

    /**
     * @deprecated use getScaled()
     */
    public static Tree scale(Tree mutationRateTree, MutationRateModel muModel, int newUnits) {
        return getScaled(mutationRateTree, muModel, newUnits);
    }

    /**
     * Takes a tree and returns a scaled version of it.
     * @param rate scale factor. If the original tree is in generations
     * and the desired units are expected substitutions then this scale
     * factor should be equal to the mutation rate.
     * @param newUnits the new units of the tree. (Such as the mutationTree is measured in expected substitutions/newUnits)
     */
    public static Tree getScaled(Tree mutationRateTree, MutationRateModel muModel, int newUnits) {
        if (muModel.getMutationRate(0.0) <= 0.0) {
            throw new IllegalArgumentException("Non-positive mutation rate is not permitted!");
        }

        SimpleTree tree = new SimpleTree(mutationRateTree);
        if (newUnits == Units.EXPECTED_SUBSTITUTIONS) {
            //Changed for what I think is the correct behaviour for converting to Expected Substitutions
            for (int i = 0; i < tree.getExternalNodeCount(); i++) {
                double oldHeight = tree.getExternalNode(i).getNodeHeight();
                tree.getExternalNode(i).setNodeHeight(muModel.getExpectedSubstitutions(oldHeight));
            }
            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                double oldHeight = tree.getInternalNode(i).getNodeHeight();
                tree.getInternalNode(i).setNodeHeight(muModel.getExpectedSubstitutions(oldHeight));
            }
        } else {
            for (int i = 0; i < tree.getExternalNodeCount(); i++) {
                double oldHeight = tree.getExternalNode(i).getNodeHeight();
                tree.getExternalNode(i).setNodeHeight(muModel.getTime(oldHeight));
            }
            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                double oldHeight = tree.getInternalNode(i).getNodeHeight();
                tree.getInternalNode(i).setNodeHeight(muModel.getTime(oldHeight));
            }
        }
        NodeUtils.heights2Lengths(tree.getRoot());
        tree.setUnits(newUnits);
        return tree;
    }

    /**
     * Given a translation table where the keys are the current
     * identifier names and the values are the new identifier names,
     * this method replaces the current identifiers in the tree with new
     * identifiers.
     */
    public static void renameNodes(Tree tree, Hashtable table) {

        tree.createNodeList();

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            String newName =
                    (String) table.get(tree.getExternalNode(i).getIdentifier().getName());

            if (newName != null) {
                tree.getExternalNode(i).setIdentifier(new Identifier(newName));
            }
        }
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {



            String newName =
                    (String) table.get(tree.getInternalNode(i).getIdentifier().getName());

            if (newName != null) {
                tree.getInternalNode(i).setIdentifier(new Identifier(newName));
            }
        }
    }

    /**
     * Rotates branches by leaf count.
     * WARNING: assumes binary tree!
     */
    public static void rotateByLeafCount(Tree tree) {
        rotateByLeafCount(tree.getRoot());
    }

    /**
     * get list of the identifiers of the external nodes
     *
     * @return leaf identifier group
     */
    public static final IdGroup getLeafIdGroup(Tree tree) {
        tree.createNodeList();

        IdGroup labelList =
                new SimpleIdGroup(tree.getExternalNodeCount());

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            labelList.setIdentifier(i, tree.getExternalNode(i).getIdentifier());
        }

        return labelList;
    }

    /**
     * map external identifiers in the tree to a set of given identifiers
     * (which can be larger than the set of external identifiers but
     * must contain all of them)
     * <b> NOTE: </b> for efficiency it is assumed that the node lists of the tree are
     * correctly maintained. This should take effect everywhere soon.
     * Only methods that actually change a tree need to call createNodeList().
     *
     * @param idGroup an ordered group of identifiers
     *
     * @return list of links
     */
    public static final int[] mapExternalIdentifiers(IdGroup idGroup, Tree tree)
            throws IllegalArgumentException {

        int[] alias = new int[tree.getExternalNodeCount()];

        // Check whether for each label in tree there is
        // a correspondence in the given set of labels
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            alias[i] = idGroup.whichIdNumber(tree.getExternalNode(i).getIdentifier().getName());

            if (alias[i] == -1) {
                throw new IllegalArgumentException("Tree label " + tree.getExternalNode(i).getIdentifier()
                        + " not present in given set of labels");
            }
        }

        return alias;
    }

    /**
     * Labels the internal nodes of the tree using numbers starting from 0.
     * Skips numbers already used by external leaves.
     */
    public static final void labelInternalNodes(Tree tree) {

        int counter = 0;
        String pos = "0";

        IdGroup ids = getLeafIdGroup(tree);

        for (int i = 0; i < tree.getInternalNodeCount(); i++) {

            //if label already used find a better one
            while (ids.whichIdNumber(pos) >= 0) {
                counter += 1;
                pos = "" + counter;
            }
            tree.getInternalNode(i).setIdentifier(new Identifier(pos));
            counter += 1;
            pos = "" + counter;
        }
    }

    /**
     * Extracts a time order character data from a tree.
     */
    public static TimeOrderCharacterData extractTimeOrderCharacterData(Tree tree, int units) {

        tree.createNodeList();
        IdGroup identifiers = getLeafIdGroup(tree);
        TimeOrderCharacterData tocd = new TimeOrderCharacterData(identifiers, units);
        double[] times = new double[tree.getExternalNodeCount()];

        // WARNING: following code assumes that getLeafIdGroup
        //has same order as external node list.
        for (int i = 0; i < times.length; i++) {
            times[i] = tree.getExternalNode(i).getNodeHeight();
        }

        // this sets the ordinals as well
        tocd.setTimes(times, units);

        return tocd;
    }

    /**
     * print a this tree in New Hampshire format
     * (including distances and internal labels)
     *
     * @param out output stream
     */
    public static void printNH(Tree tree, PrintWriter out) {
        printNH(tree, out, true, true);
    }

    /**
     * print this tree in New Hampshire format
     *
     * @param out output stream
     * @param printLengths boolean variable determining whether
     *		branch lengths should be included in output
     * @param printInternalLabels boolean variable determining whether
     *		internal labels should be included in output
     */
    public static void printNH(Tree tree, PrintWriter out,
            boolean printLengths, boolean printInternalLabels) {

        NodeUtils.printNH(out, tree.getRoot(),
                printLengths, printInternalLabels);
        out.println(";");
    }

    /**
     * Roots a tree (that was previously unroot - ie 3 or more children at the
     * compsci tree root)
     * @param outgroupMembers the names of the nodes that form the outgroup.
     *      Multiple nodes will make the clade covering all outgroup nodes (and
     *      any others that fall with in that clade) form the outgroup.
     * @note if none of the outgroup members are actually in the tree, or the outgroup clade is
     * the whole tree, the result is just an unrooted clone of the input tree.
     */
//	public static final Tree getRooted(Tree unrooted, String[] outgroupMembers) {
//		Tree t2 = new SimpleTree(unrooted);
//		Node[] nodes = NodeUtils.findByIdentifier(t2.getRoot(),outgroupMembers);
//		if(nodes==null) {
//			return t2;
//		}
//		Node common = (nodes.length==1 ? nodes[0] : NodeUtils.getFirstCommonAncestor(nodes));
//		if(common==null) {	return t2;	}
//		if(common==t2.getRoot()) { return t2; }
//		//TreeUtils.reroot(t2,common);
//		Node newRoot = NodeUtils.rootAbove(common);
//		NodeUtils.exchangeInfo(newRoot,t2.getRoot());
//		NodeUtils.lengths2Heights(newRoot);
//		t2.setRoot(newRoot);
//
//		return t2;
//	}
    /**
     * @return a tree that has been re-rooted at the given
     * internal node. The original root of the tree must have had at least 3 children.
     */
    public static void reroot(Tree tree, Node node) {
        reroot(node);
        tree.setRoot(node);
    }

    /**
     * Return node with the highest average minimum and maximum distance.
     * Ed thinks this is roughly the midpoint - need to look this up
     * @param tree
     * @return
     */
    public static Node findMidpointNode(Tree tree) {
        Node midNode = null;
        double maxAvgDist = -1;
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            Node iNode = tree.getInternalNode(i);
            double[] d = NodeUtils.getPathLengthInfo(iNode);
            double avg = (d[0] + d[1]) / 2.0;
            if (maxAvgDist < avg) {
                maxAvgDist = avg;
                midNode = iNode;
            }
        }
        return midNode;
    }

    /*
     * compute distance of external node a to all other leaves
     * (computational complexity of this method is only O(n), following
     * D.Bryant and P. Wadell. 1998. MBE 15:1346-1359)
     *
     * @param tree tree
     * @param a node
     * @param dist array for the node-to-node distance distances
     * @param idist array for the distance between a and all internal nodes
     * @param countEdges boolean variable deciding whether the actual
     *                   branch lengths are used in computing the distance
     *                   or whether simply all edges larger or equal a certain
     *                   threshold length are counted (each with weight 1.0)
     * @param epsilon    minimum branch length for a which an edge is counted
     */
    public static void computeAllDistances(Tree tree,
            int a, double[] dist, double[] idist,
            boolean countEdges, double epsilon) {
        tree.createNodeList();

        dist[a] = 0.0;

        Node node = tree.getExternalNode(a);

        computeNodeDist(node, node.getParent(), dist, idist, countEdges, epsilon);
    }

    private static void computeNodeDist(Node origin, Node center,
            double[] dist, double[] idist,
            boolean countEdges, double epsilon) {
        int indexCenter = center.getNumber();
        int indexOrigin = origin.getNumber();
        double[] distCenter;
        double[] distOrigin;
        if (center.isLeaf()) {
            distCenter = dist;
        } else {
            distCenter = idist;
        }
        if (origin.isLeaf()) {
            distOrigin = dist;
        } else {
            distOrigin = idist;
        }

        double len;
        double tmp;
        if (origin.getParent() == center) {
            // center is parent of origin
            tmp = origin.getBranchLength();
        } else {
            // center is child of origin
            tmp = center.getBranchLength();
        }


        if (countEdges) // count all edges >= epsilon
        {
            if (tmp < epsilon) {
                len = 0.0;
            } else {
                len = 1.0;
            }
        } else // use branch lengths
        {
            len = tmp;
        }


        distCenter[indexCenter] = distOrigin[indexOrigin] + len;

        if (!center.isLeaf()) {
            for (int i = 0; i < center.getChildCount(); i++) {
                Node c = center.getChild(i);

                if (c != origin) {
                    computeNodeDist(center, c, dist, idist, countEdges, epsilon);
                }
            }

            if (!center.isRoot()) {
                Node p = center.getParent();

                if (p != origin) {
                    computeNodeDist(center, p, dist, idist, countEdges, epsilon);
                }
            }
        }
    }
    private static Node[] path;

    /**
     * compute distance between two external nodes
     *
     * @param tree tree
     * @param a external node 1
     * @param b external node 2
     *
     * @return distance between node a and b
     */
    public static final double computeDistance(Tree tree, int a, int b) {
        tree.createNodeList();
        int maxLen = tree.getInternalNodeCount() + 1;
        if (path == null || path.length < maxLen) {
            path = new Node[maxLen];
        }

        // len might be different from path.length
        int len = findPath(tree, a, b);

        double dist = 0.0;
        for (int i = 0; i < len; i++) {
            dist += path[i].getBranchLength();
        }

        return dist;
    }

    // Find path between external nodes a and b
    // After calling this method path contains all nodes
    // with edges lying between a and b (including a and b)
    // (note that the node lying on the intersection of a-root
    // and b-root is NOT contained because this node does
    // not contain a branch of the path)
    // The length of the path is also returned
    private static final int findPath(Tree tree, int a, int b) {
        // clean path
        for (int i = 0; i < path.length; i++) {
            path[i] = null;
        }
        // path from node a to root
        Node node = tree.getExternalNode(a);
        int len = 0;
        path[len] = node;
        len++;
        while (!node.isRoot()) {
            node = node.getParent();
            path[len] = node;
            len++;
        }

        // find intersection with path from node b to root
        Node stopNode = null;
        node = tree.getExternalNode(b);
        while (!node.isRoot()) {
            node = node.getParent();
            int pos = findInPath(node);

            if (pos != -1) {
                len = pos;
                stopNode = node;
                break;
            }
        }

        // fill rest of path
        node = tree.getExternalNode(b);
        path[len] = node;
        len++;
        node = node.getParent();
        while (node != stopNode) {
            path[len] = node;
            len++;
            node = node.getParent();
        }

        // clean rest
        for (int i = len; i < path.length; i++) {
            path[i] = null;
        }

        return len;
    }

    private static final int findInPath(Node node) {
        for (int i = 0; i < path.length; i++) {
            if (path[i] == node) {
                return i;
            } else if (path[i] == null) {
                return -1;
            }
        }

        return -1;
    }

    /**
     * Rotates branches by leaf count.
     * WARNING: assumes binary tree!
     */
    private static void rotateByLeafCount(Node node) {

        if (!node.isLeaf()) {
            if (NodeUtils.getLeafCount(node.getChild(0))
                    > NodeUtils.getLeafCount(node.getChild(1))) {
                Node temp = node.getChild(0);
                node.removeChild(0);
                node.addChild(temp);
            }
            for (int i = 0; i < node.getChildCount(); i++) {
                rotateByLeafCount(node.getChild(i));
            }
        }
    }

    public static void report(Tree tree, PrintWriter out) {
        printASCII(tree, out);
        out.println();
        branchInfo(tree, out);
        out.println();
        heightInfo(tree, out);
    }
    private static FormattedOutput format;
    private static double proportion;
    private static int minLength;
    private static boolean[] umbrella;
    private static int[] position;
    private static int numExternalNodes;
    private static int numInternalNodes;
    private static int numBranches;

    // Print picture of current tree in ASCII
    private static void printASCII(Tree tree, PrintWriter out) {
        format = FormattedOutput.getInstance();

        tree.createNodeList();

        numExternalNodes = tree.getExternalNodeCount();
        numInternalNodes = tree.getInternalNodeCount();
        numBranches = numInternalNodes + numExternalNodes - 1;

        umbrella = new boolean[numExternalNodes];
        position = new int[numExternalNodes];

        minLength = (Integer.toString(numBranches)).length() + 1;

        int MAXCOLUMN = 40;
        Node root = tree.getRoot();
        if (root.getNodeHeight() == 0.0) {
            NodeUtils.lengths2Heights(root);
        }
        proportion = (double) MAXCOLUMN / root.getNodeHeight();

        for (int n = 0; n < numExternalNodes; n++) {
            umbrella[n] = false;
        }

        position[0] = 1;
        for (int i = root.getChildCount() - 1; i > -1; i--) {
            printNodeInASCII(out, root.getChild(i), 1, i, root.getChildCount());
            if (i != 0) {
                putCharAtLevel(out, 0, '|');
                out.println();
            }
        }
    }

    // Print branch information
    private static void branchInfo(Tree tree, PrintWriter out) {

        //
        // CALL PRINTASCII FIRST !!!
        //

        // check if some SE values differ from the default zero
        boolean showSE = false;
        for (int i = 0; i < numExternalNodes && showSE == false; i++) {
            if (tree.getExternalNode(i).getBranchLengthSE() != 0.0) {
                showSE = true;
            }
            if (i < numInternalNodes - 1) {
                if (tree.getInternalNode(i).getBranchLengthSE() != 0.0) {
                    showSE = true;
                }
            }
        }

        format.displayIntegerWhite(out, numExternalNodes);
        out.print("   Length    ");
        if (showSE) {
            out.print("S.E.      ");
        }
        out.print("Label     ");
        if (numInternalNodes > 1) {
            format.displayIntegerWhite(out, numBranches);
            out.print("        Length    ");
            if (showSE) {
                out.print("S.E.      ");
            }
            out.print("Label");
        }
        out.println();

        for (int i = 0; i < numExternalNodes; i++) {
            format.displayInteger(out, i + 1, numExternalNodes);
            out.print("   ");
            format.displayDecimal(out, tree.getExternalNode(i).getBranchLength(), 5);
            out.print("   ");
            if (showSE) {
                format.displayDecimal(out, tree.getExternalNode(i).getBranchLengthSE(), 5);
                out.print("   ");
            }
            format.displayLabel(out, tree.getExternalNode(i).getIdentifier().getName(), 10);

            if (i < numInternalNodes - 1) {
                format.multiplePrint(out, ' ', 5);
                format.displayInteger(out, i + 1 + numExternalNodes, numBranches);
                out.print("   ");
                format.displayDecimal(out, tree.getInternalNode(i).getBranchLength(), 5);
                out.print("   ");
                if (showSE) {
                    format.displayDecimal(out, tree.getInternalNode(i).getBranchLengthSE(), 5);
                    out.print("   ");
                }
                format.displayLabel(out, tree.getInternalNode(i).getIdentifier().getName(), 10);
            }

            out.println();
        }
    }

    // Print height information
    private static void heightInfo(Tree tree, PrintWriter out) {
        //
        // CALL PRINTASCII FIRST
        //

        if (tree.getRoot().getNodeHeight() == 0.0) {
            NodeUtils.lengths2Heights(tree.getRoot());
        }

        // check if some SE values differ from the default zero
        //boolean showSE = false;
        //for (int i = 0; i < numInternalNodes && showSE == false; i++)
        //{
        //	if (tree.getInternalNode(i).getNodeHeightSE() != 0.0)
        //	{
        //		showSE = true;
        //	}
        //}

        format.displayIntegerWhite(out, numExternalNodes);
        out.print("   Height    ");
        format.displayIntegerWhite(out, numBranches);
        out.print("        Height    ");
        //if (showSE) out.print("S.E.");

        out.println();

        for (int i = 0; i < numExternalNodes; i++) {
            format.displayInteger(out, i + 1, numExternalNodes);
            out.print("   ");
            format.displayDecimal(out, tree.getExternalNode(i).getNodeHeight(), 7);
            out.print("   ");

            if (i < numInternalNodes) {
                format.multiplePrint(out, ' ', 5);

                if (i == numInternalNodes - 1) {
                    out.print("R");
                    format.multiplePrint(out, ' ', Integer.toString(numBranches).length() - 1);
                } else {
                    format.displayInteger(out, i + 1 + numExternalNodes, numBranches);
                }

                out.print("   ");
                format.displayDecimal(out, tree.getInternalNode(i).getNodeHeight(), 7);
                out.print("   ");
                //if (showSE)
                //{
                //	format.displayDecimal(out, tree.getInternalNode(i).getNodeHeightSE(), 7);
                //}
            }

            out.println();
        }
    }

    private static void printNodeInASCII(PrintWriter out, Node node, int level, int m, int maxm) {
        position[level] = (int) (node.getBranchLength() * proportion);

        if (position[level] < minLength) {
            position[level] = minLength;
        }

        if (node.isLeaf()) // external branch
        {
            if (m == maxm - 1) {
                umbrella[level - 1] = true;
            }

            printlnNodeWithNumberAndLabel(out, node, level);

            if (m == 0) {
                umbrella[level - 1] = false;
            }
        } else // internal branch
        {
            for (int n = node.getChildCount() - 1; n > -1; n--) {
                printNodeInASCII(out, node.getChild(n), level + 1, n, node.getChildCount());

                if (m == maxm - 1 && n == node.getChildCount() / 2) {
                    umbrella[level - 1] = true;
                }

                if (n != 0) {
                    if (n == node.getChildCount() / 2) {
                        printlnNodeWithNumberAndLabel(out, node, level);
                    } else {
                        for (int i = 0; i < level + 1; i++) {
                            if (umbrella[i]) {
                                putCharAtLevel(out, i, '|');
                            } else {
                                putCharAtLevel(out, i, ' ');
                            }
                        }
                        out.println();
                    }
                }

                if (m == 0 && n == node.getChildCount() / 2) {
                    umbrella[level - 1] = false;
                }
            }
        }
    }

    private static void printlnNodeWithNumberAndLabel(PrintWriter out, Node node, int level) {
        for (int i = 0; i < level - 1; i++) {
            if (umbrella[i]) {
                putCharAtLevel(out, i, '|');
            } else {
                putCharAtLevel(out, i, ' ');
            }
        }

        putCharAtLevel(out, level - 1, '+');

        int branchNumber;
        if (node.isLeaf()) {
            branchNumber = node.getNumber() + 1;
        } else {
            branchNumber = node.getNumber() + 1 + numExternalNodes;
        }

        String numberAsString = Integer.toString(branchNumber);

        int numDashs = position[level] - numberAsString.length();
        for (int i = 0; i < numDashs; i++) {
            out.print('-');
        }
        out.print(numberAsString);

        if (node.isLeaf()) {
            out.println(" " + node.getIdentifier());
        } else {
            if (!node.getIdentifier().equals(Identifier.ANONYMOUS)) {
                out.print("(" + node.getIdentifier() + ")");
            }
            out.println();
        }
    }

    private static void putCharAtLevel(PrintWriter out, int level, char c) {
        int n = position[level] - 1;
        for (int i = 0; i < n; i++) {
            out.print(' ');
        }
        out.print(c);
    }

    /**
     * make given node the root node
     *
     * @param node new root node
     */
    private static void reroot(Node node) {
        if (node.isRoot() || node.isLeaf()) {
            return;
        }

        if (!node.getParent().isRoot()) {
            reroot(node.getParent());
        }

        // Now the parent of node is root

        if (node.getParent().getChildCount() < 3) {
            // Rerooting not possible
            return;
        }

        // Exchange branch label, length et cetera
        NodeUtils.exchangeInfo(node.getParent(), node);

        // Rearrange topology
        Node parent = node.getParent();
        NodeUtils.removeChild(parent, node);
        node.addChild(parent);
    }

    /**
     * Generates a tree which is identical to baseTree but has attributes (defined by attributeName)
     * at all internal nodes excluding the root node signifying (as a value between 0 and 100) the bootstrap
     * support by clade (that is the proportion of replicates that produce the sub clade under that node)
     * @note assumes all alternative trees have the exact same set of labels
     * @deprecated Use getReplicateCladeSupport instead
     */
    public static final Tree getBootstrapSupportByCladeTree(String attributeName, Tree baseTree, Tree[] alternativeTrees) {
        SimpleTree result = new SimpleTree(baseTree);
        IdGroup ids = TreeUtils.getLeafIdGroup(baseTree);
        SplitSystem baseSystem = SplitUtils.getSplits(ids, baseTree);
        boolean[][] baseVector = baseSystem.getSplitVector();
        int[] supportCount = new int[baseVector.length];
        for (int i = 0; i < alternativeTrees.length; i++) {
            SplitSystem alternativeSystem = SplitUtils.getSplits(ids, alternativeTrees[i]);
            for (int j = 0; j < baseVector.length; j++) {
                if (alternativeSystem.hasSplit(baseVector[j])) {
                    supportCount[j]++;
                }
            }
        }
        for (int i = 0; i < supportCount.length; i++) {
            int support = (int) (supportCount[i] * 100 / (double) alternativeTrees.length);
            result.setAttribute(
                    result.getInternalNode(i),
                    attributeName,
                    new Integer(support));
        }
        return result;
    }

    /**
     * Create a new tree such that the labels are redifined from a base tree in such a manner:
     *  For each leaf label
     * <ol>
     *  <li> If the base label is not a number the new label is just the original label </i>
     *  <li> If the base label is a number the new label appropriately index label from a set identifiers </i>
     * </ol>
     * @param baseTree The base tree
     * @param ids The set of identifiers
     * @return A relabelled tree, or the input tree if no numbered leaves.
     */
    public static final Tree getNumberRelabelledTree(Tree baseTree, IdGroup ids) {
        LabelMapping lm = new LabelMapping();
        int minIndex = Integer.MAX_VALUE;
        for (int i = 0; i < baseTree.getIdCount(); i++) {
            String treeLabel = baseTree.getIdentifier(i).getName();
            try {
                int number = Integer.parseInt(treeLabel);
                if (number < minIndex) {
                    minIndex = number;
                }
            } catch (NumberFormatException e) {
            }
        }
        int changes = 0;
        for (int i = 0; i < baseTree.getIdCount(); i++) {
            String treeLabel = baseTree.getIdentifier(i).getName();
            try {
                int number = Integer.parseInt(treeLabel) - minIndex;
                if (number < ids.getIdCount()) {
                    lm.addMapping(treeLabel, ids.getIdentifier(number).getName());
                    changes++;
                }
            } catch (NumberFormatException e) {
            }
        }
        if (changes == 0) {
            return baseTree;
        }
        return new SimpleTree(baseTree, lm);
    }
}
