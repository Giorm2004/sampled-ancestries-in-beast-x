/*
 * SampledAncestorTrees.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

/**
 * Sampled-ancestor trees in the representation of BEAST 2's SA package (Gavryushkina et al. 2014): a sampled
 * ancestor is a tip at height > 0 on a zero-length branch, and its parent (at the same height) is a "fake" node.
 * Operators keep the heights of a sampled ancestor and its fake node exactly equal.
 *
 * @author Giovanni Mazzeo
 */
public final class SampledAncestorTrees {

    private SampledAncestorTrees() { }

    public static boolean isSampledAncestor(Tree tree, NodeRef node) {
        return tree.isExternal(node) && !tree.isRoot(node) && tree.getNodeHeight(node) > 0.0
                && tree.getNodeHeight(tree.getParent(node)) == tree.getNodeHeight(node);
    }

    /** An internal node with a sampled-ancestor child. */
    public static boolean isFake(Tree tree, NodeRef node) {
        if (tree.isExternal(node)) {
            return false;
        }
        for (int i = 0; i < tree.getChildCount(node); i++) {
            if (isSampledAncestor(tree, tree.getChild(node, i))) {
                return true;
            }
        }
        return false;
    }

    public static int countSampledAncestors(Tree tree) {
        int count = 0;
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            if (isSampledAncestor(tree, tree.getExternalNode(i))) {
                count++;
            }
        }
        return count;
    }

    public static NodeRef otherChild(Tree tree, NodeRef parent, NodeRef child) {
        NodeRef first = tree.getChild(parent, 0);
        return first == child ? tree.getChild(parent, 1) : first;
    }

    /** For a fake node, its child that is not the sampled ancestor. */
    public static NodeRef nonSampledAncestorChild(Tree tree, NodeRef fake) {
        NodeRef first = tree.getChild(fake, 0);
        return isSampledAncestor(tree, first) ? tree.getChild(fake, 1) : first;
    }
}
