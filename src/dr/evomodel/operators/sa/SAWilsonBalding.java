/*
 * SAWilsonBalding.java
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

package dr.evomodel.operators.sa;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.SampledAncestorTrees;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

/**
 * Wilson-Balding (subtree prune and regraft) on sampled-ancestor trees: the pruned subtree may be regrafted onto a
 * branch (new parent height uniform on it, or Exp(1) above the root) or onto a tip, which then becomes a sampled
 * ancestor; a fake parent likewise frees its sampled ancestor. Port of BEAST 2 SA's SAWilsonBalding; log Hastings
 * ratio log(|oldDimension / newDimension * newRange / oldRange|), dimension = nodes - sampled ancestors - 1.
 *
 * @author Giovanni Mazzeo
 */
public class SAWilsonBalding extends SATreeOperator {

    public SAWilsonBalding(TreeModel tree, Parameter removalProbability, double weight) {
        super(tree, removalProbability, weight);
    }

    @Override
    public double doOperation() {
        final int nodeCount = tree.getNodeCount();
        final int leafCount = tree.getExternalNodeCount();

        NodeRef i;
        do {
            i = tree.getNode(MathUtils.nextInt(nodeCount));
        } while (tree.isRoot(i) || isSA(i));
        final NodeRef iP = tree.getParent(i);
        final NodeRef CiP = SampledAncestorTrees.otherChild(tree, iP, i);
        if (tree.isRoot(iP) && h(CiP) <= h(i)) {
            return Double.NEGATIVE_INFINITY;
        }

        NodeRef j, jP;
        double newParentHeight;
        boolean attachingToLeaf, adjacentEdge;
        do {
            adjacentEdge = false;
            int nodeNumber = MathUtils.nextInt(nodeCount + leafCount);
            if (nodeNumber < nodeCount) {
                j = tree.getNode(nodeNumber);
                jP = tree.isRoot(j) ? null : tree.getParent(j);
                newParentHeight = jP != null ? h(jP) : Double.POSITIVE_INFINITY;
                if (!isSA(CiP)) {
                    adjacentEdge = (CiP == j || iP == j);
                }
                attachingToLeaf = false;
            } else {
                j = tree.getExternalNode(nodeNumber - nodeCount);
                jP = tree.getParent(j);
                newParentHeight = h(j);
                attachingToLeaf = true;
            }
        } while (isSA(j) || newParentHeight <= h(i) || i == j || adjacentEdge);

        if (attachingToLeaf && iP == j) {
            return Double.NEGATIVE_INFINITY;
        }
        if (jP != null && jP == i) {
            return Double.NEGATIVE_INFINITY;
        }

        final int oldDimension = nodeCount - SampledAncestorTrees.countSampledAncestors(tree) - 1;

        double newRange, newAge;
        if (attachingToLeaf) {
            newRange = 1.0;
            newAge = h(j);
        } else if (jP != null) {
            double newMinAge = Math.max(h(i), h(j));
            newRange = h(jP) - newMinAge;
            newAge = newMinAge + MathUtils.nextDouble() * newRange;
        } else {
            double x = MathUtils.nextExponential(1.0);
            newRange = Math.exp(x);
            newAge = h(j) + x;
        }

        final NodeRef PiP = tree.isRoot(iP) ? null : tree.getParent(iP);
        double oldRange;
        if (isSA(CiP)) {
            oldRange = 1.0;
        } else {
            double oldMinAge = Math.max(h(i), h(CiP));
            oldRange = PiP != null ? h(PiP) - oldMinAge : Math.exp(h(iP) - oldMinAge);
        }

        tree.beginTreeEdit();
        if (iP != j && CiP != j) {
            if (jP == null) {                         // j is the root (so iP is not): iP becomes the root
                tree.removeChild(iP, CiP);
                tree.removeChild(PiP, iP);
                tree.addChild(iP, j);
                tree.addChild(PiP, CiP);
                tree.setRoot(iP);
            } else if (PiP == null) {                 // iP is the root: CiP becomes the root
                tree.removeChild(jP, j);
                tree.removeChild(iP, CiP);
                tree.addChild(iP, j);
                tree.addChild(jP, iP);
                tree.setRoot(CiP);
            } else {
                tree.removeChild(jP, j);
                tree.removeChild(iP, CiP);
                tree.removeChild(PiP, iP);
                tree.addChild(iP, j);
                tree.addChild(jP, iP);
                tree.addChild(PiP, CiP);
            }
        }
        tree.setNodeHeight(iP, newAge);
        tree.endTreeEdit();

        if (removalCheck() == Double.NEGATIVE_INFINITY) {
            return Double.NEGATIVE_INFINITY;
        }
        final int newDimension = nodeCount - SampledAncestorTrees.countSampledAncestors(tree) - 1;
        return Math.log(Math.abs((double) oldDimension / newDimension * newRange / oldRange));
    }

    @Override
    public String getOperatorName() {
        return "saWilsonBalding(" + tree.getId() + ")";
    }
}
