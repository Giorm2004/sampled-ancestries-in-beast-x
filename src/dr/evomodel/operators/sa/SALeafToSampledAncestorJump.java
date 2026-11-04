/*
 * SALeafToSampledAncestorJump.java
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
 * Turns a random tip into a sampled ancestor (its parent drops onto it) or a sampled ancestor into an ordinary tip
 * (its fake parent moves up uniformly on the branch above, or by Exp(1) above the root). Port of BEAST 2 SA's
 * LeafToSampledAncestorJump; log Hastings ratio log(newRange / oldRange).
 *
 * @author Giovanni Mazzeo
 */
public class SALeafToSampledAncestorJump extends SATreeOperator {

    public SALeafToSampledAncestorJump(TreeModel tree, Parameter removalProbability, double weight) {
        super(tree, removalProbability, weight);
    }

    @Override
    public double doOperation() {
        NodeRef leaf = tree.getExternalNode(MathUtils.nextInt(tree.getExternalNodeCount()));
        NodeRef parent = tree.getParent(leaf);
        double newHeight, newRange, oldRange;

        if (isSA(leaf)) {
            oldRange = 1.0;
            if (tree.isRoot(parent)) {
                double x = MathUtils.nextExponential(1.0);
                newHeight = h(parent) + x;
                newRange = Math.exp(x);
            } else {
                newRange = h(tree.getParent(parent)) - h(parent);
                newHeight = h(parent) + MathUtils.nextDouble() * newRange;
            }
        } else {
            newRange = 1.0;
            if (h(SampledAncestorTrees.otherChild(tree, parent, leaf)) >= h(leaf)) {
                return Double.NEGATIVE_INFINITY;    // the sibling clade must lie below the new sampled ancestor
            }
            oldRange = tree.isRoot(parent) ? Math.exp(h(parent) - h(leaf)) : h(tree.getParent(parent)) - h(leaf);
            newHeight = h(leaf);
        }
        tree.setNodeHeight(parent, newHeight);
        if (removalCheck() == Double.NEGATIVE_INFINITY) {
            return Double.NEGATIVE_INFINITY;
        }
        return Math.log(newRange / oldRange);
    }

    @Override
    public String getOperatorName() {
        return "saLeafToSampledAncestorJump(" + tree.getId() + ")";
    }
}
