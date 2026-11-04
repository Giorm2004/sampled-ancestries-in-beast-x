/*
 * SATreeOperator.java
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
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.tree.SampledAncestorTrees;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;

/**
 * Base for the sampled-ancestor tree operators, ports of BEAST 2's SA package (v2.1.1, sa.evolution.operators) with
 * the same proposals and Hastings ratios. Sampled ancestors are tips on zero-length branches (SampledAncestorTrees);
 * every operator keeps a sampled ancestor and its fake parent at exactly the same height. TreeModel.endTreeEdit does
 * not check heights, so each operator checks the proposed heights itself and returns -infinity (before editing)
 * for an invalid tree.
 *
 * @author Giovanni Mazzeo
 */
public abstract class SATreeOperator extends AbstractTreeOperator {

    protected final TreeModel tree;
    /** Optional removal probability: no sampled ancestors when it is 1 (as in BEAST 2's operators). */
    protected final Parameter removalProbability;

    protected SATreeOperator(TreeModel tree, Parameter removalProbability, double weight) {
        this.tree = tree;
        this.removalProbability = removalProbability;
        setWeight(weight);
    }

    protected final double h(NodeRef node) {
        return tree.getNodeHeight(node);
    }

    protected final boolean isSA(NodeRef node) {
        return SampledAncestorTrees.isSampledAncestor(tree, node);
    }

    protected final boolean isFake(NodeRef node) {
        return SampledAncestorTrees.isFake(tree, node);
    }

    /** -infinity if the tree has sampled ancestors but the removal probability is 1, else 0. */
    protected final double removalCheck() {
        return removalProbability != null && removalProbability.getParameterValue(0) == 1.0
                && SampledAncestorTrees.countSampledAncestors(tree) > 0 ? Double.NEGATIVE_INFINITY : 0.0;
    }

    public String getPerformanceSuggestion() {
        return "";
    }
}
