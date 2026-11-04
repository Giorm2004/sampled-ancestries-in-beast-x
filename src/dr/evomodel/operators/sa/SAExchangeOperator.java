/*
 * SAExchangeOperator.java
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
import dr.math.MathUtils;

/**
 * Narrow or wide exchange that never moves a sampled ancestor off its fake node. Port of BEAST 2 SA's SAExchange;
 * Hastings ratio 1.
 *
 * @author Giovanni Mazzeo
 */
public class SAExchangeOperator extends SATreeOperator {

    private final boolean narrow;

    public SAExchangeOperator(TreeModel tree, boolean narrow, double weight) {
        super(tree, null, weight);
        this.narrow = narrow;
    }

    @Override
    public double doOperation() {
        final int nodeCount = tree.getNodeCount();
        if (nodeCount == 3 && isFake(tree.getRoot())) {
            return Double.NEGATIVE_INFINITY;
        }
        return narrow ? narrow(nodeCount) : wide(nodeCount);
    }

    private double narrow(int nodeCount) {
        boolean any = false;                              // guard the rejection loop
        for (int k = 0; k < nodeCount && !any; k++) {
            NodeRef n = tree.getNode(k);
            any = !tree.isRoot(n) && !tree.isRoot(tree.getParent(n)) && !isSA(n);
        }
        if (!any) {
            return Double.NEGATIVE_INFINITY;
        }
        NodeRef i;
        do {
            i = tree.getNode(MathUtils.nextInt(nodeCount));
        } while (tree.isRoot(i) || tree.isRoot(tree.getParent(i)) || isSA(i));
        final NodeRef iParent = tree.getParent(i);
        final NodeRef iGrandParent = tree.getParent(iParent);
        final NodeRef iUncle = SampledAncestorTrees.otherChild(tree, iGrandParent, iParent);
        if (h(iUncle) < h(iParent)) {
            exchangeNodes(tree, i, iUncle, iParent, iGrandParent);
            return 0.0;
        }
        return Double.NEGATIVE_INFINITY;
    }

    private double wide(int nodeCount) {
        NodeRef i, j;
        do {
            i = tree.getNode(MathUtils.nextInt(nodeCount));
        } while (tree.isRoot(i) || isSA(i));
        do {
            j = tree.getNode(MathUtils.nextInt(nodeCount));
        } while (j == i || tree.isRoot(j) || isSA(j));
        final NodeRef iP = tree.getParent(i);
        final NodeRef jP = tree.getParent(j);
        if (iP != jP && i != jP && j != iP && h(j) < h(iP) && h(i) < h(jP)) {
            exchangeNodes(tree, i, j, iP, jP);
            return 0.0;
        }
        return Double.NEGATIVE_INFINITY;
    }

    @Override
    public String getOperatorName() {
        return (narrow ? "saNarrowExchange(" : "saWideExchange(") + tree.getId() + ")";
    }
}
