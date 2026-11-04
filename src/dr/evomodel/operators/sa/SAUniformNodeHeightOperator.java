/*
 * SAUniformNodeHeightOperator.java
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
import dr.evomodel.tree.TreeModel;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * Draws the height of a random non-root, non-fake internal node uniformly between its oldest child and its parent.
 * Port of BEAST 2 SA's SAUniform; Hastings ratio 1.
 *
 * @author Giovanni Mazzeo
 */
public class SAUniformNodeHeightOperator extends SATreeOperator {

    public SAUniformNodeHeightOperator(TreeModel tree, double weight) {
        super(tree, null, weight);
    }

    @Override
    public double doOperation() {
        List<NodeRef> eligible = new ArrayList<>();
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getInternalNode(i);
            if (!tree.isRoot(node) && !isFake(node)) {
                eligible.add(node);
            }
        }
        if (eligible.isEmpty()) {
            return Double.NEGATIVE_INFINITY;
        }
        NodeRef node = eligible.get(MathUtils.nextInt(eligible.size()));
        double upper = h(tree.getParent(node));
        double lower = Math.max(h(tree.getChild(node, 0)), h(tree.getChild(node, 1)));
        tree.setNodeHeight(node, lower + MathUtils.nextDouble() * (upper - lower));
        return 0.0;
    }

    @Override
    public String getOperatorName() {
        return "saUniformNodeHeight(" + tree.getId() + ")";
    }
}
