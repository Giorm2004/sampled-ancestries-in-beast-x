/*
 * SAScaleOperator.java
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
 * Scales every non-fake internal node height (fake nodes stay on their sampled ancestors; tips do not move), or the
 * root only. The scale is uniform on [f, 1/f] as in BEAST 2's ScaleOperator. Port of BEAST 2 SA's SAScaleOperator:
 * log Hastings ratio (dof - 2) log(s) with dof the number of non-fake internal nodes, or -log(s) for the root.
 *
 * @author Giovanni Mazzeo
 */
public class SAScaleOperator extends SATreeOperator {

    private final double scaleFactor;
    private final boolean rootOnly;

    public SAScaleOperator(TreeModel tree, double scaleFactor, boolean rootOnly, double weight) {
        super(tree, null, weight);
        if (!(scaleFactor > 0.0 && scaleFactor < 1.0)) {
            throw new IllegalArgumentException("scaleFactor must be in (0, 1)");
        }
        this.scaleFactor = scaleFactor;
        this.rootOnly = rootOnly;
    }

    static double drawScale(double f) {
        return f + MathUtils.nextDouble() * (1.0 / f - f);
    }

    @Override
    public double doOperation() {
        final double scale = drawScale(scaleFactor);
        if (rootOnly) {
            NodeRef root = tree.getRoot();
            if (isFake(root)) {
                return Double.NEGATIVE_INFINITY;
            }
            double newHeight = h(root) * scale;
            if (newHeight < Math.max(h(tree.getChild(root, 0)), h(tree.getChild(root, 1)))) {
                return Double.NEGATIVE_INFINITY;
            }
            tree.setNodeHeight(root, newHeight);
            return -Math.log(scale);
        }
        int dof = scaleTree(tree, scale);
        return dof < 0 ? Double.NEGATIVE_INFINITY : Math.log(scale) * (dof - 2);
    }

    /**
     * Multiplies the height of every non-fake internal node by scale, if the result is a valid tree (fake nodes and
     * tips held). Returns the number of nodes scaled, or -1 (tree unchanged) if some node would end below a child.
     */
    static int scaleTree(TreeModel tree, double scale) {
        final int n = tree.getNodeCount();
        double[] height = new double[n];
        boolean[] scaled = new boolean[n];
        int dof = 0;
        for (int k = 0; k < n; k++) {
            NodeRef node = tree.getNode(k);
            height[k] = tree.getNodeHeight(node);
            if (!tree.isExternal(node) && !SampledAncestorTrees.isFake(tree, node)) {
                height[k] *= scale;
                scaled[k] = true;
                dof++;
            }
        }
        for (int k = 0; k < n; k++) {
            NodeRef node = tree.getNode(k);
            if (!tree.isRoot(node) && height[k] > height[tree.getParent(node).getNumber()]) {
                return -1;
            }
        }
        tree.beginTreeEdit();
        for (int k = 0; k < n; k++) {
            if (scaled[k]) {
                tree.setNodeHeight(tree.getNode(k), height[k]);
            }
        }
        tree.endTreeEdit();
        return dof;
    }

    @Override
    public String getOperatorName() {
        return (rootOnly ? "saRootScale(" : "saTreeScale(") + tree.getId() + ")";
    }
}
