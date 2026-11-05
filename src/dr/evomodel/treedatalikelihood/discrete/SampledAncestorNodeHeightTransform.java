/*
 * SampledAncestorNodeHeightTransform.java
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.SampledAncestorTrees;
import dr.evomodel.tree.TreeModel;
import dr.util.Transform;
import dr.xml.Reportable;

/**
 * Maps the heights of all internal nodes (root included, in NodeHeightProxyParameter order: index i is node
 * number externalCount + i) of a possibly sampled-ancestor tree to the real line, for HMC on node heights.
 *
 * With a(v) the largest tip height below v (a sampled-ancestor tip counts at its own height):
 *
 *     root:                    x = log(h_root - a(root))
 *     other free node v:       x = logit(r_v),   r_v = (h_v - a(v)) / (h_parent - a(v))
 *     fake node (parent of a sampled ancestor): pinned to the sampled ancestor's height; x = 0 is a dummy
 *
 * This is the ratio parameterisation of NodeHeightToRatiosFullTransformDelegate composed with log/logit, except
 * that fake nodes are held fixed rather than given ratio 0 (where the logit is infinite). Fake coordinates have
 * zero gradient, so their momentum is constant and their kinetic energy cancels; the inverse ignores them.
 * Which nodes are fake is read from the tree on every call, so the transform follows sampled-ancestor moves
 * made by other operators between trajectories; along a trajectory the set cannot change, because every free
 * node stays strictly above all tips below it (the inverse throws ArithmeticException otherwise, which the HMC
 * operator treats as a numerical instability and rejects).
 *
 * logJacobian is log|dx/dh| (forward), as Transform expects; updateGradientLogDensity returns the gradient
 * w.r.t. x of log pi(h(x)) + log|dh/dx|.
 *
 * @author Giovanni Mazzeo
 */
public class SampledAncestorNodeHeightTransform extends Transform.MultivariateTransform implements Reportable {

    private final TreeModel tree;
    private final int external;
    private final int internal;

    // per internal index, recomputed on every call
    private final double[] anchor;
    private final boolean[] fake;
    private final int[] parentIndex;     // -1 for the root
    private final int[] preOrder;        // internal indices, parents first

    public SampledAncestorNodeHeightTransform(TreeModel tree) {
        super(tree.getInternalNodeCount());
        this.tree = tree;
        this.external = tree.getExternalNodeCount();
        this.internal = tree.getInternalNodeCount();
        this.anchor = new double[internal];
        this.fake = new boolean[internal];
        this.parentIndex = new int[internal];
        this.preOrder = new int[internal];
    }

    public TreeModel getTree() {
        return tree;
    }

    private int index(NodeRef node) {
        return node.getNumber() - external;
    }

    /** Topology, anchors and fake nodes from the current tree. */
    private void refresh() {
        int[] count = {0};
        visit(tree.getRoot(), -1, count);
        assert count[0] == internal;
    }

    /** Returns the largest tip height below node. */
    private double visit(NodeRef node, int parent, int[] count) {
        if (tree.isExternal(node)) {
            return tree.getNodeHeight(node);
        }
        int i = index(node);
        parentIndex[i] = parent;
        preOrder[count[0]++] = i;
        fake[i] = SampledAncestorTrees.isFake(tree, node);
        double a = Double.NEGATIVE_INFINITY;
        for (int c = 0; c < tree.getChildCount(node); c++) {
            a = Math.max(a, visit(tree.getChild(node, c), i, count));
        }
        anchor[i] = a;
        return a;
    }

    /** The sampled-ancestor tip below a fake node. */
    private double sampledAncestorHeight(int i) {
        NodeRef node = tree.getNode(i + external);
        for (int c = 0; c < tree.getChildCount(node); c++) {
            NodeRef child = tree.getChild(node, c);
            if (SampledAncestorTrees.isSampledAncestor(tree, child)) {
                return tree.getNodeHeight(child);
            }
        }
        throw new IllegalStateException("node " + (i + external) + " is not fake");
    }

    private static double logistic(double x) {
        return x >= 0 ? 1.0 / (1.0 + Math.exp(-x)) : Math.exp(x) / (1.0 + Math.exp(x));
    }

    @Override
    protected double[] transform(double[] h) {
        refresh();
        double[] x = new double[internal];
        for (int i = 0; i < internal; i++) {
            if (fake[i]) {
                x[i] = 0.0;
            } else if (parentIndex[i] < 0) {
                x[i] = Math.log(h[i] - anchor[i]);
            } else {
                x[i] = Math.log(h[i] - anchor[i]) - Math.log(h[parentIndex[i]] - h[i]);
            }
        }
        return x;
    }

    @Override
    protected double[] inverse(double[] x) {
        refresh();
        double[] h = new double[internal];
        for (int k = 0; k < internal; k++) {
            int i = preOrder[k];
            int p = parentIndex[i];
            if (fake[i]) {
                h[i] = sampledAncestorHeight(i);
            } else if (p < 0) {
                h[i] = anchor[i] + Math.exp(x[i]);
            } else {
                h[i] = anchor[i] + logistic(x[i]) * (h[p] - anchor[i]);
            }
            if (!fake[i] && !(h[i] > anchor[i] && (p < 0 || h[i] < h[p]) && h[i] < Double.POSITIVE_INFINITY)) {
                throw new ArithmeticException("node " + (i + external) + " left its interval");
            }
        }
        return h;
    }

    /** log|dh/dx|, the log Jacobian of the inverse, at heights h. */
    private double logJacobianInverse(double[] h) {
        double sum = 0.0;
        for (int i = 0; i < internal; i++) {
            if (fake[i]) {
                continue;
            }
            int p = parentIndex[i];
            if (p < 0) {
                sum += Math.log(h[i] - anchor[i]);
            } else {
                sum += Math.log(h[i] - anchor[i]) + Math.log(h[p] - h[i]) - Math.log(h[p] - anchor[i]);
            }
        }
        return sum;
    }

    @Override
    protected double getLogJacobian(double[] h) {
        refresh();
        return -logJacobianInverse(h);
    }

    /**
     * Reverse-mode chain rule. heightGradient is d/dh of the target; with withJacobian the gradient of
     * log|dh/dx| is added.
     */
    private double[] gradientWrtX(double[] heightGradient, double[] h, boolean withJacobian) {
        refresh();
        double[] g = heightGradient.clone();
        if (withJacobian) {
            // the parts of log|dh/dx| not written directly in x:  sum_v log(h_parent(v) - a(v))
            for (int i = 0; i < internal; i++) {
                int p = parentIndex[i];
                if (!fake[i] && p >= 0) {
                    g[p] += 1.0 / (h[p] - anchor[i]);
                }
            }
        }
        // post-order: total derivative w.r.t. h_i with the descendants' x held fixed
        for (int k = internal - 1; k >= 0; k--) {
            int i = preOrder[k];
            int p = parentIndex[i];
            if (!fake[i] && p >= 0) {
                double r = (h[i] - anchor[i]) / (h[p] - anchor[i]);
                g[p] += r * g[i];
            }
        }
        double[] gx = new double[internal];
        for (int i = 0; i < internal; i++) {
            if (fake[i]) {
                continue;
            }
            int p = parentIndex[i];
            if (p < 0) {
                gx[i] = g[i] * (h[i] - anchor[i]) + (withJacobian ? 1.0 : 0.0);
            } else {
                double r = (h[i] - anchor[i]) / (h[p] - anchor[i]);
                gx[i] = g[i] * r * (1.0 - r) * (h[p] - anchor[i]) + (withJacobian ? 1.0 - 2.0 * r : 0.0);
            }
        }
        return gx;
    }

    @Override
    protected double[] updateGradientLogDensity(double[] gradient, double[] h) {
        return gradientWrtX(gradient, h, true);
    }

    @Override
    public double[] updateGradientUnWeightedLogDensity(double[] gradient, double[] h, int from, int to) {
        return gradientWrtX(gradient, h, false);
    }

    @Override
    protected double[] updateGradientInverseUnWeightedLogDensity(double[] gradient, double[] x) {
        return gradientWrtX(gradient, inverse(x), false);
    }

    @Override
    protected double[] getGradientLogJacobianInverse(double[] x) {
        return gradientWrtX(new double[internal], inverse(x), true);
    }

    @Override
    public double[][] computeJacobianMatrixInverse(double[] x) {
        throw new RuntimeException("Not implemented; use updateGradientLogDensity");
    }

    @Override
    protected boolean isInInteriorDomain(double[] h) {
        refresh();
        for (int i = 0; i < internal; i++) {
            if (fake[i]) {
                continue;
            }
            int p = parentIndex[i];
            if (!(h[i] > anchor[i]) || (p >= 0 && !(h[i] < h[p]))) {
                return false;
            }
        }
        return true;
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public double[] gradient(double[] values, int from, int to) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public double[] gradientInverse(double[] values, int from, int to) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public String getTransformName() {
        return "saNodeHeightTransform";
    }

    @Override
    public String getReport() {
        refresh();
        int count = 0;
        for (boolean f : fake) {
            if (f) count++;
        }
        return "saNodeHeightTransform: " + internal + " internal nodes, " + count + " pinned (fake)";
    }
}
