/*
 * SampledAncestorNodeHeightTransformTest.java
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

package test.dr.evomodel.treedatalikelihood;

import dr.app.beast.BeastParser;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.SampledAncestorTrees;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.discrete.SampledAncestorNodeHeightTransform;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.GeneralOperator;
import dr.inference.operators.MCMCOperator;
import dr.math.MathUtils;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.FileReader;
import java.util.Map;

/**
 * SampledAncestorNodeHeightTransform on a tree with two sampled ancestors
 * (examples/TestXML/testSampledAncestorNodeHeightTransform.xml): round trip, pinned fake nodes, log Jacobian and
 * gradient against finite differences, and HMC trajectories that keep the sampled ancestors in place.
 */
public class SampledAncestorNodeHeightTransformTest extends TestCase {

    private static final String XML = "examples/TestXML/testSampledAncestorNodeHeightTransform.xml";

    private Map<String, XMLObject> store;
    private TreeModel tree;
    private Parameter heights;
    private SampledAncestorNodeHeightTransform transform;
    private boolean[] fake;
    private int n;

    public void setUp() throws Exception {
        super.setUp();
        BeastParser parser = new BeastParser(new String[]{XML}, null, false, false, false, null);
        store = parser.parse(new FileReader(XML), true);
        tree = (TreeModel) get("treeModel");
        heights = (Parameter) get("heights");
        transform = (SampledAncestorNodeHeightTransform) get("saTransform");
        n = heights.getDimension();
        fake = new boolean[n];
        for (int i = 0; i < n; i++) {
            fake[i] = SampledAncestorTrees.isFake(tree, tree.getNode(tree.getExternalNodeCount() + i));
        }
    }

    private Object get(String id) {
        return store.get(id).getNativeObject();
    }

    private double[] forward(double[] h) {
        return transform.transform(h, 0, h.length);
    }

    private double[] inverse(double[] x) {
        return transform.inverse(x, 0, x.length);
    }

    private void setHeights(double[] h) {
        for (int i = 0; i < n; i++) {
            heights.setParameterValueQuietly(i, h[i]);
        }
        heights.fireParameterChangedEvent();
    }

    public void testPinsTheParentsOfSampledAncestors() {
        int count = 0;
        for (boolean f : fake) {
            if (f) count++;
        }
        assertEquals(2, count);
        assertEquals(2, SampledAncestorTrees.countSampledAncestors(tree));
        assertEquals(n, tree.getInternalNodeCount());
    }

    public void testRoundTrip() {
        double[] h = heights.getParameterValues();
        double[] x = forward(h);
        double[] back = inverse(x);
        for (int i = 0; i < n; i++) {
            assertEquals("node " + i, h[i], back[i], 1e-12 * Math.max(1.0, h[i]));
            if (fake[i]) {
                assertEquals(0.0, x[i]);
            }
        }
        // the fake coordinates are dummies
        for (int i = 0; i < n; i++) {
            if (fake[i]) x[i] = 3.7;
        }
        double[] again = inverse(x);
        for (int i = 0; i < n; i++) {
            assertEquals(back[i], again[i], 0.0);
        }
    }

    /** -logJacobian(h) = log|det dh/dx| over the free coordinates, by finite differences of the inverse. */
    public void testLogJacobian() {
        double[] h = heights.getParameterValues();
        double[] x = forward(h);
        int[] free = freeIndices();
        int m = free.length;
        double[][] jac = new double[m][m];
        double step = 1e-6;
        for (int b = 0; b < m; b++) {
            double[] xp = x.clone();
            double[] xm = x.clone();
            xp[free[b]] += step;
            xm[free[b]] -= step;
            double[] hp = inverse(xp);
            double[] hm = inverse(xm);
            for (int a = 0; a < m; a++) {
                jac[a][b] = (hp[free[a]] - hm[free[a]]) / (2 * step);
            }
        }
        double numeric = logAbsDet(jac);
        double analytic = -transform.logJacobian(h, 0, n);
        System.out.println("log|dh/dx|: analytic " + analytic + "  numeric " + numeric);
        assertEquals(numeric, analytic, 1e-6);
    }

    /** d/dx [log pi(h(x)) + log|dh/dx|] against central differences, at several points. */
    public void testGradient() {
        GradientWrtParameterProvider grad = (GradientWrtParameterProvider) get("grad1");
        Likelihood likelihood = grad.getLikelihood();
        assertSame(heights.getDimension(), grad.getParameter().getDimension());

        double[] start = heights.getParameterValues();
        double[] x0 = forward(start);
        MathUtils.setSeed(17);
        for (int trial = 0; trial < 4; trial++) {
            double[] x = x0.clone();
            if (trial > 0) {
                for (int i = 0; i < n; i++) {
                    // the root stays put: it must stay below the origin
                    if (!fake[i] && !tree.isRoot(tree.getNode(tree.getExternalNodeCount() + i))) {
                        x[i] += 0.5 * MathUtils.nextGaussian();
                    }
                }
            }
            double[] h = inverse(x);
            setHeights(h);
            likelihood.makeDirty();
            double[] analytic = transform.updateGradientLogDensity(grad.getGradientLogDensity(), h, 0, n);

            double worst = 0.0;
            for (int j = 0; j < n; j++) {
                double step = 1e-5;
                double[] xp = x.clone();
                double[] xm = x.clone();
                xp[j] += step;
                xm[j] -= step;
                double numeric = (target(likelihood, xp) - target(likelihood, xm)) / (2 * step);
                double err = Math.abs(analytic[j] - numeric) / (1.0 + Math.abs(numeric));
                assertTrue("trial " + trial + " coordinate " + j + (fake[j] ? " (fake)" : "") + ": analytic "
                        + analytic[j] + " vs numeric " + numeric, err < 1e-6);
                if (fake[j]) {
                    assertEquals(0.0, analytic[j]);
                }
                worst = Math.max(worst, err);
            }
            System.out.println("trial " + trial + ": gradient w.r.t. x, worst relative error " + worst);
        }
        setHeights(start);
    }

    private double target(Likelihood likelihood, double[] x) {
        double[] h = inverse(x);
        setHeights(h);
        likelihood.makeDirty();
        return likelihood.getLogLikelihood() - transform.logJacobian(h, 0, n);
    }

    public void testInverseRejectsADegenerateNode() {
        double[] x = forward(heights.getParameterValues());
        int j = -1;
        for (int i = 0; i < n; i++) {
            if (!fake[i] && !tree.isRoot(tree.getNode(tree.getExternalNodeCount() + i))) {
                j = i;
                break;
            }
        }
        x[j] = 800.0;     // ratio rounds to 1: the node would meet its parent
        try {
            inverse(x);
            fail("expected ArithmeticException");
        } catch (ArithmeticException expected) {
            // ok
        }
        x[j] = -800.0;    // ratio rounds to 0: the node would meet its highest tip
        try {
            inverse(x);
            fail("expected ArithmeticException");
        } catch (ArithmeticException expected) {
            // ok
        }
    }

    /** HMC trajectories, kept whether accepted or not, never move a fake node or create or remove a sampled ancestor. */
    public void testHmcKeepsSampledAncestors() {
        MCMCOperator hmc = (MCMCOperator) get("hmc");
        Likelihood likelihood = (Likelihood) get("bd1");
        double[] saHeights = sampledAncestorHeights();
        double[] before = heights.getParameterValues();
        boolean moved = false;
        MathUtils.setSeed(5);
        for (int k = 0; k < 200; k++) {
            ((GeneralOperator) hmc).operate(likelihood);
            hmc.accept(0.0);
            assertEquals("step " + k, 2, SampledAncestorTrees.countSampledAncestors(tree));
            for (int i = 0; i < n; i++) {
                NodeRef node = tree.getNode(tree.getExternalNodeCount() + i);
                if (fake[i]) {
                    assertTrue("step " + k + " node " + i, SampledAncestorTrees.isFake(tree, node));
                    assertEquals(saHeights[i], tree.getNodeHeight(node), 0.0);
                }
            }
            likelihood.makeDirty();
            assertFalse("step " + k + " root " + tree.getNodeHeight(tree.getRoot()) + " heights " + java.util.Arrays.toString(heights.getParameterValues()), Double.isNaN(likelihood.getLogLikelihood()));
        }
        double[] after = heights.getParameterValues();
        for (int i = 0; i < n; i++) {
            if (!fake[i] && after[i] != before[i]) moved = true;
        }
        assertTrue("HMC never moved a free node", moved);
    }

    private double[] sampledAncestorHeights() {
        double[] out = new double[n];
        for (int i = 0; i < n; i++) {
            out[i] = fake[i] ? tree.getNodeHeight(tree.getNode(tree.getExternalNodeCount() + i)) : Double.NaN;
        }
        return out;
    }

    private int[] freeIndices() {
        int m = 0;
        for (boolean f : fake) {
            if (!f) m++;
        }
        int[] free = new int[m];
        int k = 0;
        for (int i = 0; i < n; i++) {
            if (!fake[i]) free[k++] = i;
        }
        return free;
    }

    private static double logAbsDet(double[][] matrix) {
        int m = matrix.length;
        double[][] a = new double[m][];
        for (int i = 0; i < m; i++) a[i] = matrix[i].clone();
        double logDet = 0.0;
        for (int c = 0; c < m; c++) {
            int pivot = c;
            for (int r = c + 1; r < m; r++) {
                if (Math.abs(a[r][c]) > Math.abs(a[pivot][c])) pivot = r;
            }
            double[] t = a[c]; a[c] = a[pivot]; a[pivot] = t;
            logDet += Math.log(Math.abs(a[c][c]));
            for (int r = c + 1; r < m; r++) {
                double f = a[r][c] / a[c][c];
                for (int k = c; k < m; k++) a[r][k] -= f * a[c][k];
            }
        }
        return logDet;
    }
}
