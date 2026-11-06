/*
 * TMRCAPriorGradientTest.java
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

package test.dr.evomodel.tree;

import dr.app.beast.BeastParser;
import dr.evomodel.tree.TMRCAPriorGradient;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.FileReader;
import java.util.Map;

/**
 * TMRCAPriorGradient (examples/TestXML/testTMRCAPriorGradient.xml): the node-height gradient of clade-age priors
 * (normal, offset exponential stacked on the same clade, offset log-normal, uniform, root) against central
 * differences of the priors' joint log density, at three tree scalings.
 */
public class TMRCAPriorGradientTest extends TestCase {

    private static final String XML = "examples/TestXML/testTMRCAPriorGradient.xml";

    private Map<String, XMLObject> store;

    public void setUp() throws Exception {
        super.setUp();
        BeastParser parser = new BeastParser(new String[]{XML}, null, false, false, false, null);
        store = parser.parse(new FileReader(XML), true);
    }

    private Object get(String id) {
        return store.get(id).getNativeObject();
    }

    public void testGradientMatchesFiniteDifferences() {
        TMRCAPriorGradient grad = (TMRCAPriorGradient) get("grad");
        TreeModel tree = (TreeModel) get("treeModel");
        Parameter heights = grad.getParameter();
        Likelihood priors = grad.getLikelihood();
        int n = heights.getDimension();
        assertEquals(tree.getInternalNodeCount(), n);
        double[] start = heights.getParameterValues();

        for (double scale : new double[]{1.0, 0.92, 1.07}) {
            double[] h = new double[n];
            for (int i = 0; i < n; i++) h[i] = start[i] * scale;     // tips at 0: scaling keeps the tree valid
            set(heights, h);
            priors.makeDirty();
            double[] analytic = grad.getGradientLogDensity();
            double worst = 0.0;
            for (int j = 0; j < n; j++) {
                double step = 1e-5;
                double[] hp = h.clone();
                double[] hm = h.clone();
                hp[j] += step;
                hm[j] -= step;
                double numeric = (density(heights, priors, hp) - density(heights, priors, hm)) / (2 * step);
                double err = Math.abs(analytic[j] - numeric) / (1.0 + Math.abs(numeric));
                worst = Math.max(worst, err);
                assertTrue("scale " + scale + " node " + j + ": analytic " + analytic[j] + " numeric " + numeric,
                        err < 1e-6);
            }
            set(heights, h);
            System.out.println("scale " + scale + ": log prior " + priors.getLogLikelihood()
                    + ", worst relative gradient error " + worst);
        }
        set(heights, start);
    }

    /** The two priors stacked on one clade add, and a uniform contributes nothing inside its support. */
    public void testContributionsAdd() {
        TMRCAPriorGradient grad = (TMRCAPriorGradient) get("grad");
        TreeModel tree = (TreeModel) get("treeModel");
        double[] g = grad.getGradientLogDensity();
        double sum = 0.0;
        for (double x : g) sum += Math.abs(x);
        assertTrue(sum > 0.0);
        // node heights: AB 1.0, DE 0.5, CDE 1.5, root 3.0
        int ext = tree.getExternalNodeCount();
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            double height = tree.getNodeHeight(tree.getNode(ext + i));
            double expected;
            if (Math.abs(height - 1.0) < 1e-12) {
                expected = -(1.0 - 1.1) / (0.2 * 0.2) - 1.0 / 2.0;          // normal + exponential(mean 2)
            } else if (Math.abs(height - 1.5) < 1e-12) {
                expected = 0.0;                                              // uniform
            } else if (Math.abs(height - 3.0) < 1e-12) {
                expected = -(3.0 - 3.5) / (0.4 * 0.4);                       // normal on the root
            } else {
                double x = height - 0.2;                                     // log-normal, offset 0.2
                expected = -(1.0 + (Math.log(x) + 1.2) / 0.25) / x;
            }
            assertEquals("height " + height, expected, g[i], 1e-6 * (1.0 + Math.abs(expected)));
        }
    }

    private static void set(Parameter heights, double[] h) {
        for (int i = 0; i < h.length; i++) heights.setParameterValueQuietly(i, h[i]);
        heights.fireParameterChangedEvent();
    }

    private static double density(Parameter heights, Likelihood priors, double[] h) {
        set(heights, h);
        priors.makeDirty();
        return priors.getLogLikelihood();
    }
}
