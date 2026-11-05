/*
 * ProductBranchRateModelTest.java
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

package test.dr.evomodel.branchratemodel;

import dr.app.beast.BeastParser;
import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.ProductBranchRateModel;
import dr.evomodel.tree.SampledAncestorTrees;
import dr.evomodel.tree.TreeModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.FileReader;
import java.util.Map;

/**
 * ProductBranchRateModel (relaxed x time-varying) in a binary-covarion TreeDataLikelihood on a tree
 * with two sampled ancestors (examples/TestXML/testProductBranchRates.xml): the likelihood equals
 * BeagleTreeLikelihood, the rates of the zero-length sampled-ancestor branches are finite, and the gradients w.r.t.
 * the relaxed rates, the time-varying log-rates and the node heights match central differences. The node-height
 * gradient needs the instantaneous rate at each node (TimeDependentBranchRates), not the branch average.
 */
public class ProductBranchRateModelTest extends TestCase {

    private static final String XML = "examples/TestXML/testProductBranchRates.xml";

    private Map<String, XMLObject> store;
    private TreeModel tree;
    private Likelihood tdl;

    public void setUp() throws Exception {
        super.setUp();
        BeastParser parser = new BeastParser(new String[]{XML}, null, false, false, false, null);
        store = parser.parse(new FileReader(XML), true);
        tree = (TreeModel) get("treeModel");
        tdl = (Likelihood) get("tdl");
    }

    private Object get(String id) {
        return store.get(id).getNativeObject();
    }

    private double logL() {
        tdl.makeDirty();
        return tdl.getLogLikelihood();
    }

    public void testLikelihoodMatchesBeagleTreeLikelihood() {
        Likelihood reference = (Likelihood) get("reference");
        reference.makeDirty();
        double a = logL(), b = reference.getLogLikelihood();
        System.out.println("TreeDataLikelihood " + a + "  BeagleTreeLikelihood " + b);
        assertFalse(Double.isNaN(a));
        assertEquals(b, a, 1e-9 * Math.abs(b));
    }

    /** The zero-length branches above a1 (0.9) and a4 (3.0) get relaxed x exp(log-rate of their epoch). */
    public void testSampledAncestorBranchRates() {
        ProductBranchRateModel clock = (ProductBranchRateModel) get("clock");
        Parameter logRates = (Parameter) get("timeVarying.logRates");
        double[] grid = ((Parameter) get("timeVarying.grid")).getParameterValues();
        int seen = 0;
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef tip = tree.getExternalNode(i);
            if (!SampledAncestorTrees.isSampledAncestor(tree, tip)) continue;
            double h = tree.getNodeHeight(tip);
            int epoch = 0;
            while (epoch < grid.length && grid[epoch] <= h) epoch++;
            double expected = clock.getFactors().get(0).getBranchRate(tree, tip) * Math.exp(logRates.getParameterValue(epoch));
            double rate = clock.getBranchRate(tree, tip);
            System.out.println(tree.getTaxonId(i) + " at " + h + ": rate " + rate + " (expected " + expected + ")");
            assertEquals(expected, rate, 1e-12);
            seen++;
        }
        assertEquals(2, seen);
    }

    public void testRelaxedRateGradient() {
        checkGradient("relaxedGradient", false);
    }

    public void testTimeVaryingLogRateGradient() {
        checkGradient("timeVaryingGradient", false);
    }

    public void testNodeHeightGradient() {
        checkGradient("heightGradient", true);
    }

    private void checkGradient(String id, boolean heights) {
        GradientWrtParameterProvider provider = (GradientWrtParameterProvider) get(id);
        Parameter x = provider.getParameter();
        logL();
        double[] analytic = provider.getGradientLogDensity();
        assertEquals(x.getDimension(), analytic.length);
        double worst = 0.0;
        int checked = 0;
        for (int i = 0; i < x.getDimension(); i++) {
            if (heights && SampledAncestorTrees.isFake(tree, tree.getNode(tree.getExternalNodeCount() + i))) {
                continue;   // pinned to its sampled ancestor: not a free coordinate
            }
            double xi = x.getParameterValue(i);
            double h = 1e-5 * Math.max(1.0, Math.abs(xi));
            x.setParameterValue(i, xi + h);
            double fp = logL();
            x.setParameterValue(i, xi - h);
            double fm = logL();
            x.setParameterValue(i, xi);
            double numeric = (fp - fm) / (2 * h);
            double err = Math.abs(analytic[i] - numeric) / (1.0 + Math.abs(numeric));
            assertTrue(id + " [" + i + "]: analytic " + analytic[i] + " vs numeric " + numeric, err < 1e-5);
            worst = Math.max(worst, err);
            checked++;
        }
        logL();
        System.out.println(id + ": " + checked + " coordinates, worst relative error " + worst);
    }
}
