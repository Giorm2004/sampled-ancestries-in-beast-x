/*
 * AscertainedTreeDataLikelihoodTest.java
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
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.FileReader;
import java.util.Map;

/**
 * Ascertainment correction in TreeDataLikelihood (examples/TestXML/testAscertainedTreeDataLikelihood.xml):
 * the corrected log likelihood equals BeagleTreeLikelihood's, and the pre-order branch-rate and
 * node-height gradients and diagonal Hessians equal central differences of the corrected log
 * likelihood. Two partitions (one with gamma rate categories), dated tips, missing data.
 */
public class AscertainedTreeDataLikelihoodTest extends TestCase {

    private static final String XML = "examples/TestXML/testAscertainedTreeDataLikelihood.xml";

    private Map<String, XMLObject> store;

    public void setUp() throws Exception {
        super.setUp();
        BeastParser parser = new BeastParser(new String[]{XML}, null, false, false, false, null);
        store = parser.parse(new FileReader(XML), true);
    }

    private Object get(String id) {
        return store.get(id).getNativeObject();
    }

    public void testLikelihoodMatchesBeagleTreeLikelihood() {
        for (int p = 1; p <= 2; p++) {
            double tdl = ((Likelihood) get("tdl" + p)).getLogLikelihood();
            double old = ((Likelihood) get("old" + p)).getLogLikelihood();
            System.out.println("partition " + p + ": TreeDataLikelihood " + tdl + "  BeagleTreeLikelihood " + old);
            assertEquals("partition " + p, old, tdl, 1e-9 * Math.abs(old));
        }
    }

    public void testBranchRateDerivatives() {
        for (int p = 1; p <= 2; p++) {
            checkDerivatives((GradientWrtParameterProvider) get("rateGradient" + p), "rates, partition " + p);
        }
    }

    public void testNodeHeightDerivatives() {
        for (int p = 1; p <= 2; p++) {
            checkDerivatives((GradientWrtParameterProvider) get("heightGradient" + p), "heights, partition " + p);
        }
    }

    private static void checkDerivatives(GradientWrtParameterProvider provider, String label) {
        Parameter x = provider.getParameter();
        Likelihood likelihood = provider.getLikelihood();
        double[] gradient = provider.getGradientLogDensity();
        double[] hessian = provider instanceof HessianWrtParameterProvider
                ? ((HessianWrtParameterProvider) provider).getDiagonalHessianLogDensity() : null;
        double f0 = likelihood.getLogLikelihood();
        double worstGradient = 0.0;
        double worstHessian = 0.0;
        for (int i = 0; i < x.getDimension(); i++) {
            double xi = x.getParameterValue(i);
            double h1 = 1e-6 * Math.max(1.0, Math.abs(xi));
            double fp = valueAt(x, i, xi + h1, likelihood);
            double fm = valueAt(x, i, xi - h1, likelihood);
            double numeric = (fp - fm) / (2 * h1);
            worstGradient = Math.max(worstGradient, check(gradient[i], numeric, 1e-5, label + " gradient[" + i + "]"));
            if (hessian != null) {
                double h2 = 1e-4 * Math.max(1.0, Math.abs(xi));
                double fp2 = valueAt(x, i, xi + h2, likelihood);
                double fm2 = valueAt(x, i, xi - h2, likelihood);
                double numeric2 = (fp2 - 2 * f0 + fm2) / (h2 * h2);
                worstHessian = Math.max(worstHessian, check(hessian[i], numeric2, 1e-3, label + " hessian[" + i + "]"));
            }
            valueAt(x, i, xi, likelihood);
        }
        System.out.println(label + ": " + x.getDimension() + " coordinates, worst gradient error " + worstGradient
                + (hessian != null ? ", worst diagonal-Hessian error " + worstHessian : ", no Hessian"));
    }

    private static double valueAt(Parameter x, int i, double value, Likelihood likelihood) {
        x.setParameterValue(i, value);
        likelihood.makeDirty();
        return likelihood.getLogLikelihood();
    }

    /** |analytic - numeric| / (1 + |numeric|), asserted below tol; returns it. */
    private static double check(double analytic, double numeric, double tol, String label) {
        double err = Math.abs(analytic - numeric) / (1.0 + Math.abs(numeric));
        assertTrue(label + ": analytic " + analytic + " vs numeric " + numeric, err < tol);
        return err;
    }
}
