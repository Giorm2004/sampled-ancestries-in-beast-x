/*
 * EpisodicBDSKYSampledAncestorTest.java
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

package test.dr.evomodel.birthdeath;

import dr.app.beast.BeastParser;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.FileReader;
import java.util.Map;

/**
 * The episodic birth-death sampling model with BDSKY's removal-aware epidemiological parameterisation and
 * sampled sampling-window break points (examples/TestXML/testEpisodicBDSKYSampledAncestors.xml).
 *
 * Reference values: BEAST 2 BDSKY 1.5.1, BirthDeathSkylineModel with reproductiveNumber, becomeUninfectiousRate,
 * samplingProportion "0 s 0", samplingRateChangeTimes "b2 b1 0" with reverseTimeArrays (false false true false
 * false), rho at present, removalProbability 0, conditionOnSurvival (default), on the same tree. BDSKY's
 * sampled-ancestor likelihood adds log(2) per internal node (13 here), a labelling constant BEAST X omits.
 */
public class EpisodicBDSKYSampledAncestorTest extends TestCase {

    private static final String XML = "examples/TestXML/testEpisodicBDSKYSampledAncestors.xml";
    private static final double[] BDSKY = {-75.64849004553956, -60.90356377782266, -90.40905857121427};
    private static final double LABELLING = 0.0;

    private Map<String, XMLObject> store;

    public void setUp() throws Exception {
        super.setUp();
        BeastParser parser = new BeastParser(new String[]{XML}, null, false, false, false, null);
        store = parser.parse(new FileReader(XML), true);
    }

    private Object get(String id) {
        return store.get(id).getNativeObject();
    }

    private double logL(int set) {
        Likelihood l = (Likelihood) get("bd" + set);
        l.makeDirty();
        return l.getLogLikelihood();
    }

    public void testMatchesBDSKY() {
        for (int k = 1; k <= 3; k++) {
            double value = logL(k);
            System.out.println("set " + k + ": BEAST X " + value + "  BDSKY " + BDSKY[k - 1]);
            assertEquals("set " + k, BDSKY[k - 1] - LABELLING, value, 1e-10);
        }
    }

    /** Moving set 1's parameters, break points included, onto set 3's values reproduces set 3. */
    public void testParametersAndBreakPointsAreLive() {
        String[][] ids = {{"R0", "delta", "s", "rho", "origin", "window"}};
        for (String id : ids[0]) {
            Parameter from = (Parameter) get(id + ".3");
            Parameter to = (Parameter) get(id + ".1");
            for (int i = 0; i < from.getDimension(); i++) {
                to.setParameterValue(i, from.getParameterValue(i));
            }
        }
        assertEquals(logL(3), logL(1), 1e-10);
    }

    public void testBreakPointsOutOfOrder() {
        Parameter window = (Parameter) get("window.1");
        window.setParameterValue(0, 5.0);          // above the second break point (4.401)
        assertEquals(Double.NEGATIVE_INFINITY, logL(1));
        window.setParameterValue(0, 0.02);
        assertEquals(BDSKY[0] - LABELLING, logL(1), 1e-10);
    }

    public void testNodeHeightGradient() {
        for (int k = 1; k <= 3; k++) {
            GradientWrtParameterProvider provider = (GradientWrtParameterProvider) get("grad" + k);
            Parameter x = provider.getParameter();
            Likelihood likelihood = provider.getLikelihood();
            double[] gradient = provider.getGradientLogDensity();
            double worst = 0.0;
            dr.evomodel.tree.TreeModel tree = (dr.evomodel.tree.TreeModel) get("treeModel");
            int fake = 0;
            for (int i = 0; i < x.getDimension(); i++) {
                if (dr.evomodel.birthdeath.EpisodicBirthDeathSamplingModel.isFakeNode(tree,
                        tree.getNode(i + tree.getExternalNodeCount()))) {
                    assertEquals("fake node " + i + " is pinned", 0.0, gradient[i]);
                    fake++;
                    continue;
                }
                double xi = x.getParameterValue(i);
                double h = 1e-6 * Math.max(1.0, Math.abs(xi));
                x.setParameterValue(i, xi + h);
                likelihood.makeDirty();
                double fp = likelihood.getLogLikelihood();
                x.setParameterValue(i, xi - h);
                likelihood.makeDirty();
                double fm = likelihood.getLogLikelihood();
                x.setParameterValue(i, xi);
                double numeric = (fp - fm) / (2 * h);
                double err = Math.abs(gradient[i] - numeric) / (1.0 + Math.abs(numeric));
                assertTrue("set " + k + " node " + i + ": analytic " + gradient[i] + " vs numeric " + numeric, err < 1e-5);
                worst = Math.max(worst, err);
            }
            assertEquals("fake nodes", 2, fake);
            System.out.println("set " + k + ": node-height gradient, " + (x.getDimension() - fake) + " free nodes ("
                    + fake + " pinned), worst error " + worst);
        }
    }
}
