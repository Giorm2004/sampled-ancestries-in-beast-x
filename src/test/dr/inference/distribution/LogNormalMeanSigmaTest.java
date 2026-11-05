/*
 * LogNormalMeanSigmaTest.java
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

package test.dr.inference.distribution;

import dr.inference.distribution.LogNormalDistributionModel;
import dr.inference.distribution.LogNormalDistributionModel.Parameterization;
import dr.inference.model.Parameter;
import dr.inferencexml.distribution.LogNormalDistributionModelParser;
import dr.xml.XMLObject;
import dr.xml.XMLParser;
import junit.framework.TestCase;

import java.io.StringReader;
import java.util.Map;

/**
 * The MEAN_SIGMA parameterisation (real-space mean, log-space sigma: mu = log(mean) - sigma^2 / 2), BEAST 2's
 * meanInRealSpace convention, selected in XML by meanInRealSpace="true" stdevInRealSpace="false".
 */
public class LogNormalMeanSigmaTest extends TestCase {

    private static final double MEAN = 0.08, SIGMA = 0.3;

    private static LogNormalDistributionModel reference() {
        return new LogNormalDistributionModel(Parameterization.MU_SIGMA,
                new Parameter.Default(Math.log(MEAN) - 0.5 * SIGMA * SIGMA), new Parameter.Default(SIGMA), 0.0);
    }

    private static void assertSameDistribution(LogNormalDistributionModel model) {
        LogNormalDistributionModel ref = reference();
        assertEquals(ref.getMu(), model.getMu(), 1e-15);
        assertEquals(SIGMA, model.getSigma(), 0.0);
        assertEquals(MEAN, model.getMean(), 1e-15);
        assertEquals(ref.getStdev(), model.getStdev(), 1e-15);
        for (double p : new double[]{0.0047, 0.3, 0.5, 0.9}) {
            assertEquals(ref.quantile(p), model.quantile(p), 1e-15);
        }
        for (double x : new double[]{0.02, 0.08, 0.2}) {
            assertEquals(ref.logPdf(x), model.logPdf(x), 1e-12);
        }
        // the real-space mean is the mean
        assertEquals(MEAN, Math.exp(model.getMu() + 0.5 * SIGMA * SIGMA), 1e-15);
    }

    public void testConstructor() {
        assertSameDistribution(new LogNormalDistributionModel(Parameterization.MEAN_SIGMA,
                new Parameter.Default(MEAN), new Parameter.Default(SIGMA), 0.0));
    }

    public void testParser() throws Exception {
        String xml = "<beast><logNormalDistributionModel id=\"ln\" meanInRealSpace=\"true\" stdevInRealSpace=\"false\">"
                + "<mean><parameter value=\"" + MEAN + "\"/></mean><stdev><parameter value=\"" + SIGMA + "\"/></stdev>"
                + "</logNormalDistributionModel></beast>";
        XMLParser parser = new XMLParser(false, false, false, null);
        parser.addXMLObjectParser(new LogNormalDistributionModelParser());
        parser.addXMLObjectParser(new dr.inference.model.ParameterParser());
        Map<String, XMLObject> store = parser.parse(new StringReader(xml), true);
        LogNormalDistributionModel model = (LogNormalDistributionModel) store.get("ln").getNativeObject();
        assertEquals(Parameterization.MEAN_SIGMA, model.getParameterization());
        assertSameDistribution(model);
    }
}
