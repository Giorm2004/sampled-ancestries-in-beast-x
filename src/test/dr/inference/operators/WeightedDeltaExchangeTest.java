/*
 * WeightedDeltaExchangeTest.java
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

package test.dr.inference.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.DeltaExchangeOperator;
import dr.math.MathUtils;
import junit.framework.TestCase;

/**
 * The weighted delta exchange must sample the target restricted to the plane sum_k w_k x_k = const.
 * Target: independent Gaussians whose means lie on the plane, so the exact restricted mean is the vector of means.
 * Metropolis-Hastings with the operator alone, at two window sizes: the old move (x1 -= d, x2 += d w1/w2 with a
 * zero Hastings ratio) misses the mean by up to ~0.6 sd here, by a delta-dependent amount; the symmetric move must
 * hit it to Monte Carlo error.
 */
public class WeightedDeltaExchangeTest extends TestCase {

    // weights, means and sds of the IE-CoR per-class rate multipliers (posterior-like scales)
    private static final int[] W = {1497, 1293, 839, 522, 224, 219, 131, 95};
    private static final double[] M0 = {1.035, 0.7386, 1.255, 0.428, 2.042, 1.493, 1.863, 0.1034};
    private static final double[] S = {0.032, 0.027, 0.059, 0.023, 0.21, 0.15, 0.27, 0.023};

    private static double[] means() {
        double wsum = 0.0, wm = 0.0;
        for (int k = 0; k < W.length; k++) {
            wsum += W[k];
            wm += W[k] * M0[k];
        }
        double[] m = new double[W.length];
        for (int k = 0; k < W.length; k++) {
            m[k] = M0[k] * wsum / wm;           // on the plane sum w x = sum w
        }
        return m;
    }

    private static double logTarget(double[] x, double[] m) {
        double lp = 0.0;
        for (int k = 0; k < x.length; k++) {
            if (x[k] <= 0.0) return Double.NEGATIVE_INFINITY;
            double z = (x[k] - m[k]) / S[k];
            lp -= 0.5 * z * z;
        }
        return lp;
    }

    private double[] sampleMean(double delta, long seed, int n) {
        double[] m = means();
        Parameter x = new Parameter.Default("x", m.clone());
        x.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, m.length));
        DeltaExchangeOperator op = new DeltaExchangeOperator(x, W, delta, 1.0, false,
                AdaptationMode.ADAPTATION_OFF);
        MathUtils.setSeed(seed);
        double lp = logTarget(x.getParameterValues(), m);
        double[] sum = new double[m.length];
        int burn = n / 10;
        for (int t = 0; t < n; t++) {
            double[] old = x.getParameterValues();
            double hr = op.doOperation();
            double lpNew = hr == Double.NEGATIVE_INFINITY ? Double.NEGATIVE_INFINITY : logTarget(x.getParameterValues(), m);
            if (Math.log(MathUtils.nextDouble()) < lpNew - lp + hr) {
                lp = lpNew;
            } else {
                for (int k = 0; k < m.length; k++) x.setParameterValueQuietly(k, old[k]);
            }
            if (t >= burn) {
                for (int k = 0; k < m.length; k++) sum[k] += x.getParameterValue(k);
            }
        }
        for (int k = 0; k < m.length; k++) sum[k] /= (n - burn);
        return sum;
    }

    public void testPreservesTheWeightedSum() {
        double[] m = means();
        Parameter x = new Parameter.Default("x", m.clone());
        x.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, m.length));
        DeltaExchangeOperator op = new DeltaExchangeOperator(x, W, 0.3, 1.0, false,
                AdaptationMode.ADAPTATION_OFF);
        MathUtils.setSeed(3);
        double before = 0.0;
        for (int k = 0; k < W.length; k++) before += W[k] * x.getParameterValue(k);
        for (int t = 0; t < 1000; t++) op.doOperation();
        double after = 0.0;
        for (int k = 0; k < W.length; k++) after += W[k] * x.getParameterValue(k);
        assertEquals(before, after, 1e-9 * before);
    }

    public void testSamplesTheRestrictedTarget() {
        double[] m = means();
        for (double delta : new double[]{0.23, 0.37}) {
            double[] mean = sampleMean(delta, 11, 3_000_000);
            double worst = 0.0;
            for (int k = 0; k < m.length; k++) {
                double z = (mean[k] - m[k]) / S[k];
                worst = Math.max(worst, Math.abs(z));
                assertTrue("delta " + delta + " component " + k + ": mean " + mean[k] + " vs exact " + m[k]
                        + " (" + z + " sd)", Math.abs(z) < 0.05);
            }
            System.out.println("delta " + delta + ": worst |mean - exact| = " + worst + " sd");
        }
    }

    public void testEqualWeightsUnchanged() {
        Parameter x = new Parameter.Default("x", new double[]{0.2, 0.3, 0.5});
        x.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 3));
        DeltaExchangeOperator op = new DeltaExchangeOperator(x, new int[]{4, 4, 4}, 0.1, 1.0, false,
                AdaptationMode.ADAPTATION_OFF);
        MathUtils.setSeed(7);
        op.doOperation();
        // the same draws as the original equal-weight move: dim1, dim2, then d
        MathUtils.setSeed(7);
        int dim1 = MathUtils.nextInt(3);
        int dim2 = dim1;
        while (dim1 == dim2) dim2 = MathUtils.nextInt(3);
        double d = MathUtils.nextDouble() * 0.1;
        double[] expect = {0.2, 0.3, 0.5};
        expect[dim1] -= d;
        expect[dim2] += d;
        for (int k = 0; k < 3; k++) assertEquals(expect[k], x.getParameterValue(k), 0.0);
    }
}
