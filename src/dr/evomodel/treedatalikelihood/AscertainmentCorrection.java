/*
 * AscertainmentCorrection.java
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

package dr.evomodel.treedatalikelihood;

import dr.evolution.alignment.AscertainedSitePatterns;

/**
 * Ascertainment correction for a TreeDataLikelihood, as in BeagleTreeLikelihood.
 *
 * The corrected log likelihood is sum_k w_k (l_k - log R), where l_k is the log likelihood of
 * pattern k and R is the probability of the ascertainment event,
 *
 *     R = I - E,   I = sum over included patterns of p_k (1 if there are none),
 *                  E = sum over excluded patterns of p_k,
 *
 * with p_k = exp(l_k); ascertainment patterns themselves carry weight 0. Writing a_k = +1 for an
 * included pattern, -1 for an excluded one and 0 otherwise, and c_k = a_k p_k / R, any derivative
 * of the correction follows from the per-pattern derivatives d_k = dl_k/dt and s_k = (d^2 L_k/dt^2) / L_k:
 *
 *     d/dt    [-W log R] = -W sum_k c_k d_k,
 *     d^2/dt^2[-W log R] = -W (sum_k c_k s_k - (sum_k c_k d_k)^2),     W = sum_k w_k,
 *
 * so the corrected gradient is sum_k (w_k - W c_k) d_k: an ordinary weighted sum with adjusted
 * weights.
 *
 * @author Giovanni Mazzeo
 */
public final class AscertainmentCorrection {

    private AscertainmentCorrection() { }

    /** log R; NaN when R is not positive. */
    public static double logProbability(AscertainedSitePatterns patterns, double[] patternLogLikelihoods) {
        double r = probability(patterns, patternLogLikelihoods);
        return r > 0.0 ? Math.log(r) : Double.NaN;
    }

    /** c_k = a_k p_k / R for every pattern (0 for patterns outside the ascertainment set). */
    public static double[] coefficients(AscertainedSitePatterns patterns, double[] patternLogLikelihoods) {
        double r = probability(patterns, patternLogLikelihoods);
        double[] c = new double[patternLogLikelihoods.length];
        int[] include = patterns.getIncludePatternIndices();
        for (int i = 0; i < patterns.getIncludePatternCount(); i++) {
            c[include[i]] += Math.exp(patternLogLikelihoods[include[i]]) / r;
        }
        int[] exclude = patterns.getExcludePatternIndices();
        for (int i = 0; i < patterns.getExcludePatternCount(); i++) {
            c[exclude[i]] -= Math.exp(patternLogLikelihoods[exclude[i]]) / r;
        }
        return c;
    }

    public static double totalWeight(double[] patternWeights) {
        double w = 0.0;
        for (double x : patternWeights) {
            w += x;
        }
        return w;
    }

    private static double probability(AscertainedSitePatterns patterns, double[] patternLogLikelihoods) {
        double r = patterns.getIncludePatternCount() > 0 ? 0.0 : 1.0;
        int[] include = patterns.getIncludePatternIndices();
        for (int i = 0; i < patterns.getIncludePatternCount(); i++) {
            r += Math.exp(patternLogLikelihoods[include[i]]);
        }
        int[] exclude = patterns.getExcludePatternIndices();
        for (int i = 0; i < patterns.getExcludePatternCount(); i++) {
            r -= Math.exp(patternLogLikelihoods[exclude[i]]);
        }
        return r;
    }
}
