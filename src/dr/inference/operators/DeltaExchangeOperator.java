/*
 * DeltaExchangeOperator.java
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

package dr.inference.operators;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

/**
 * A generic operator for use with a sum-constrained (possibly weighted) vector parameter.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class DeltaExchangeOperator extends AbstractAdaptableOperator {

    public DeltaExchangeOperator(Parameter parameter, double delta) {

        super(AdaptationMode.ADAPTATION_ON);

        this.parameter = parameter;
        this.delta = delta;
        setWeight(1.0);
        this.isIntegerOperator = false;

        parameterWeights = new int[parameter.getDimension()];
        for (int i = 0; i < parameterWeights.length; i++) {
            parameterWeights[i] = 1;
        }
        this.stepScale = stepScales(parameterWeights);
    }

    public DeltaExchangeOperator(Parameter parameter, int[] parameterWeights, double delta, double weight, boolean isIntegerOperator, AdaptationMode mode) {

        super(mode);

        this.parameter = parameter;
        this.delta = delta;
        setWeight(weight);
        this.isIntegerOperator = isIntegerOperator;
        this.parameterWeights = parameterWeights;
        this.stepScale = stepScales(parameterWeights);

        if (isIntegerOperator && delta != Math.round(delta)) {
            throw new IllegalArgumentException("Can't be an integer operator if delta is not integer");
        }
    }

    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return parameter;
    }

    /**
     * change the parameter and return the hastings ratio.
     * performs a delta exchange operation between two scalars in the vector
     * and return the hastings ratio.
     */
    public final double doOperation() {

        // get two dimensions
        final int dim = parameter.getDimension();
        final int dim1 = MathUtils.nextInt(dim);
        int dim2 = dim1;
        while (dim1 == dim2) {
            dim2 = MathUtils.nextInt(dim);
        }

        double scalar1 = parameter.getParameterValue(dim1);
        double scalar2 = parameter.getParameterValue(dim2);

        if (isIntegerOperator) {
            int d = MathUtils.nextInt((int) Math.round(delta)) + 1;

            if (parameterWeights[dim1] != parameterWeights[dim2]) throw new RuntimeException();
            scalar1 = Math.round(scalar1 - d);
            scalar2 = Math.round(scalar2 + d);
        } else {

            // exchange a random delta of the weighted sum: x1 -= d * wbar / w1, x2 += d * wbar / w2. The reverse move
            // (dim2, dim1) uses the same d, so the proposal is symmetric for any weights. (The previous form,
            // x1 -= d, x2 += d * w1 / w2, is not: its reverse needs d * w1 / w2 <= delta and has a different density,
            // so with unequal weights it sampled a delta-dependent distortion of the target.) With equal weights the
            // step scales are 1 and the move is unchanged.
            final double d = MathUtils.nextDouble() * delta;
            scalar1 -= d * stepScale[dim1];
            scalar2 += d * stepScale[dim2];

        }
        Bounds<Double> bounds = parameter.getBounds();

        if (scalar1 < bounds.getLowerLimit(dim1) ||
                scalar1 > bounds.getUpperLimit(dim1) ||
                scalar2 < bounds.getLowerLimit(dim2) ||
                scalar2 > bounds.getUpperLimit(dim2)) {
            return Double.NEGATIVE_INFINITY;
        }
        parameter.setParameterValue(dim1, scalar1);
        parameter.setParameterValue(dim2, scalar2);

        // symmetrical move so return a zero hasting ratio
        return 0.0;
    }

    // Interface MCMCOperator
    public final String getOperatorName() {
        return "deltaExchange(" + parameter.getParameterName() + ")";
    }

    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(delta);
    }

    public void setAdaptableParameterValue(double value) {
        delta = Math.exp(value);
    }

    public double getRawParameter() {
        return delta;
    }

    public String getAdaptableParameterName() {
        return "delta";
    }

    public String toString() {
        return getOperatorName() + "(windowsize=" + delta + ")";
    }

    /** wbar / w_k, so that a step of d in the weighted sum moves each component by d * wbar / w_k. */
    private static double[] stepScales(int[] weights) {
        double mean = 0.0;
        for (int w : weights) {
            mean += w;
        }
        mean /= weights.length;
        double[] scale = new double[weights.length];
        for (int k = 0; k < weights.length; k++) {
            scale[k] = weights[k] == mean ? 1.0 : mean / weights[k];
        }
        return scale;
    }

    // Private instance variables

    private Parameter parameter = null;
    private final int[] parameterWeights;
    private final double[] stepScale;
    private double delta = 0.02;
    private boolean isIntegerOperator = false;
}
