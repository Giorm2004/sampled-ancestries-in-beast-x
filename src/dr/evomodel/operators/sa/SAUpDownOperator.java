/*
 * SAUpDownOperator.java
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

package dr.evomodel.operators.sa;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter;

import java.util.List;

/**
 * Scales the non-fake internal node heights of a tree (and any "up" parameters) by s and the "down" parameters by
 * 1/s, s uniform on [f, 1/f]. Port of BEAST 2's UpDownOperator for sampled-ancestor trees (BEAST 2's Tree.scale
 * leaves fake nodes on their sampled ancestors): log Hastings ratio (up - down - 2) log(s), counting scaled
 * dimensions (the tree counts its non-fake internal nodes).
 *
 * @author Giovanni Mazzeo
 */
public class SAUpDownOperator extends SATreeOperator {

    private final double scaleFactor;
    private final List<Parameter> up;
    private final List<Parameter> down;

    public SAUpDownOperator(TreeModel tree, List<Parameter> up, List<Parameter> down, double scaleFactor, double weight) {
        super(tree, null, weight);
        this.up = up;
        this.down = down;
        this.scaleFactor = scaleFactor;
    }

    @Override
    public double doOperation() {
        final double scale = SAScaleOperator.drawScale(scaleFactor);
        int goingUp = 0, goingDown = 0;
        for (Parameter p : up) {
            if (!inBounds(p, scale)) {
                return Double.NEGATIVE_INFINITY;
            }
        }
        for (Parameter p : down) {
            if (!inBounds(p, 1.0 / scale)) {
                return Double.NEGATIVE_INFINITY;
            }
        }
        if (tree != null) {
            int dof = SAScaleOperator.scaleTree(tree, scale);
            if (dof < 0) {
                return Double.NEGATIVE_INFINITY;
            }
            goingUp += dof;
        }
        for (Parameter p : up) {
            goingUp += scale(p, scale);
        }
        for (Parameter p : down) {
            goingDown += scale(p, 1.0 / scale);
        }
        return (goingUp - goingDown - 2) * Math.log(scale);
    }

    private static boolean inBounds(Parameter p, double scale) {
        Bounds<Double> bounds = p.getBounds();
        for (int i = 0; i < p.getDimension(); i++) {
            double v = p.getParameterValue(i) * scale;
            if (v < bounds.getLowerLimit(i) || v > bounds.getUpperLimit(i)) {
                return false;
            }
        }
        return true;
    }

    private static int scale(Parameter p, double scale) {
        for (int i = 0; i < p.getDimension(); i++) {
            p.setParameterValue(i, p.getParameterValue(i) * scale);
        }
        return p.getDimension();
    }

    @Override
    public String getOperatorName() {
        return "saUpDown(" + (tree != null ? tree.getId() : "") + ")";
    }
}
