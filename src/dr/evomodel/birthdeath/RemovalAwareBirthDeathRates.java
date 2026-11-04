/*
 * RemovalAwareBirthDeathRates.java
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

package dr.evomodel.birthdeath;

import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.VariableListener;

/**
 * Birth, death and serial-sampling rates from the epidemiological parameterisation of the birth-death
 * skyline with a removal probability r, as in BEAST 2's BDSKY (BirthDeathSkylineModel.transformParameters,
 * sampled-ancestor case):
 *
 *     delta = mu + r psi      (become-uninfectious rate)
 *     s     = psi / (mu + psi)  (sampling proportion)
 *     R0    = lambda / delta
 *
 * so  lambda = R0 delta,  psi = s delta / (1 + (r - 1) s),  mu = delta - r psi.
 * With r = 1 this is CompoundBirthDeathParameters' transform (psi = D S, mu = D (1 - S)); with r < 1 the
 * two differ, so priors written for BDSKY carry over only through this class.
 *
 * Each input may have dimension 1 (shared by every interval) or the number of intervals; the three
 * rates are read-only proxies with one value per interval.
 *
 * @author Giovanni Mazzeo
 */
public class RemovalAwareBirthDeathRates {

    public final Parameter birthRate;
    public final Parameter deathRate;
    public final Parameter samplingRate;

    private final Parameter reproductiveNumber;
    private final Parameter becomeUninfectiousRate;
    private final Parameter samplingProportion;
    private final Parameter removalProbability;

    public RemovalAwareBirthDeathRates(Parameter reproductiveNumber, Parameter becomeUninfectiousRate,
                                       Parameter samplingProportion, Parameter removalProbability,
                                       int numIntervals) {
        for (Parameter p : new Parameter[]{reproductiveNumber, becomeUninfectiousRate, samplingProportion,
                removalProbability}) {
            if (p.getDimension() != 1 && p.getDimension() != numIntervals) {
                throw new IllegalArgumentException("Parameter " + p.getId() + " has dimension " + p.getDimension()
                        + "; expected 1 or the number of intervals (" + numIntervals + ")");
            }
        }
        this.reproductiveNumber = reproductiveNumber;
        this.becomeUninfectiousRate = becomeUninfectiousRate;
        this.samplingProportion = samplingProportion;
        this.removalProbability = removalProbability;
        reproductiveNumber.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, reproductiveNumber.getDimension()));
        becomeUninfectiousRate.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, becomeUninfectiousRate.getDimension()));
        samplingProportion.addBounds(new Parameter.DefaultBounds(1.0, 0.0, samplingProportion.getDimension()));

        this.birthRate = new Rate("lambda", numIntervals) {
            double value(int i) {
                return at(reproductiveNumber, i) * at(becomeUninfectiousRate, i);
            }
        };
        this.samplingRate = new Rate("psi", numIntervals) {
            double value(int i) {
                return psi(i);
            }
        };
        this.deathRate = new Rate("mu", numIntervals) {
            double value(int i) {
                return at(becomeUninfectiousRate, i) - at(removalProbability, i) * psi(i);
            }
        };
    }

    private double psi(int i) {
        double s = at(samplingProportion, i);
        return s * at(becomeUninfectiousRate, i) / (1.0 + (at(removalProbability, i) - 1.0) * s);
    }

    private static double at(Parameter p, int i) {
        return p.getDimension() == 1 ? p.getParameterValue(0) : p.getParameterValue(i);
    }

    private abstract class Rate extends Parameter.Proxy implements VariableListener {

        Rate(String name, int dimension) {
            super(name, dimension);
            reproductiveNumber.addVariableListener(this);
            becomeUninfectiousRate.addVariableListener(this);
            samplingProportion.addVariableListener(this);
            removalProbability.addVariableListener(this);
        }

        abstract double value(int i);

        @Override
        public double getParameterValue(int i) {
            return value(i);
        }

        @Override
        public void setParameterValue(int i, double value) {
            throw new UnsupportedOperationException("derived rate; set R0, becomeUninfectiousRate or samplingProportion");
        }

        @Override
        public void setParameterValueQuietly(int i, double value) {
            throw new UnsupportedOperationException("derived rate");
        }

        @Override
        public void setParameterValueNotifyChangedAll(int i, double value) {
            throw new UnsupportedOperationException("derived rate");
        }

        @Override
        public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            fireParameterChangedEvent();
        }
    }
}
