/*
 * TMRCAPriorGradient.java
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

package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightProxyParameter;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.util.Attribute;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.List;

/**
 * Gradient, w.r.t. the internal node heights (root included, node-number order), of a set of priors on clade ages:
 * each prior is a DistributionLikelihood whose single datum is a TMRCAStatistic (node height of the MRCA, not an
 * absolute time and not the parent). The prior's derivative at the MRCA node's height is added to that node.
 * <p>
 * The derivative is a central difference of the prior's own log density (DistributionLikelihood.logDensityAt, so
 * offsets, scales and any distribution are handled the same way as in the likelihood), one-sided next to a support
 * bound. It is a deterministic function of the heights, which is all HMC needs to stay exact; it is accurate to
 * ~1e-8 relative for smooth densities.
 *
 * @author Giovanni Mazzeo
 */
public class TMRCAPriorGradient implements GradientWrtParameterProvider, Reportable {

    private final TreeModel tree;
    private final List<DistributionLikelihood> priors;
    private final List<TMRCAStatistic> statistics;
    private final NodeHeightProxyParameter parameter;
    private final Likelihood likelihood;

    public TMRCAPriorGradient(String name, TreeModel tree, List<DistributionLikelihood> priors) {
        this.tree = tree;
        this.priors = priors;
        this.statistics = new ArrayList<>();
        for (DistributionLikelihood prior : priors) {
            List<Attribute<double[]>> data = prior.getDataList();
            if (data.size() != 1 || !(data.get(0) instanceof TMRCAStatistic)) {
                throw new IllegalArgumentException("each prior must have exactly one tmrcaStatistic as its data");
            }
            TMRCAStatistic statistic = (TMRCAStatistic) data.get(0);
            if (statistic.getTree() != tree) {
                throw new IllegalArgumentException("tmrcaStatistic " + statistic.getStatisticName() + " is on another tree");
            }
            if (statistic.getLeafSet() == null) {
                throw new IllegalArgumentException("tmrcaStatistic " + statistic.getStatisticName() + " has no taxa");
            }
            statistics.add(statistic);
        }
        this.parameter = new NodeHeightProxyParameter(name + ".heights", tree, true);
        this.likelihood = new CompoundLikelihood(new ArrayList<Likelihood>(priors));
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] gradient = new double[getDimension()];
        int external = tree.getExternalNodeCount();
        for (int k = 0; k < priors.size(); k++) {
            TMRCAStatistic statistic = statistics.get(k);
            NodeRef node = TreeUtils.getCommonAncestorNode(tree, statistic.getLeafSet());
            double height = tree.getNodeHeight(node);
            if (Math.abs(statistic.getStatisticValue(0) - height) > 1e-9 * Math.max(1.0, height)) {
                throw new IllegalStateException("tmrcaStatistic " + statistic.getStatisticName()
                        + " is not the MRCA node height (absolute or forParent statistics are not supported)");
            }
            gradient[node.getNumber() - external] += derivative(priors.get(k), height);
        }
        return gradient;
    }

    static double derivative(DistributionLikelihood prior, double x) {
        double h = 1e-6 * Math.max(1.0, Math.abs(x));
        double up = prior.logDensityAt(x + h);
        double down = prior.logDensityAt(x - h);
        if (isFinite(up) && isFinite(down)) {
            return (up - down) / (2.0 * h);
        }
        double here = prior.logDensityAt(x);
        if (!isFinite(here)) {
            return 0.0;                 // outside the support: the joint density rejects the state anyway
        }
        if (isFinite(up)) {
            return (up - here) / h;
        }
        if (isFinite(down)) {
            return (here - down) / h;
        }
        return 0.0;
    }

    private static boolean isFinite(double x) {
        return !Double.isNaN(x) && !Double.isInfinite(x);
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, 1e-3);
    }
}
