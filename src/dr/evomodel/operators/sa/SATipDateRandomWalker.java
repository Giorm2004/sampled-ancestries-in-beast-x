/*
 * SATipDateRandomWalker.java
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

import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.SampledAncestorTrees;
import dr.evomodel.tree.TreeModel;
import dr.math.MathUtils;

/**
 * Random walk on the height of a tip drawn from a taxon set; a sampled ancestor moves together with its fake parent
 * (between its fake node's other child and the node above), an ordinary tip between 0 and its parent. Tips at
 * height 0 are not moved. Port of BEAST 2 SA's SampledNodeDateRandomWalker (Gaussian or uniform window); Hastings
 * ratio 1.
 *
 * @author Giovanni Mazzeo
 */
public class SATipDateRandomWalker extends SATreeOperator {

    private final int[] tips;
    private final double windowSize;
    private final boolean gaussian;

    public SATipDateRandomWalker(TreeModel tree, TaxonList taxa, double windowSize, boolean gaussian, double weight) {
        super(tree, null, weight);
        this.windowSize = windowSize;
        this.gaussian = gaussian;
        this.tips = new int[taxa.getTaxonCount()];
        for (int k = 0; k < tips.length; k++) {
            Taxon taxon = taxa.getTaxon(k);
            int index = tree.getTaxonIndex(taxon);
            if (index < 0) {
                throw new IllegalArgumentException("taxon " + taxon.getId() + " is not in tree " + tree.getId());
            }
            tips[k] = index;
        }
    }

    @Override
    public double doOperation() {
        NodeRef node = tree.getExternalNode(tips[MathUtils.nextInt(tips.length)]);
        final double value = h(node);
        if (value == 0.0) {
            return Double.NEGATIVE_INFINITY;
        }
        final double newValue = value + (gaussian ? MathUtils.nextGaussian() * windowSize
                : (2.0 * MathUtils.nextDouble() - 1.0) * windowSize);
        NodeRef fake = null;
        double lower, upper;
        if (isSA(node)) {
            fake = tree.getParent(node);
            lower = h(SampledAncestorTrees.otherChild(tree, fake, node));
            upper = tree.isRoot(fake) ? Double.POSITIVE_INFINITY : h(tree.getParent(fake));
        } else {
            lower = 0.0;
            upper = h(tree.getParent(node));
        }
        if (newValue < lower || newValue > upper || newValue == value) {
            return Double.NEGATIVE_INFINITY;
        }
        if (fake != null) {
            tree.setNodeHeight(fake, newValue);
        }
        tree.setNodeHeight(node, newValue);
        return 0.0;
    }

    @Override
    public String getOperatorName() {
        return "saTipDateRandomWalk(" + tree.getId() + ")";
    }
}
