/*
 * TMRCAPriorGradientParser.java
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

package dr.evomodelxml.tree;

import dr.evomodel.tree.TMRCAPriorGradient;
import dr.evomodel.tree.TreeModel;
import dr.inference.distribution.DistributionLikelihood;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * {@code <tmrcaPriorGradient id="..."><treeModel idref="..."/><normalPrior idref="..."/>...</tmrcaPriorGradient>}:
 * node-height gradient of priors on tmrcaStatistics (clade-age calibrations), for HMC on node heights.
 */
public class TMRCAPriorGradientParser extends AbstractXMLObjectParser {

    public static final String NAME = "tmrcaPriorGradient";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        List<DistributionLikelihood> priors = new ArrayList<>();
        for (int i = 0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);
            if (child instanceof DistributionLikelihood) {
                priors.add((DistributionLikelihood) child);
            }
        }
        try {
            return new TMRCAPriorGradient(xo.getId(), tree, priors);
        } catch (IllegalArgumentException e) {
            throw new XMLParseException(e.getMessage());
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(TreeModel.class),
                new ElementRule(DistributionLikelihood.class, 1, Integer.MAX_VALUE),
        };
    }

    @Override
    public String getParserDescription() {
        return "Gradient w.r.t. the internal node heights of priors on clade ages (tmrcaStatistics)";
    }

    @Override
    public Class getReturnType() {
        return TMRCAPriorGradient.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
