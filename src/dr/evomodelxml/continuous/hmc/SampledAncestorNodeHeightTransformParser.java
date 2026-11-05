/*
 * SampledAncestorNodeHeightTransformParser.java
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

package dr.evomodelxml.continuous.hmc;

import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.discrete.SampledAncestorNodeHeightTransform;
import dr.xml.*;

/**
 * XML parser for SampledAncestorNodeHeightTransform:
 *
 *     &lt;saNodeHeightTransform id="..."&gt; &lt;treeModel idref="..."/&gt; &lt;/saNodeHeightTransform&gt;
 *
 * Use with a hamiltonianMonteCarloOperator on a nodeHeightProxyParameter with rootNode="true".
 */
public class SampledAncestorNodeHeightTransformParser extends AbstractXMLObjectParser {

    public String getParserName() {
        return "saNodeHeightTransform";
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        return new SampledAncestorNodeHeightTransform((TreeModel) xo.getChild(TreeModel.class));
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{new ElementRule(TreeModel.class)};
    }

    public String getParserDescription() {
        return "Real-line node-height transform for sampled-ancestor trees: log/logit ratios, fake nodes pinned";
    }

    public Class getReturnType() {
        return SampledAncestorNodeHeightTransform.class;
    }
}
