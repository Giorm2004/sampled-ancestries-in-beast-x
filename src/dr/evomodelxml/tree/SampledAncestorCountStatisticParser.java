/*
 * SampledAncestorCountStatisticParser.java
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

import dr.evolution.tree.Tree;
import dr.evomodel.tree.SampledAncestorCountStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/** XML parser for SampledAncestorCountStatistic. */
public class SampledAncestorCountStatisticParser extends AbstractXMLObjectParser {

    public String getParserName() {
        return "sampledAncestorCount";
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        return new SampledAncestorCountStatistic(xo.getAttribute(Statistic.NAME, xo.getId()), (Tree) xo.getChild(Tree.class));
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{AttributeRule.newStringRule(Statistic.NAME, true), new ElementRule(Tree.class)};
    }

    public String getParserDescription() {
        return "The number of sampled ancestors (tips on zero-length branches) in a tree";
    }

    public Class getReturnType() {
        return SampledAncestorCountStatistic.class;
    }
}
